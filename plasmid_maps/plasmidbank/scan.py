"""Scan one map folder into the SQLite index, mirror the readable maps, write a hygiene report.

    python -m plasmidbank.scan <root> [--db PATH] [--mirror DIR] [--report PATH]

Read-only on <root>: the index, mirror and report must all live outside it (checked). Walks only
<root>; skips dot files/folders, CLC workspace clutter (bank.SKIP_SUFFIXES, `segments*`) and files
over bank.MAX_BYTES before reading them. Unchanged files (same size + mtime) are not re-read.
A killed run is simply re-run: everything committed so far is reused from the cache.

Exit 0 on success (warnings and per-file errors included), 2 on a fatal error.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

from .bank import MAX_BYTES, SKIP_SUFFIXES
from .formats import MapRecord, read_map
from .ids import identify
from .index import DEFAULT_DB, DEFAULT_HOME, Index, nfc, now
from .mirror import DEFAULT_MIRROR, MIRRORED_FORMATS, mirror_bytes

COMMIT_EVERY = 100
MAX_LISTED = 200  # skipped-large paths and walk errors kept in the run row
# Folders in the bank that hold no plasmid maps and are slow to read over SMB (for
# example a genome folder with very large CLC files). Matched on the name with spaces trimmed, at
# any depth; their files are counted as skipped, and a complete run marks earlier rows missing.
# Workspaces/: My workspace files mirrored from the app (uploads.sync_workspaces, 2026-10-27), not bank maps.
# _replaced/: files overwritten from the app, kept for undo, not maps of the bank any more.
EXCLUDE_DIRS = {"Morex genome v1.0 data files", "Workspaces", "_replaced"}


class FatalError(Exception):
    pass


def _inside(child: Path, parent: Path) -> bool:
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def skip_reason(name: str) -> str | None:
    low = name.lower()
    if name.startswith("."):
        return "dotfile"
    suffix = os.path.splitext(low)[1]
    if suffix in SKIP_SUFFIXES:
        return suffix
    if low.startswith("segments"):  # Lucene segment files in CLC workspaces
        return "segments"
    return None


def _error_record(path: Path, msg: str) -> MapRecord:
    ids = identify(path.name)
    return MapRecord(
        str(path),
        ids.plasmid_id,
        ids.kind,
        ids.source,
        ids.xrefs,
        "unknown",
        "error",
        [msg],
    )


def scan(
    root: str | Path,
    db: str | Path = DEFAULT_DB,
    mirror_dir: str | Path = DEFAULT_MIRROR,
    log=print,
) -> dict:
    root = Path(root).expanduser().resolve()
    mirror_dir = Path(mirror_dir).expanduser()
    db = Path(db).expanduser()
    if not root.is_dir():
        raise FatalError(f"not a folder: {root}")
    for what, p in (("index", db), ("mirror", mirror_dir)):
        if _inside(p, root):
            raise FatalError(
                f"refusing to write the {what} inside the scanned folder: {p}"
            )

    idx = Index(db)
    run_id, started = idx.start_run(str(root))
    c: Counter = Counter()
    skipped: Counter = Counter()
    large: list[str] = []
    errors: list[str] = []

    def onerror(e: OSError) -> None:
        errors.append(f"walk: {e.filename}: {e.strerror or e}")

    try:
        for dirpath, dirnames, filenames in os.walk(root, onerror=onerror):
            excluded = [d for d in dirnames if d.strip() in EXCLUDE_DIRS]
            if excluded:
                skipped["excluded folder"] += len(excluded)
            dirnames[:] = sorted(
                d
                for d in dirnames
                if not d.startswith(".") and d.strip() not in EXCLUDE_DIRS
            )
            for fn in sorted(filenames):
                reason = skip_reason(fn)
                if reason:
                    skipped[reason] += 1
                    continue
                p = Path(dirpath) / fn
                try:
                    st = p.lstat()
                except OSError as e:
                    errors.append(f"stat: {p}: {e}")
                    continue
                if not p.is_file() or p.is_symlink():
                    skipped["not a regular file"] += 1
                    continue
                if st.st_size > MAX_BYTES:
                    skipped["over MAX_BYTES"] += 1
                    if len(large) < MAX_LISTED:
                        large.append(nfc(p))
                    continue
                c["seen"] += 1
                _one(idx, p, st, started, mirror_dir, c, errors)
                if c["seen"] % COMMIT_EVERY == 0:
                    idx.commit()
                    log(
                        f"  {c['seen']} files ({c['parsed']} parsed, {c['cached']} cached)"
                    )
    except KeyboardInterrupt:
        idx.commit()
        raise
    complete = not any(e.startswith("walk:") for e in errors)
    counts = {
        **dict(c),
        "skipped": dict(skipped),
        "skipped_large": large,
        "complete": complete,
    }
    for k in ("seen", "parsed", "cached", "mirrored_new", "read_errors"):
        counts.setdefault(k, 0)
    missing = idx.finish_run(run_id, counts, errors[:MAX_LISTED], mark_missing=complete)
    counts["newly_missing"] = missing
    counts["run_id"] = run_id
    counts["db"] = str(db)
    counts["root"] = str(root)
    idx.close()
    return counts


def _one(
    idx: Index, p: Path, st, seen: str, mirror_dir: Path, c: Counter, errors
) -> None:
    row = idx.cached(str(p), st.st_size, st.st_mtime_ns)
    if row is not None:
        idx.touch(row["id"], seen)
        c["cached"] += 1
        mp = row["mirror_path"]
        needs_mirror = row["format"] in MIRRORED_FORMATS and not (
            mp and Path(mp).exists() and _inside(Path(mp), mirror_dir)
        )
        if not needs_mirror:
            return
        try:
            data = p.read_bytes()
        except OSError as e:
            errors.append(f"read: {p}: {e}")
            return
        if hashlib.sha256(data).hexdigest() == row["sha256"]:
            dest, new = mirror_bytes(
                data,
                kind=row["kind"],
                plasmid_id=row["plasmid_id"],
                sha256=row["sha256"],
                fmt=row["format"],
                mirror_dir=mirror_dir,
            )
            idx.set_mirror(row["id"], str(dest))
            c["mirrored_new"] += new
            return
        c["cached"] -= (
            1  # content changed under an unchanged size+mtime: re-parse below
        )

    c["parsed"] += 1
    try:
        data = p.read_bytes()
    except OSError as e:
        c["read_errors"] += 1
        errors.append(f"read: {p}: {e}")
        idx.upsert(
            _error_record(p, f"read error: {e}"),
            size=st.st_size,
            mtime_ns=None,
            sha256=None,
            seen=seen,
        )
        return
    sha = hashlib.sha256(data).hexdigest()
    try:
        rec = read_map(p)
    except Exception as e:  # one bad file must never stop the scan
        rec = _error_record(p, f"read_map failed: {type(e).__name__}: {e}")
    mtime = st.st_mtime_ns
    try:
        after = p.stat()
        if (after.st_size, after.st_mtime_ns) != (st.st_size, st.st_mtime_ns):
            rec.messages.append("file changed during the scan; re-read next run")
            mtime = None
    except OSError:
        mtime = None
    file_id = idx.upsert(rec, size=st.st_size, mtime_ns=mtime, sha256=sha, seen=seen)
    dest, new = mirror_bytes(
        data,
        kind=rec.kind,
        plasmid_id=rec.plasmid_id,
        sha256=sha,
        fmt=rec.format,
        mirror_dir=mirror_dir,
    )
    idx.set_mirror(file_id, str(dest) if dest else None)
    c["mirrored_new"] += new


# ── report ────────────────────────────────────────────────────────────────────

CSS = """
:root{--bg:#fbfaf7;--fg:#1d1f21;--muted:#5f6368;--card:#ffffff;--line:#e3e1db;--accent:#1f6f5c;
--warn:#9a6700;--bad:#b3261e;--ok:#1f6f5c;--chip:#eef3f1;color-scheme:light dark}
@media (prefers-color-scheme: dark){:root{--bg:#151718;--fg:#e8e6e3;--muted:#9aa0a6;--card:#1d2022;
--line:#2f3336;--accent:#5fc1a5;--warn:#e3b341;--bad:#f28b82;--ok:#5fc1a5;--chip:#22302c}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",
Roboto,Helvetica,Arial,sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:1.6rem;margin:0 0 4px}h2{font-size:1.15rem;margin:32px 0 8px;border-bottom:1px solid var(--line);
padding-bottom:4px}.muted{color:var(--muted)}
.tiles{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:10px;margin:16px 0}
.tile{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px 12px}
.tile b{display:block;font-size:1.4rem;font-variant-numeric:tabular-nums}
.tile span{color:var(--muted);font-size:.85rem}
.tile.warn b{color:var(--warn)}.tile.bad b{color:var(--bad)}
.grid3{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:16px}
table{border-collapse:collapse;width:100%;background:var(--card);border:1px solid var(--line);font-size:.88rem}
th,td{text-align:left;padding:5px 8px;border-bottom:1px solid var(--line);vertical-align:top}
th{font-weight:600;background:var(--chip)}td.n{text-align:right;font-variant-numeric:tabular-nums}
.wrap{overflow-x:auto}code,.path{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
font-size:.82rem;word-break:break-all}
.ids{background:var(--card);border:1px solid var(--line);padding:8px 10px;font-family:ui-monospace,Menlo,monospace;
font-size:.8rem;line-height:1.7;word-spacing:6px}
.note{font-size:.88rem;color:var(--muted);margin:4px 0 8px}
.chip{display:inline-block;background:var(--chip);border-radius:4px;padding:0 6px;margin-right:4px;font-size:.8rem}
.empty{color:var(--ok)}
@media print{body{background:#fff;color:#000;font-size:10pt}main{max-width:none;padding:0}
.tile,table,.ids{border-color:#999;background:#fff}th{background:#eee}h2{break-after:avoid}
tr{break-inside:avoid}a{color:#000}}
"""

FAVICON = (
    "data:image/svg+xml,"
    "%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E"
    "%3Ccircle cx='16' cy='16' r='11' fill='none' stroke='%231f6f5c' stroke-width='4'/%3E"
    "%3Crect x='13' y='2' width='6' height='8' rx='2' fill='%23e3b341'/%3E%3C/svg%3E"
)

MAX_ROWS = 2000


def _e(s) -> str:
    return html.escape("" if s is None else str(s))


def _rel(path: str, root: str) -> str:
    prefix = root.rstrip("/") + "/"
    return path[len(prefix) :] if path.startswith(prefix) else path


def _table(
    headers: list[str], rows: list[list[str]], numeric: set[int] = frozenset()
) -> str:
    if not rows:
        return '<p class="empty">None.</p>'
    more = ""
    if len(rows) > MAX_ROWS:
        more = f'<p class="note">Showing {MAX_ROWS} of {len(rows)}.</p>'
        rows = rows[:MAX_ROWS]
    head = "".join(f"<th>{_e(h)}</th>" for h in headers)
    body = "".join(
        "<tr>"
        + "".join(
            f'<td class="n">{cell}</td>' if i in numeric else f"<td>{cell}</td>"
            for i, cell in enumerate(r)
        )
        + "</tr>"
        for r in rows
    )
    return f'<div class="wrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>{more}'


def _counts_table(title: str, d: dict) -> str:
    rows = [[_e(k), str(v)] for k, v in d.items()]
    return f"<div><h3>{_e(title)}</h3>{_table([title, 'files'], rows, {1})}</div>"


def _msgs(row) -> str:
    try:
        return "; ".join(json.loads(row["messages"] or "[]"))
    except ValueError:
        return row["messages"] or ""


def render_report(db: str | Path, root: str | Path) -> str:
    idx = Index(db)
    root = nfc(Path(root).expanduser().resolve())
    run = idx.last_run(root)
    counts = json.loads(run["counts"]) if run and run["counts"] else {}
    run_errors = json.loads(run["errors"]) if run and run["errors"] else []
    files = idx.files(root)
    unident = idx.unidentified(root)
    locus = idx.from_locus(root)
    amb = idx.ambiguous_ids(root)
    dups = idx.duplicates(root)
    perr = idx.parse_errors(root)
    oob = idx.out_of_bin(root)
    clc_only = idx.clc_only(root)
    missing = idx.missing(root)
    r = lambda p: f'<span class="path">{_e(_rel(p, root))}</span>'  # noqa: E731

    def tile(n, label, cls=""):
        return f'<div class="tile {cls}"><b>{n:,}</b><span>{_e(label)}</span></div>'

    tiles = "".join(
        [
            tile(len(files), "files indexed"),
            tile(
                len({f["plasmid_id"] for f in files if f["plasmid_id"]}), "distinct IDs"
            ),
            tile(counts.get("parsed", 0), "parsed this run"),
            tile(counts.get("cached", 0), "reused from cache"),
            tile(counts.get("mirrored_new", 0), "new mirror copies"),
            tile(len(unident), "unidentified", "warn" if unident else ""),
            tile(len(amb), "ambiguous IDs", "warn" if amb else ""),
            tile(len(dups), "duplicate groups", "warn" if dups else ""),
            tile(len(perr), "parse errors", "bad" if perr else ""),
            tile(len(oob), "out of bin", "warn" if oob else ""),
            tile(len(clc_only), "CLC-only IDs", "warn" if clc_only else ""),
            tile(len(missing), "missing files", "warn" if missing else ""),
        ]
    )

    amb_rows = [
        [
            _e(pid),
            str(len(rows)),
            _e(", ".join(sorted({x["format"] for x in rows}))),
            "<br>".join(r(x["path"]) for x in rows),
        ]
        for pid, rows in amb.items()
    ]
    text_amb = sum(
        1 for rows in amb.values() if sum(x["format"] != "clc" for x in rows) > 1
    )
    dup_rows = [
        [
            str(i + 1),
            _e(g[0]["length"]),
            "<br>".join(f"{_e(x['plasmid_id'] or '—')} · {r(x['path'])}" for x in g),
        ]
        for i, g in enumerate(dups)
    ]
    skipped = counts.get("skipped", {})
    started = run["started"] if run else "never"
    finished = run["finished"] if run and run["finished"] else "not finished"
    complete = counts.get("complete")
    status_line = (
        (
            "Walk complete; files not seen were marked missing."
            if complete
            else "Walk INCOMPLETE (see walk errors): nothing was marked missing this run."
        )
        if run
        else ""
    )

    sections = [
        f"<h2>Counts</h2><div class='grid3'>{_counts_table('status', idx.counts('status', root))}"
        f"{_counts_table('format', idx.counts('format', root))}"
        f"{_counts_table('kind', idx.counts('kind', root))}</div>",
        "<h2>Unidentified</h2><p class='note'>No PMC/EC number in the filename or the LOCUS line.</p>"
        + _table(
            ["file", "format", "status"],
            [[r(x["path"]), _e(x["format"]), _e(x["status"])] for x in unident],
        ),
        "<h2>ID from LOCUS line only</h2><p class='note'>The filename has no ID; the LOCUS line supplied it. "
        "Check, then rename the file.</p>"
        + _table(["ID", "file"], [[_e(x["plasmid_id"]), r(x["path"])] for x in locus]),
        f"<h2>Ambiguous IDs</h2><p class='note'>More than one file holds the ID, so it is never synced. "
        f"{text_amb} of {len(amb)} have more than one readable (non-CLC) file; the rest are "
        f".clc + text pairs.</p>"
        + _table(["ID", "files", "formats", "paths"], amb_rows, {1}),
        "<h2>Duplicate sequences</h2><p class='note'>Identical sequence, regardless of rotation and "
        "strand, under more than one file.</p>"
        + _table(["group", "length", "ID · file"], dup_rows, {0, 1}),
        "<h2>Parse errors</h2>"
        + _table(
            ["file", "format", "messages"],
            [[r(x["path"]), _e(x["format"]), _e(_msgs(x))] for x in perr],
        ),
        "<h2>Out of bin</h2><p class='note'>PMC number outside the range of the "
        "<code>PMC-xxxxx - PMC-yyyyy</code> folder it sits in (folder names matched loosely).</p>"
        + _table(
            ["ID", "bin folder", "range", "file"],
            [
                [_e(x["plasmid_id"]), _e(f), f"{lo}–{hi}", r(x["path"])]
                for x, f, lo, hi in oob
            ],
        ),
        "<h2>CLC-only IDs</h2><p class='note'>Every file for these IDs is a CLC binary: nothing readable "
        "to mirror or sync until they are exported to GenBank.</p>"
        + (
            f"<div class='ids'>{' '.join(_e(i) for i in clc_only)}</div>"
            if clc_only
            else '<p class="empty">None.</p>'
        ),
        "<h2>Missing files</h2><p class='note'>Indexed earlier, not found since. Rows are kept.</p>"
        + _table(
            ["file", "ID", "missing since"],
            [
                [r(x["path"]), _e(x["plasmid_id"]), _e(x["missing_since"])]
                for x in missing
            ],
        ),
        "<h2>Skipped</h2>"
        + _table(
            ["reason", "files"],
            [[_e(k), str(v)] for k, v in sorted(skipped.items())],
            {1},
        )
        + (
            (
                "<h3>Over the size limit</h3>"
                + _table(["file"], [[r(p)] for p in counts.get("skipped_large", [])])
            )
            if counts.get("skipped_large")
            else ""
        ),
        "<h2>Walk and read errors</h2>"
        + _table(["error"], [[f"<code>{_e(e)}</code>"] for e in run_errors]),
    ]
    idx.close()
    title = "Plasmid bank scan"
    return (
        '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<meta name="color-scheme" content="light dark">\n'
        f'<title>{title}</title>\n<link rel="icon" href="{FAVICON}">\n'
        f"<style>{CSS}</style>\n</head>\n<body>\n<main>\n"
        f"<h1>{title}</h1>\n<p class='muted'>Folder <code>{_e(root)}</code><br>"
        f"Last run {_e(started)} → {_e(finished)} · report written "
        f"{datetime.now():%Y-%m-%d %H:%M} · index <code>{_e(db)}</code></p>\n"
        f"<p>{_e(status_line)}</p>\n<div class='tiles'>{tiles}</div>\n"
        + "\n".join(f"<section>{s}</section>" for s in sections)
        + "\n</main>\n</body>\n</html>\n"
    )


def write_report(db: str | Path, root: str | Path, out: str | Path) -> Path:
    out = Path(out).expanduser()
    if _inside(out, Path(root).expanduser()):
        raise FatalError(
            f"refusing to write the report inside the scanned folder: {out}"
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f".{out.name}.tmp")
    tmp.write_text(render_report(db, root), encoding="utf-8")
    os.replace(tmp, out)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m plasmidbank.scan", description=__doc__.split("\n")[0]
    )
    ap.add_argument("root", help="the map folder to scan (read-only)")
    ap.add_argument(
        "--db", default=str(DEFAULT_DB), help=f"index (default {DEFAULT_DB})"
    )
    ap.add_argument(
        "--mirror",
        default=str(DEFAULT_MIRROR),
        help=f"mirror folder (default {DEFAULT_MIRROR})",
    )
    ap.add_argument(
        "--report",
        default=None,
        help=f"HTML report (default {DEFAULT_HOME}/reports/scan_<date>.html)",
    )
    a = ap.parse_args(argv)
    report = a.report or str(
        DEFAULT_HOME / "reports" / f"scan_{datetime.now():%Y-%m-%d}.html"
    )
    t0 = datetime.now()
    try:
        if _inside(Path(report).expanduser(), Path(a.root).expanduser()):
            raise FatalError(
                f"refusing to write the report inside the scanned folder: {report}"
            )
        print(f"{now()}  scanning {a.root}")
        counts = scan(a.root, a.db, a.mirror)
        out = write_report(a.db, a.root, report)
    except FatalError as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 2
    secs = (datetime.now() - t0).total_seconds()
    summary = {
        k: counts[k]
        for k in (
            "seen",
            "parsed",
            "cached",
            "mirrored_new",
            "read_errors",
            "newly_missing",
            "complete",
        )
    }
    print(f"{now()}  done in {secs:.1f}s: {json.dumps(summary)}")
    print(f"skipped: {json.dumps(counts['skipped'])}")
    print(f"report: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
