"""Cleanup overview of the map bank on the share, from the index: what to delete, rename or move, and why.

    python -m plasmidbank.cleanup_report [--db PATH] [--out PATH]

Read-only: it reads the index the scan wrote (no share access) and writes an HTML page plus a TSV
with the same rows (default ~/PlasmidBank/reports/cleanup_<date>.html / .tsv). It changes nothing on
O; a person acts on the list. The scan report (`scan_<date>.html`) stays the full hygiene view; this
page groups only what someone has to act on, with a suggested action per file.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path, PurePosixPath

from .index import DEFAULT_DB, DEFAULT_HOME, Index
from .scan import CSS, FAVICON, _e, _rel, _table

READABLE = {"genbank", "fasta", "raw"}
STUB_BYTES = 64  # a CLC file this small is the 32-byte header with nothing in it
NOT_MAP_FOLDERS = {
    "qPCR amplicons": "qPCR amplicon, not a plasmid map"
}
COPY = re.compile(r"(?i)(\bcopy\b|\bkopi\b|\s\(\d+\)$|_2$)")

# (key, heading, what it means) in the order the page shows them; most urgent first
CATEGORIES = [
    ("empty", "Empty maps", "The file holds no sequence."),
    (
        "unreadable",
        "Unreadable files",
        "Not a map format, or a GenBank file the parser rejects.",
    ),
    (
        "no_id",
        "No PMC/EC number",
        "The filename carries no PMC or EC number, so the bank cannot place it.",
    ),
    (
        "copy",
        "Copies and clashes",
        "'copy', '(1)' or '_2' in the name: compared with the other files for the same ID.",
    ),
    (
        "odd_name",
        "Odd filenames",
        "Names that trip up scripts, Windows or exports. Low priority.",
    ),
    (
        "clc_only",
        "Needs a CLC export",
        "Every file for this ID is a CLC binary: nothing outside CLC can read it.",
    ),
    (
        "no_features",
        "No annotations",
        "Readable sequence with no features at all. Information only.",
    ),
]


def _stem(path: str) -> str:
    name = PurePosixPath(path).name
    while PurePosixPath(name).suffix.lower() in {
        ".clc",
        ".gb",
        ".gbk",
        ".fa",
        ".fasta",
        ".seq",
    }:
        name = name[: -len(PurePosixPath(name).suffix)]
    return name


def _no_id_reason(rel: str) -> tuple[str, str]:
    parts = PurePosixPath(rel).parts
    name = _stem(rel)
    if parts[0] in NOT_MAP_FOLDERS:
        return NOT_MAP_FOLDERS[
            parts[0]
        ], "move the folder out of the bank, or exclude it from the scan"
    if len(parts) == 1:
        return (
            "CLC data at the top of the bank (enzyme list, BLAST database…)",
            "move out of the bank",
        )
    if re.match(r"(?i)PMS[-_ ]?\d+", name):
        return (
            "PMS number: the index does not read PMS numbers yet",
            "none for the file (index change)",
        )
    if m := re.match(r"PM-0*(\d{3,5})(?!\d)", name):
        return (
            "'PM-' where 'PMC-' is expected (PM- is the oligo prefix)",
            f"check, then rename to PMC-{int(m.group(1)):05d}…",
        )
    if parts[0] == "ENSA" and (m := re.match(r"(\d{5})(?!\d)", name)):
        return "EC number without the 'EC' prefix", f"rename to EC{m.group(1)}…"
    return (
        "no PMC/EC number in the name",
        "rename with its number, or move out of the bank",
    )


def _odd_name(rel: str) -> list[str]:
    name = PurePosixPath(rel).name
    out = []
    if re.search(r"\.(clc|gbk?|dna)\.\w+$", name, re.I):
        out.append("double extension")
    if name != name.strip() or re.search(r"\s\.\w+$", name):
        out.append("space before the extension or at an end")
    if "�" in name:
        out.append("broken character (�): the name was mangled in a copy")
    elif re.search(r"[^\x00-\x7f]", name):
        out.append(
            "special character ("
            + "".join(sorted({c for c in name if ord(c) > 127}))
            + ")"
        )
    return out


def classify(
    rows, feature_counts: dict[int, int], clc_only: set[str], root: str
) -> list[dict]:
    """One dict per (file, problem): category, path, id, reason, action. Pure: rows are mappings."""
    by_id = defaultdict(list)
    for r in rows:
        if r["plasmid_id"]:
            by_id[r["plasmid_id"]].append(r)
    out = []

    def add(cat, r, reason, action):
        out.append(
            {
                "category": cat,
                "path": _rel(r["path"], root),
                "id": r["plasmid_id"] or "",
                "format": r["format"],
                "size": r["size"],
                "reason": reason,
                "action": action,
            }
        )

    for r in rows:
        rel = _rel(r["path"], root)
        pid = r["plasmid_id"]
        others = [o for o in by_id.get(pid, []) if o["path"] != r["path"]]
        holder = next(
            (o for o in others if o["format"] in READABLE and o["length"]), None
        ) or next(
            (o for o in others if o["format"] == "clc" and o["size"] > STUB_BYTES), None
        )
        if r["format"] == "clc" and (r["size"] or 0) <= STUB_BYTES:
            if holder:
                add(
                    "empty",
                    r,
                    "empty CLC file (header only)",
                    f"delete; {_rel(holder['path'], root)} holds the map",
                )
            else:
                add(
                    "empty",
                    r,
                    "empty CLC file (header only), and no other file for this ID",
                    "map lost on the share: look in Benchling or CLC, else delete",
                )
        elif r["format"] in READABLE and r["status"] != "error" and not r["length"]:
            add("empty", r, "no sequence", "re-export from the source, or delete")
        if r["status"] == "error":
            msgs = "; ".join(json.loads(r["messages"] or "[]"))
            if r["format"] == "genbank":
                add(
                    "unreadable",
                    r,
                    msgs.split("\n")[0],
                    "open in SnapGene/Benchling, fix and save again",
                )
            else:
                add(
                    "unreadable",
                    r,
                    f"not a map format ({PurePosixPath(rel).suffix or 'no extension'})",
                    "move out of the bank; export as GenBank first if the map is needed",
                )
        if not pid:
            add("no_id", r, *_no_id_reason(rel))
        elif COPY.search(_stem(rel)):
            same = [
                o for o in others if (o["format"] == "clc") == (r["format"] == "clc")
            ]
            key = "sha256" if r["format"] == "clc" else "content_hash"
            twin = next((o for o in same if r[key] and o[key] == r[key]), None)
            if twin:
                add(
                    "copy",
                    r,
                    f"identical to {_rel(twin['path'], root)}",
                    "delete the copy",
                )
            elif same:
                add(
                    "copy",
                    r,
                    "differs from " + ", ".join(_rel(o["path"], root) for o in same),
                    "someone who knows the plasmid picks one",
                )
            else:
                add(
                    "copy",
                    r,
                    "the only file of its kind for this ID",
                    "rename without the copy marker",
                )
        for why in _odd_name(rel):
            add("odd_name", r, why, "rename")
        if (
            r["format"] in READABLE
            and r["length"]
            and not feature_counts.get(r["id"], 0)
        ):
            add(
                "no_features",
                r,
                f"{r['length']:,} bp, no features",
                "annotate when next used",
            )
    for pid in sorted(clc_only):
        files = by_id.get(pid, [])
        if files and not all((f["size"] or 0) <= STUB_BYTES for f in files):
            for f in files:
                add(
                    "clc_only",
                    f,
                    "only CLC files for this ID",
                    "export to GenBank from CLC",
                )
    return out


def render(items: list[dict], scanned: str, db: str) -> str:
    by_cat = defaultdict(list)
    for it in items:
        by_cat[it["category"]].append(it)
    summary = [
        [_e(head), str(len(by_cat[k])), _e(what)] for k, head, what in CATEGORIES
    ]
    sections = [
        f"<h2>Summary</h2>{_table(['category', 'files', 'meaning'], summary, {1})}"
    ]
    for key, head, what in CATEGORIES:
        rows = by_cat[key]
        actions = defaultdict(int)
        for it in rows:
            actions[it["action"].split(";")[0].split(":")[0]] += 1
        tally = " · ".join(
            f"{_e(a)}: {n}" for a, n in sorted(actions.items(), key=lambda x: -x[1])
        )
        table = _table(
            ["file", "ID", "why", "suggested action"],
            [
                [
                    f'<span class="path">{_e(it["path"])}</span>',
                    _e(it["id"]),
                    _e(it["reason"]),
                    _e(it["action"]),
                ]
                for it in rows
            ],
        )
        sections.append(
            f"<h2 id='{key}'>{_e(head)} ({len(rows):,})</h2><p class='note'>{_e(what)}</p>"
            + (f"<p class='muted'>{tally}</p>" if rows else "")
            + table
        )
    title = "Map bank cleanup"
    return (
        '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<meta name="color-scheme" content="light dark">\n'
        f'<title>{title}</title>\n<link rel="icon" href="{FAVICON}">\n'
        f"<style>{CSS}</style>\n</head>\n<body>\n<main>\n"
        f"<h1>{title}</h1>\n<p class='muted'>From the index <code>{_e(db)}</code>, scan of {_e(scanned)} · "
        f"written {datetime.now():%Y-%m-%d %H:%M}. Read-only: nothing on the share was changed. "
        "A file can appear in more than one section.</p>\n"
        + "\n".join(f"<section>{s}</section>" for s in sections)
        + "\n</main>\n</body>\n</html>\n"
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m plasmidbank.cleanup_report", description=__doc__.split("\n")[0]
    )
    ap.add_argument(
        "--db", default=str(DEFAULT_DB), help=f"index (default {DEFAULT_DB})"
    )
    ap.add_argument(
        "--out",
        default=None,
        help=f"HTML (default {DEFAULT_HOME}/reports/cleanup_<date>.html)",
    )
    a = ap.parse_args(argv)
    idx = Index(a.db)
    run = idx.last_run()
    if not run:
        print("FATAL: the index has no scan run", file=sys.stderr)
        return 2
    root = run["root"]
    rows = idx.files(root)
    counts = dict(
        idx.db.execute(
            "SELECT file_id, COUNT(*) FROM feature GROUP BY file_id"
        ).fetchall()
    )
    items = classify(rows, counts, set(idx.clc_only(root)), root)
    idx.close()
    out = Path(
        a.out or DEFAULT_HOME / "reports" / f"cleanup_{datetime.now():%Y-%m-%d}.html"
    ).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        render(items, f"{root} ({run['started'][:16]})", a.db), encoding="utf-8"
    )
    with open(out.with_suffix(".tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            ["category", "path", "id", "format", "size", "reason", "action"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(items)
    per = defaultdict(int)
    for it in items:
        per[it["category"]] += 1
    print(", ".join(f"{k} {per[k]}" for k, _, _ in CATEGORIES))
    print(f"report: {out}\ntsv:    {out.with_suffix('.tsv')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
