"""Build a genome for the browser from a FASTA plus annotation tracks. Stdlib only.

One genome = one folder `<GENOMES_DIR>/<slug>/`:
  genome.fa       the FASTA, normalised on import: plain text, 60 bases per line, names cut at
                  the first whitespace (so any input, gzipped or ragged, can be read by seeking)
  genome.sqlite   index + features (schema below)
  meta.json       name, description, who added it and when

From the command line (an admin with files on the server):
    python -m genomebrowser.build --genomes-dir DIR --name "Barley Morex V3" --fasta genome.fa.gz \
        --track genes=genes.gff3.gz --track repeats=repeats.gff3 --track peaks=peaks.bed \
        [--deflines deflines.tsv]
The web app runs the same code as a subprocess:  python -m genomebrowser.build --job JOB.json

Schema (1-based, closed coordinates, like GFF):
  meta(key, value)
  seq(num, name, length, offset, linebases, linewidth)   a samtools-style .fai
  tracks(name, kind, format, source, n, added_by, added) kind: genes (gene models) | features
  feat(id, track, seqid, start, end, strand, type, fid, parent, name, attrs)
  feat_rt         rtree (32-bit float) on (seq number, start, end). Floats so chromosomes longer
                  than 2^31 bp fit; SQLite rounds stored boxes outward, so a query never misses a
                  feature, and every caller filters exactly on feat.start/end afterwards.
  names(name NOCASE, feat_id)    search
  gene_info(gene, defline)       descriptions: GFF description/Note/product, or a deflines TSV
  density(track, seqid, bin, n)  feature starts per BIN bp
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
import shutil
import sqlite3
import sys
import time
from datetime import datetime
from pathlib import Path
from urllib.parse import unquote

BIN = 100_000
LINE = 60
SLUG_RE = re.compile(r"^[a-z0-9][a-z0-9_-]{1,40}$")
TRACK_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 _.()+-]{0,39}$")
GENE_TYPES = ("gene", "mRNA", "transcript", "pseudogene", "ncRNA_gene")
DESC_KEYS = ("description", "Note", "product", "function")
log = logging.getLogger("genomebrowser.build")


class BuildError(Exception):
    """The input cannot be built; the message is shown to the person who uploaded it."""


def slugify(name: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")[:40].strip("-")
    if len(s) < 2:
        raise BuildError(f"name {name!r} is too short")
    return s


def _open_text(path: Path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    if gz:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, encoding="utf-8", errors="replace", newline="")


# ── FASTA ─────────────────────────────────────────────────────────────────────


def normalise_fasta(src: Path, dst: Path) -> list[tuple[str, int, int, int, int]]:
    """Copy any FASTA (gzip or plain, any line length) to `dst` with 60-column lines, and return
    its .fai rows (name, length, offset, linebases, linewidth)."""
    rows: list[tuple[str, int, int, int, int]] = []
    names: set[str] = set()
    cur: dict = {}
    buf = ""
    pos = 0
    keep = re.compile(r"[^A-Za-z*-]")
    with _open_text(src) as fin, open(dst, "wb") as out:

        def write(text: str) -> None:
            nonlocal pos
            b = text.encode("ascii")
            out.write(b)
            pos += len(b)

        def flush(final: bool) -> None:
            nonlocal buf
            n = len(buf) if final else len(buf) // LINE * LINE
            if n:
                write("".join(buf[i : i + LINE] + "\n" for i in range(0, n, LINE)))
            buf = buf[n:]

        def close_record() -> None:
            if cur:
                flush(True)
                rows.append((cur["name"], cur["length"], cur["offset"], LINE, LINE + 1))

        for lineno, line in enumerate(fin, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                close_record()
                parts = line[1:].split()
                if not parts:
                    raise BuildError(f"FASTA line {lineno}: a '>' header without a name")
                name = parts[0]
                if name in names:
                    raise BuildError(f"FASTA: the sequence name {name!r} appears twice")
                names.add(name)
                write(f">{name}\n")
                cur = {"name": name, "length": 0, "offset": pos}
            else:
                if not cur:
                    raise BuildError("not a FASTA file: it does not start with a '>' header")
                if not line.isascii():
                    raise BuildError(f"FASTA line {lineno}: not sequence text")
                seq = keep.sub("", line)
                buf += seq
                cur["length"] += len(seq)
                if len(buf) >= 1024 * LINE:
                    flush(False)
        close_record()
    if not rows:
        raise BuildError("the FASTA holds no sequences")
    return rows


# ── annotation files ──────────────────────────────────────────────────────────


def parse_attrs(col: str) -> dict[str, str]:
    d = {}
    for part in col.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k.strip()] = unquote(v.strip())
    return d


def sniff_format(path: Path) -> str:
    """'gff3' or 'bed', from the first data line."""
    with _open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c = line.rstrip("\r\n").split("\t")
            if len(c) >= 9 and c[3].isdigit() and c[4].isdigit():
                return "gff3"
            if len(c) >= 3 and c[1].isdigit() and c[2].isdigit():
                return "bed"
            break
    raise BuildError(
        f"{path.name}: not GFF3 or BED (tab-separated) as far as the first line shows"
    )


def feature_rows(path: Path, track: str, fmt: str):
    """Yield (track, seqid, start, end, strand, type, fid, parent, name, attrs), 1-based closed."""
    with _open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c = line.rstrip("\r\n").split("\t")
            if fmt == "gff3":
                if len(c) < 9:
                    continue
                a = parse_attrs(c[8])
                fid = a.get("ID")
                parent = (a.get("Parent") or "").split(",")[0] or None
                try:
                    start, end = int(c[3]), int(c[4])
                except ValueError:
                    continue
                yield (
                    track,
                    c[0],
                    start,
                    end,
                    c[6],
                    c[2],
                    fid,
                    parent,
                    a.get("Name") or fid,
                    c[8].strip(),
                )
            else:  # BED: 0-based half-open → 1-based closed
                if len(c) < 3:
                    continue
                try:
                    start, end = int(c[1]) + 1, int(c[2])
                except ValueError:
                    continue
                name = c[3] if len(c) > 3 and c[3] not in ("", ".") else None
                strand = c[5] if len(c) > 5 and c[5] in "+-" else "."
                yield (track, c[0], start, end, strand, track, name, None, name, "")


# ── database ──────────────────────────────────────────────────────────────────

SCHEMA = """
    CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);
    CREATE TABLE seq(num INTEGER PRIMARY KEY, name TEXT UNIQUE, length INT, offset INT, linebases INT, linewidth INT);
    CREATE TABLE tracks(name TEXT PRIMARY KEY, kind TEXT, format TEXT, source TEXT, n INT, added_by TEXT, added TEXT);
    CREATE TABLE feat(id INTEGER PRIMARY KEY, track TEXT, seqid TEXT, start INT, "end" INT, strand TEXT,
                      type TEXT, fid TEXT, parent TEXT, name TEXT, attrs TEXT);
    CREATE INDEX feat_parent ON feat(parent);
    CREATE INDEX feat_fid ON feat(fid);
    CREATE INDEX feat_track ON feat(track);
    CREATE VIRTUAL TABLE feat_rt USING rtree(id, c0, c1, s, e);
    CREATE TABLE names(name TEXT COLLATE NOCASE, feat_id INT);
    CREATE INDEX names_i ON names(name);
    CREATE TABLE gene_info(gene TEXT PRIMARY KEY, defline TEXT);
    CREATE TABLE density(track TEXT, seqid TEXT, bin INT, n INT);
    CREATE INDEX density_i ON density(track, seqid);
"""


def add_track(
    db: sqlite3.Connection,
    name: str,
    path: Path,
    kind: str,
    user: str = "",
    fmt: str | None = None,
) -> int:
    if not TRACK_RE.match(name):
        raise BuildError(
            f"track name {name!r}: use letters, digits, spaces, _ . ( ) + - (max 40)"
        )
    if db.execute("SELECT 1 FROM tracks WHERE name = ?", (name,)).fetchone():
        raise BuildError(f"this genome already has a track called {name!r}")
    if kind not in ("genes", "features"):
        raise BuildError(f"track kind must be genes or features, not {kind!r}")
    fmt = fmt or sniff_format(path)
    num = dict(db.execute("SELECT name, num FROM seq"))
    first = db.execute("SELECT coalesce(max(id), 0) FROM feat").fetchone()[0] + 1
    n, unknown, batch, t0 = 0, {}, [], time.time()

    def flush():
        nonlocal n, batch
        rows = [(first + n + i, *r) for i, r in enumerate(batch)]
        db.executemany("INSERT INTO feat VALUES (?,?,?,?,?,?,?,?,?,?,?)", rows)
        db.executemany(
            "INSERT INTO feat_rt VALUES (?,?,?,?,?)",
            [(r[0], num[r[2]], num[r[2]], r[3], r[4]) for r in rows],
        )
        n += len(rows)
        batch = []

    for row in feature_rows(path, name, fmt):
        if row[1] not in num:
            unknown[row[1]] = unknown.get(row[1], 0) + 1
            continue
        batch.append(row)
        if len(batch) >= 50_000:
            flush()
    if batch:
        flush()
    if not n:
        example = next(iter(unknown), None)
        hint = (
            f" (the file uses e.g. {example!r}; the FASTA has e.g. {next(iter(num))!r})"
            if example
            else ""
        )
        raise BuildError(f"{path.name}: no features on sequences of this genome{hint}")
    if unknown:
        log.warning(
            "%s: skipped features on sequences not in the genome: %s",
            name,
            dict(list(unknown.items())[:10]),
        )
    lo = first
    # exact duplicate lines (some published tracks repeat a large share of theirs) are dropped
    db.execute("""CREATE TEMP TABLE dup AS SELECT id FROM feat WHERE id >= ? AND id NOT IN (
                    SELECT min(id) FROM feat WHERE id >= ?
                    GROUP BY seqid, start, "end", strand, type, fid, parent, attrs)""", (lo, lo))
    dups = db.execute("SELECT count(*) FROM dup").fetchone()[0]
    if dups:
        db.execute("DELETE FROM feat_rt WHERE id IN (SELECT id FROM dup)")
        db.execute("DELETE FROM feat WHERE id IN (SELECT id FROM dup)")
        n -= dups
        log.info("%s: %s exact duplicate lines dropped", name, f"{dups:,}")
    db.execute("DROP TABLE dup")
    if kind == "genes":
        db.execute(
            """INSERT INTO names SELECT fid, id FROM feat WHERE id >= ? AND fid IS NOT NULL
                      AND type IN ('gene','mRNA','transcript','pseudogene','ncRNA_gene')""",
            (lo,),
        )
        db.execute(
            """INSERT INTO names SELECT name, id FROM feat WHERE id >= ? AND type = 'gene'
                      AND name IS NOT NULL AND name != fid""",
            (lo,),
        )
        _descriptions_from_attrs(db, lo)
        dens = "type = 'gene'"
    else:
        db.execute(
            "INSERT INTO names SELECT fid, id FROM feat WHERE id >= ? AND fid IS NOT NULL",
            (lo,),
        )
        dens = "1"
    db.execute(
        f"""INSERT INTO density SELECT track, seqid, start / ?, count(*) FROM feat
                   WHERE id >= ? AND {dens} GROUP BY 1, 2, 3""",
        (BIN, lo),
    )
    db.execute(
        "INSERT INTO tracks VALUES (?,?,?,?,?,?,?)",
        (name, kind, fmt, path.name, n, user, _now()),
    )
    log.info(
        "%s: %s features (%s, %s)  %.0fs", name, f"{n:,}", kind, fmt, time.time() - t0
    )
    return n


def _descriptions_from_attrs(db: sqlite3.Connection, lo: int) -> None:
    """gene_info from description/Note/product on the gene, or else on its first transcript/CDS."""
    best: dict[str, str] = {}
    parent_gene = {}
    for fid, typ, parent, attrs in db.execute(
        "SELECT fid, type, parent, attrs FROM feat WHERE id >= ? ORDER BY id", (lo,)
    ):
        a = parse_attrs(attrs) if attrs else {}
        gene = (
            fid
            if typ in ("gene", "pseudogene", "ncRNA_gene")
            else parent_gene.get(parent, parent)
        )
        if typ in ("mRNA", "transcript") and fid:
            parent_gene[fid] = parent
        text = next((a[k] for k in DESC_KEYS if a.get(k)), None)
        if gene and text and gene not in best:
            best[gene] = text[:500]
    db.executemany("INSERT OR IGNORE INTO gene_info VALUES (?,?)", best.items())


def add_deflines(db: sqlite3.Connection, path: Path) -> int:
    """A TSV of descriptions: `gene_or_transcript <tab> [kind <tab>] description` (Phytozome style)."""
    best: dict[str, str] = {}
    with _open_text(path) as fh:
        for line in fh:
            c = line.rstrip("\r\n").split("\t")
            if len(c) >= 2 and c[0] and c[-1]:
                best.setdefault(c[0], c[-1])
                best.setdefault(c[0].rsplit(".", 1)[0], c[-1])  # transcript .1 → gene
    db.executemany("INSERT OR REPLACE INTO gene_info VALUES (?,?)", best.items())
    return len(best)


def _require(*paths) -> None:
    for x in paths:
        if x is not None and not Path(x).is_file():
            raise BuildError(f"input file not found: {x}")


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _check_space(target: Path, need: int) -> None:
    target.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(target).free
    if free < need:
        raise BuildError(
            f"not enough disk space: {free / 1e9:.1f} GB free, about {need / 1e9:.1f} GB needed"
        )


def new_genome(
    genomes_dir: Path,
    name: str,
    fasta: Path,
    tracks: list[dict],
    user: str = "",
    description: str = "",
    deflines: Path | None = None,
) -> str:
    """Build a new genome folder next to the others and move it into place when complete."""
    _require(fasta, deflines, *(t["path"] for t in tracks))
    slug = slugify(name)
    final = genomes_dir / slug
    if final.exists():
        raise BuildError(f"a genome called {slug!r} exists already")
    work = genomes_dir / ".building" / slug
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    insize = fasta.stat().st_size
    with open(fasta, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    _check_space(
        work,
        insize * (5 if gz else 2)
        + sum(Path(t["path"]).stat().st_size for t in tracks) * 3,
    )
    try:
        t0 = time.time()
        rows = normalise_fasta(fasta, work / "genome.fa")
        log.info(
            "genome: %d sequences, %s bp  %.0fs",
            len(rows),
            f"{sum(r[1] for r in rows):,}",
            time.time() - t0,
        )
        db = sqlite3.connect(work / "genome.sqlite")
        db.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;" + SCHEMA)
        db.executemany(
            "INSERT INTO seq(name, length, offset, linebases, linewidth) VALUES (?,?,?,?,?)",
            rows,
        )
        for t in tracks:
            add_track(
                db,
                t["name"],
                Path(t["path"]),
                t.get("kind", "features"),
                user,
                t.get("format"),
            )
        if deflines:
            log.info("descriptions: %s", f"{add_deflines(db, deflines):,}")
        db.executemany(
            "INSERT INTO meta VALUES (?,?)",
            {"name": name, "built": _now(), "bin": str(BIN)}.items(),
        )
        db.commit()
        db.close()
        (work / "meta.json").write_text(
            json.dumps(
                {
                    "slug": slug,
                    "name": name,
                    "description": description,
                    "added_by": user,
                    "added": _now(),
                    "source": fasta.name,
                },
                indent=1,
            )
        )
        if final.exists():
            raise BuildError(f"a genome called {slug!r} appeared while building")
        work.rename(final)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    log.info("done: %s  %.0fs", final, time.time() - t0)
    return slug


def add_tracks(
    genome_dir: Path, tracks: list[dict], user: str = "", deflines: Path | None = None
) -> None:
    """Add tracks to an existing genome: copy the database, add, then swap it in atomically, so
    readers never see a half-written file."""
    _require(deflines, *(t["path"] for t in tracks))
    db_path = genome_dir / "genome.sqlite"
    tmp = genome_dir / "genome.sqlite.building"
    _check_space(
        genome_dir,
        db_path.stat().st_size
        + sum(Path(t["path"]).stat().st_size for t in tracks) * 3,
    )
    shutil.copyfile(db_path, tmp)
    try:
        db = sqlite3.connect(tmp)
        db.execute("PRAGMA journal_mode=OFF")
        for t in tracks:
            add_track(
                db,
                t["name"],
                Path(t["path"]),
                t.get("kind", "features"),
                user,
                t.get("format"),
            )
        if deflines:
            add_deflines(db, deflines)
        db.commit()
        db.close()
        tmp.replace(db_path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def drop_track(genome_dir: Path, track: str) -> None:
    db_path = genome_dir / "genome.sqlite"
    tmp = genome_dir / "genome.sqlite.building"
    shutil.copyfile(db_path, tmp)
    try:
        db = sqlite3.connect(tmp)
        if not db.execute("SELECT 1 FROM tracks WHERE name = ?", (track,)).fetchone():
            raise BuildError(f"no track {track!r}")
        db.execute(
            "DELETE FROM names WHERE feat_id IN (SELECT id FROM feat WHERE track = ?)",
            (track,),
        )
        db.execute(
            "DELETE FROM feat_rt WHERE id IN (SELECT id FROM feat WHERE track = ?)",
            (track,),
        )
        db.execute("DELETE FROM feat WHERE track = ?", (track,))
        db.execute("DELETE FROM density WHERE track = ?", (track,))
        db.execute("DELETE FROM tracks WHERE name = ?", (track,))
        db.commit()
        db.execute("VACUUM")
        db.close()
        tmp.replace(db_path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def run_job(job: dict, genomes_dir: Path) -> str:
    """One queued job from the web app (see imports.py). Returns a one-line result."""
    kind, user = job["kind"], job.get("user", "")
    defl = Path(job["deflines"]) if job.get("deflines") else None
    if kind == "new_genome":
        slug = new_genome(
            genomes_dir,
            job["name"],
            Path(job["fasta"]),
            job.get("tracks", []),
            user,
            job.get("description", ""),
            defl,
        )
        return f"genome {slug} is ready"
    if (
        not SLUG_RE.match(job.get("slug", ""))
        or not (genomes_dir / job["slug"]).is_dir()
    ):
        raise BuildError(f"no genome {job.get('slug')!r}")
    if kind == "add_tracks":
        add_tracks(genomes_dir / job["slug"], job.get("tracks", []), user, defl)
        return f"{len(job.get('tracks', []))} track(s) added"
    if kind == "drop_track":
        drop_track(genomes_dir / job["slug"], job["track"])
        return f"track {job['track']} removed"
    raise BuildError(f"unknown job kind {kind!r}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genomes-dir", default=None)
    ap.add_argument("--job", help="a job file written by the web app")
    ap.add_argument("--name")
    ap.add_argument("--fasta")
    ap.add_argument(
        "--genes", action="append", default=[], help="gene models GFF3 (track 'genes')"
    )
    ap.add_argument(
        "--track",
        action="append",
        default=[],
        metavar="NAME=FILE",
        help="an extra GFF3/BED track; repeatable",
    )
    ap.add_argument("--deflines")
    ap.add_argument("--description", default="")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stdout)
    from .store import genomes_dir as default_dir

    gdir = Path(a.genomes_dir).expanduser() if a.genomes_dir else default_dir()
    try:
        if a.job:
            print(run_job(json.loads(Path(a.job).read_text()), gdir))
            return 0
        if not (a.name and a.fasta):
            ap.error("--name and --fasta are required (or --job)")
        tracks = [
            {
                "name": "genes" if i == 0 else f"genes {i + 1}",
                "path": p,
                "kind": "genes",
            }
            for i, p in enumerate(a.genes)
        ]
        for t in a.track:
            n, _, p = t.partition("=")
            if not p:
                ap.error(f"--track {t!r}: use NAME=FILE")
            tracks.append({"name": n, "path": p, "kind": "features"})
        slug = new_genome(
            gdir,
            a.name,
            Path(a.fasta).expanduser(),
            tracks,
            "cli",
            a.description,
            Path(a.deflines).expanduser() if a.deflines else None,
        )
        print(f"genome {slug} is ready in {gdir / slug}")
        return 0
    except BuildError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
