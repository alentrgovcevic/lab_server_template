"""Import a third-party feature list (a nucleotide FASTA) into the common annotation list.

Written for pLannotate's SnapGene-derived features (github.com/mmcguffi/pLannotate, GPL-3.0): its
database bundle holds them as a BLAST database, so extract a FASTA first with BLAST+, and take the
labels and types from `snapgene.csv` (columns `sseqid,Feature,Type,Description`; v2 fixtures use
`sseqid,name,type,blurb`, both are read):

    blastdbcmd -db <bundle>/BLAST_dbs/snapgene -entry all > snapgene.fa
    python -m plasmidbank.annotation_import snapgene.fa --meta snapgene.csv --source "pLannotate (SnapGene features)"
    python -m plasmidbank.annotation_import snapgene.fa --meta snapgene.csv --source ... --write --by <name>

Dry run by default. Rules: a sequence already in the common list is skipped whatever its label (the
lab's own name wins); repeats inside the file are skipped; an unknown type becomes misc_feature with
the original type kept in the note. Matching stays exact (`annotations.annotate`), so a part that
differs by one base from the listed one is not found; pLannotate itself uses BLAST for that.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import Counter
from pathlib import Path

from .annotations import TYPES, Library, clean_seq

_TYPE_ALIASES = {
    "primer": "primer_bind",
    "origin": "rep_origin",
    "ori": "rep_origin",
    "polya": "polyA_signal",
}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out, name, seq = [], None, []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0] if line[1:].strip() else "", []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def read_meta(path: Path) -> dict[str, dict]:
    """sseqid -> {label, type, description}, from either known header."""
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh))
    meta = {}
    for r in rows:
        r = {k.strip().lower(): (v or "").strip() for k, v in r.items() if k}
        if r.get("sseqid"):
            meta[r["sseqid"]] = {
                "label": r.get("feature") or r.get("name") or "",
                "type": r.get("type", ""),
                "description": r.get("description") or r.get("blurb") or "",
            }
    return meta


def _type(raw: str) -> str:
    if raw in TYPES:
        return raw
    low = raw.lower().replace(" ", "_")
    for t in TYPES:
        if t.lower() == low:
            return t
    return _TYPE_ALIASES.get(low, "misc_feature")


def _label_from_id(sseqid: str) -> str:
    return re.sub(r"_\(\d+\)$", "", sseqid).replace(
        "_", " "
    )  # snapgene ids: CmR_(2) -> CmR


def plan(
    fasta: Path, meta: dict[str, dict], known_seqs: set[str], source: str
) -> tuple[list[dict], Counter]:
    """The entries to add, and why the others were skipped."""
    skipped, seen, add = Counter(), set(), []
    for sid, raw in read_fasta(fasta):
        try:
            seq = clean_seq(raw)
        except ValueError:
            skipped["not DNA, or outside the length limits"] += 1
            continue
        m = meta.get(sid, {})
        label = (m.get("label") or _label_from_id(sid))[:120]
        if seq in known_seqs:
            skipped["sequence already in the common list"] += 1
            continue
        if (seq, label) in seen:
            skipped["repeated in the file"] += 1
            continue
        seen.add((seq, label))
        t = _type(m.get("type", ""))
        note = "; ".join(
            x
            for x in (
                source,
                m.get("description", ""),
                f"type {m['type']}" if m.get("type") and t != m["type"] else "",
            )
            if x
        )
        add.append({"label": label, "type": t, "seq": seq, "note": note[:500]})
    return add, skipped


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("fasta", type=Path)
    ap.add_argument(
        "--meta",
        type=Path,
        help="CSV with sseqid + Feature/name + Type + Description/blurb",
    )
    ap.add_argument(
        "--source",
        default="pLannotate",
        help="written at the start of each entry's note",
    )
    ap.add_argument(
        "--state",
        type=Path,
        default=Path(os.environ.get("PLASMIDBANK_STATE", Path.home() / "PlasmidBank")),
    )
    ap.add_argument("--by", default="import", help="the name shown as added_by")
    ap.add_argument("--write", action="store_true", help="add them (default: dry run)")
    a = ap.parse_args(argv)
    lib = Library(a.state)
    known = {
        r["seq"] for r in lib.db.execute("SELECT seq FROM entries WHERE owner = ''")
    }
    add, skipped = plan(a.fasta, read_meta(a.meta) if a.meta else {}, known, a.source)
    print(f"{len(add)} to add to the common list in {a.state / 'annotations.sqlite'}")
    for t, n in Counter(e["type"] for e in add).most_common():
        print(f"  {n:6d}  {t}")
    for why, n in skipped.items():
        print(f"skipped {n}: {why}")
    for e in add[:10]:
        print(f"  e.g. {e['label']} ({e['type']}, {len(e['seq'])} bp) {e['note'][:80]}")
    if not a.write:
        print("dry run: nothing written (add --write)")
        return 0
    n = 0
    for e in add:
        try:
            lib.add(a.by, "common", e["label"], e["type"], e["seq"], e["note"])
            n += 1
        except ValueError:
            pass
    print(f"added {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
