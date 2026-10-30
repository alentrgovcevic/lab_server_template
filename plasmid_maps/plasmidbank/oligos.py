"""Parse the lab oligo export (tab-separated; see `examples/oligos.tab`).

File quirks (profiled 2026-10-02, 20,199 rows):
- tab-separated, CRLF, header has 7 columns but rows carry an 8th (order date)
- sequences in mixed case; case carries no meaning we can rely on
- `[...]` blocks are modifications (e.g. phosphorylation, fluorophores) -> stripped, flagged
- degenerate IUPAC bases (N, R, Y ...) and non-nucleotide characters (X, I, O, P, mojibake)
- empty sequences; the same sequence ordered under several PM- ids
"""

from __future__ import annotations

import csv
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

from .sequtil import IUPAC_CHARS

_MOD_RE = re.compile(r"\[[^\]]*\]|/[^/]*/|\([^)]*\)")


@dataclass
class Oligo:
    oligo_id: str
    name: str
    reference: str
    oligo_type: str
    target: str
    raw_seq: str
    seq: str  # uppercase, modifications removed, non-IUPAC replaced by N
    comment: str
    date: str
    flags: set[str] = field(default_factory=set)


def normalise(raw: str) -> tuple[str, set[str]]:
    flags: set[str] = set()
    s = "".join(raw.split())
    if "[I]" in s.upper():  # inosine pairs with anything: keep the position as N
        flags.add("inosine")
        s = re.sub(r"\[I\]", "N", s, flags=re.IGNORECASE)
    if _MOD_RE.search(s):
        flags.add("modified")
        s = _MOD_RE.sub("", s)
    s = s.upper().replace("U", "T")
    out = []
    for c in s:
        if c in IUPAC_CHARS:
            out.append(c)
            if c not in "ACGT":
                flags.add("degenerate")
        elif c in "-*.":
            continue
        else:
            out.append("N")
            flags.add("nonstandard")
    seq = "".join(out)
    # free text in the sequence column ("use 69F", "sequencing", "xx")
    if (
        seq
        and (len(seq) < 8 or seq.count("N") > 0.3 * len(seq))
        and "nonstandard" in flags
    ):
        flags.add("not_a_sequence")
        flags.discard("degenerate")
        seq = ""
    if not seq:
        flags.add("empty")
    return seq, flags


def load_oligos(path: str | Path) -> list[Oligo]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    out = []
    for r in rows[1:]:
        if not any(c.strip() for c in r):
            continue
        r = r + [""] * (8 - len(r))
        seq, flags = normalise(r[5])
        out.append(
            Oligo(
                r[2].strip(),
                r[3].strip(),
                r[0].strip(),
                r[1].strip(),
                r[4].strip(),
                r[5],
                seq,
                r[6].strip(),
                r[7].strip(),
                flags,
            )
        )
    return out


def label_codes(raw: str) -> tuple[str, ...]:
    """The [..] modification codes of a sequence, in order and lowercased: '[6FAM]ACGT' -> ('[6fam]',)."""
    return tuple(c.lower() for c in re.findall(r"\[[^\]]*\]", "".join(raw.split())) if c.upper() != "[I]")


EXPORT_HEADER = ["Reference", "Oligo_type", "Oligo id ", "Oligo name", "Oligo target", "Oligo sequence", "Comment"]


def _rows(text: str) -> list[list[str]]:
    return [r for r in csv.reader(io.StringIO(text, newline=""), delimiter="\t") if any(c.strip() for c in r)]


def merge_export(current: str, upload: str) -> tuple[str, dict]:
    """Merge an uploaded oligo export (the full one or just the newest rows) into the current one.

    Rows are keyed by Oligo id. New ids are appended, an id whose sequence changed is replaced,
    and nothing is ever removed (a partial upload must not delete the rest). Returns the merged
    file text (tab-separated, CRLF, the current header) and the diff.
    """
    up = _rows(upload)
    if not up or "oligo id" not in "|".join(c.strip().lower() for c in up[0]):
        raise ValueError("not an oligo export: the first row must be the header (Reference, Oligo_type, Oligo id, ...)")
    cur = _rows(current) if current.strip() else [EXPORT_HEADER]
    out = [r + [""] * (8 - len(r)) for r in cur[1:]]
    at = {r[2].strip(): i for i, r in enumerate(out) if r[2].strip()}
    new, changed, same, no_id = [], [], 0, 0
    for r in up[1:]:
        r = r + [""] * (8 - len(r))
        oid = r[2].strip()
        if not oid:
            no_id += 1
            continue
        if oid not in at:
            at[oid] = len(out)
            out.append(r)
            new.append({"id": oid, "name": r[3].strip(), "sequence": r[5].strip()})
        elif normalise(out[at[oid]][5])[0] != normalise(r[5])[0]:
            changed.append({"id": oid, "name": r[3].strip(), "old": out[at[oid]][5].strip(), "new": r[5].strip()})
            out[at[oid]] = r
        else:
            same += 1
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\r\n")
    w.writerow(cur[0])
    w.writerows([r if r[7] else r[:7] for r in out])
    return buf.getvalue(), {"new": new, "changed": changed, "unchanged": same, "without_id": no_id,
                            "in_upload": len(up) - 1, "after": len(out)}
