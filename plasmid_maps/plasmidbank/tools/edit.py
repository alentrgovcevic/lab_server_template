"""Apply edits from the app to a GenBank file and return the new text. Never writes the original.

Edits (all optional):
  {"name": "...", "definition": "...", "topology": "circular"|"linear",
   "features": [{"index": i, "action": "delete"} |
                {"index": i, "action": "update", "type": ..., "label": ..., "start": 1-based, "end": 1-based, "strand": 1|-1|0}],
   "add": [{"type": ..., "label": ..., "start": 1-based, "end": 1-based inclusive, "strand": 1|-1|0}]}
`index` counts non-source features in file order, the same order the app shows (formats._features).

An "add" entry may carry "parts": [[start, end], …] instead of start/end: one joined feature (a CDS
over exons), parts listed 5'->3' on its strand.

Sequence edits ("edit basepairs"), applied after the feature edits, in the
coordinates of the sequence as it was opened:
   "seq": [{"start": 1-based first base replaced, "end": last base replaced (start - 1 = insert
            before start), "bases": "ACGT…" ("" = delete)}]
Features after an edit shift; a feature with the edit inside it grows or shrinks; one partly covered
is trimmed; one lying inside a longer replaced or deleted stretch is dropped (listed in the summary).
Replacing exactly a feature's bases keeps the feature over the new bases.

"reverse_complement": true flips the whole map, features included, after all other edits.
"""

from __future__ import annotations

import io
import re
import warnings
from datetime import date

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

from .clean import read_genbank


def _location(start: int, end: int, strand: int, L: int, circular: bool):
    if not (1 <= start <= L and 1 <= end <= L):
        raise ValueError(f"positions must be within 1..{L}")
    strand = strand or None
    if end >= start:
        return SimpleLocation(start - 1, end, strand)
    if not circular:
        raise ValueError("end before start on a linear sequence")
    parts = [SimpleLocation(start - 1, L, strand), SimpleLocation(0, end, strand)]
    return CompoundLocation(parts if strand != -1 else parts[::-1])


IUPAC = set("ACGTURYSWKMBDHVN")
FEATURE_TYPE = re.compile(r"^[A-Za-z][A-Za-z0-9_'\-]{0,19}$")


def _ftype(t: str) -> str:
    """A feature type: a standard INSDC key or the lab's own. Plain characters only, so
    the GenBank file stays readable by other programs (some show an unknown type as misc_feature)."""
    t = (t or "").strip()
    if not FEATURE_TYPE.match(t) or t == "source":
        raise ValueError(f"'{t}' is not a usable feature type: letters, digits, _ - ' (up to 20), starting with a letter")
    return t
MAX_SEQ_EDIT = 100_000


def _shift_part(ps: int, pe: int, a: int, b: int, n: int) -> tuple[int, int] | None:
    """One 0-based half-open part [ps, pe) after [a, b) became n bases; None when it is gone."""
    d = n - (b - a)
    if pe <= a:  # before the edit, or an insertion right after its end
        return ps, pe
    if ps >= b:
        return ps + d, pe + d
    if ps <= a and b <= pe:
        return (ps, pe + d) if pe + d > ps else None
    if a <= ps and pe <= b:
        return None
    if ps < a:
        return ps, a
    return a + n, pe + d


def _edit_sequence(rec, feats: list, edits: list, circular: bool, summary: list[str]) -> list:
    L = len(rec.seq)
    seq = str(rec.seq)
    spans = []
    for e in edits:
        a, b = int(e["start"]) - 1, int(e["end"])
        bases = str(e.get("bases") or "").upper().replace(" ", "").replace("\n", "")
        if not (0 <= a <= b <= L):
            raise ValueError(f"a sequence edit must lie within 1..{L} (start <= end + 1)")
        if set(bases) - IUPAC or len(bases) > MAX_SEQ_EDIT:
            raise ValueError("bases must be IUPAC letters (A C G T, N, …), at most 100,000")
        if a == b and not bases:
            continue
        spans.append((a, b, bases))
    spans.sort(reverse=True)
    for (a, b, _), (a2, b2, _) in zip(spans, spans[1:]):  # sorted high to low: the next one is earlier
        if b2 > a or (a2 == b2 == a == b):
            raise ValueError("sequence edits overlap")
    for a, b, bases in spans:
        seq = seq[:a] + bases + seq[b:]
        kept = []
        for f in feats:
            parts = []
            for loc in f.location.parts:
                r = _shift_part(int(loc.start), int(loc.end), a, b, len(bases))
                if r:
                    parts.append(SimpleLocation(r[0], r[1], loc.strand))
            if not parts:
                summary.append(f"dropped {f.type} '{(f.qualifiers.get('label') or [''])[0]}' (its bases were replaced)")
                continue
            f.location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
            kept.append(f)
        feats = kept
        what = f"deleted {b - a} bp at {a + 1}" if not bases else (
            f"inserted {len(bases)} bp before {a + 1}" if a == b else f"replaced {a + 1}-{b} ({b - a} bp) with {len(bases)} bp")
        summary.append(what)
    rec.letter_annotations = {}
    rec.seq = Seq(seq)
    for f in rec.features:
        if f.type == "source":
            f.location = SimpleLocation(0, len(seq), f.location.strand)
    return feats


def apply_edits(path: str, edits: dict) -> tuple[str, list[str]]:
    rec = read_genbank(path)
    seq_before = str(rec.seq)
    L = len(rec.seq)
    summary: list[str] = []
    if edits.get("topology") in ("circular", "linear") and edits[
        "topology"
    ] != rec.annotations.get("topology"):
        summary.append(
            f"topology {rec.annotations.get('topology')} -> {edits['topology']}"
        )
        rec.annotations["topology"] = edits["topology"]
    circular = rec.annotations.get("topology") == "circular"
    if edits.get("name"):
        summary.append(f"name {rec.name} -> {edits['name']}")
        rec.name = rec.id = str(edits["name"])[:60]
    if edits.get("definition") is not None and edits["definition"] != rec.description:
        rec.description = str(edits["definition"])
        summary.append("definition changed")
    feats = [f for f in rec.features if f.type != "source"]
    drop = set()
    for e in edits.get("features", []):
        i = int(e["index"])
        if not 0 <= i < len(feats):
            raise ValueError(f"no feature with index {i}")
        f = feats[i]
        if e.get("action") == "delete":
            drop.add(i)
            summary.append(
                f"deleted {f.type} '{(f.qualifiers.get('label') or [''])[0]}'"
            )
        elif e.get("action") == "update":
            if e.get("type"):
                f.type = _ftype(e["type"])
            if e.get("label") is not None:
                f.qualifiers["label"] = [e["label"]]
            if e.get("start") and e.get("end"):
                f.location = _location(
                    int(e["start"]),
                    int(e["end"]),
                    int(e.get("strand", f.location.strand or 0)),
                    L,
                    circular,
                )
            summary.append(
                f"updated feature {i}: {f.type} '{(f.qualifiers.get('label') or [''])[0]}' {f.location}"
            )
    for a in edits.get("add", []):
        strand = int(a.get("strand", 1))
        if a.get("parts"):  # a joined feature, e.g. exons as one CDS: parts in 5'->3' order
            locs = [_location(int(x), int(y), strand, L, circular) for x, y in a["parts"]]
            loc = locs[0] if len(locs) == 1 else CompoundLocation([pt for lc in locs for pt in lc.parts])
        else:
            loc = _location(int(a["start"]), int(a["end"]), strand, L, circular)
        f = SeqFeature(
            loc,
            type=_ftype(a.get("type") or "misc_feature"),
            qualifiers={"label": [a.get("label") or "new feature"]},
        )
        feats.append(f)
        summary.append(f"added {f.type} '{f.qualifiers['label'][0]}' {f.location}")
    feats = [f for i, f in enumerate(feats) if i not in drop]
    if edits.get("seq"):
        feats = _edit_sequence(rec, feats, edits["seq"], circular, summary)
        seq_before = str(rec.seq)
    source = [f for f in rec.features if f.type == "source"]
    rec.features = source + feats
    if edits.get("reverse_complement"):  # last, so every edit above used the coordinates as opened
        rc = rec.reverse_complement(id=rec.id, name=rec.name, description=rec.description, features=True,
                                    annotations=True, dbxrefs=True)
        rec = rc
        seq_before = str(rec.seq)
        summary.append("reverse complemented")
    rec.annotations["comment"] = (
        (
            rec.annotations.get("comment", "") + "\n"
            if rec.annotations.get("comment")
            else ""
        )
        + f"edited in plasmidbank {date.today().isoformat()}: "
        + "; ".join(summary)
    )
    if str(rec.seq) != seq_before:
        raise RuntimeError("sequence changed during editing - refusing to save")
    out = io.StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SeqIO.write(rec, out, "genbank")
    return out.getvalue(), summary
