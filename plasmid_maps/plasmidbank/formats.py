"""Read one plasmid file into a MapRecord. Format is sniffed from content, never the extension.

Real-world quirks this handles (all seen in a real lab map bank):
- non-UTF-8 files (cp1252 / latin-1 fallback)
- `.seq` files that are really GenBank (Vector NTI exports)
- a bare-bases file with no header at all
- origin-wrapping feature locations (Biopython repairs them and warns)
- CLC `.clc` binaries (indexed, never parsed)
"""

from __future__ import annotations

import io
import warnings
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

from .ids import identify
from .sequtil import IUPAC_CHARS, canonical_hash

CLC_MAGIC = b"CLC\x00"


@dataclass
class Feature:
    type: str
    label: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive; end < start means the feature wraps the origin
    strand: int  # 1, -1 or 0

    @property
    def wraps(self) -> bool:
        return self.end < self.start


@dataclass
class MapRecord:
    path: str
    plasmid_id: str | None
    kind: str
    id_source: str
    xrefs: list[str]
    format: str  # genbank | fasta | raw | clc | unknown
    status: str  # ok | warn | error | clc_binary
    messages: list[str] = field(default_factory=list)
    name: str = ""
    description: str = ""
    seq: str = ""
    topology: str | None = None  # circular | linear | None (unknown)
    features: list[Feature] = field(default_factory=list)
    encoding: str = "utf-8"
    content_hash: str | None = None
    display_path: str = (
        ""  # original location when `path` is a mirror copy (index mode)
    )

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def circular(self) -> bool:
        return self.topology != "linear"


def decode(data: bytes) -> tuple[str, str]:
    for enc in ("utf-8-sig", "cp1252"):
        try:
            return data.decode(enc), enc
        except UnicodeDecodeError:
            continue
    return data.decode("latin-1"), "latin-1"


def sniff(data: bytes, text: str | None = None) -> str:
    if data.startswith(CLC_MAGIC):
        return "clc"
    t = (text if text is not None else decode(data)[0]).lstrip()
    if t.startswith("LOCUS"):
        return "genbank"
    if t.startswith(">"):
        return "fasta"
    body = "".join(t.split())
    if body and set(body.upper()) <= IUPAC_CHARS:
        return "raw"
    return "unknown"


def _label(feat) -> str:
    q = feat.qualifiers
    for key in ("label", "gene", "product", "note", "standard_name"):
        if key in q and q[key]:
            return str(q[key][0])
    return feat.type


def _span(parts, n: int, linear: bool = False) -> tuple[int, int]:
    """Outer span of a (possibly joined) location on a circle.

    The feature is taken to start after the largest uncovered gap; if that gap is
    the one spanning the origin, the feature does not wrap. On a linear record nothing wraps
    (2026-10-27: a two-exon gene on a 54 bp sequence read as wrapping, its intron being the widest gap).
    """
    segs = sorted((int(p.start), int(p.end)) for p in parts)
    if len(segs) == 1 or n == 0 or linear:
        return segs[0][0], segs[-1][1]
    gaps = [(segs[i + 1][0] - segs[i][1], i + 1) for i in range(len(segs) - 1)]
    origin_gap = segs[0][0] + n - segs[-1][1]
    widest, idx = max(gaps)
    if origin_gap >= widest:
        return segs[0][0], segs[-1][1]
    return segs[idx][0], segs[idx - 1][1]


def _features(rec) -> list[Feature]:
    out = []
    n = len(rec.seq)
    linear = rec.annotations.get("topology") == "linear"
    for f in rec.features:
        if f.type == "source" or f.location is None:
            continue
        start, end = _span(f.location.parts, n, linear)
        out.append(Feature(f.type, _label(f), start, end, f.location.strand or 0))
    return out


def read_map(path: str | Path) -> MapRecord:
    path = Path(path)
    data = path.read_bytes()
    text, enc = decode(data)
    fmt = sniff(data, text)
    first_line = text.lstrip().splitlines()[0] if text.strip() else ""
    locus = (
        first_line.split()[1]
        if fmt == "genbank" and len(first_line.split()) > 1
        else None
    )
    ids = identify(path.name, locus)
    rec = MapRecord(
        str(path),
        ids.plasmid_id,
        ids.kind,
        ids.source,
        ids.xrefs,
        fmt,
        "ok",
        encoding=enc,
    )
    if enc != "utf-8-sig":
        rec.messages.append(f"decoded as {enc}")
    if ids.source == "locus":
        rec.messages.append("ID taken from LOCUS line, not filename")

    if fmt == "clc":
        rec.status = "clc_binary"
        rec.messages.append("CLC binary: needs export to GenBank")
        return rec
    if fmt == "unknown":
        rec.status = "error"
        rec.messages.append("unrecognised format")
        return rec
    if fmt == "raw":
        rec.seq = "".join(text.split()).upper()
        rec.name = path.stem
        rec.messages.append("bare sequence, no header or features")
        rec.status = "warn"
    else:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                recs = list(SeqIO.parse(io.StringIO(text), fmt))
            except Exception as e:  # malformed files must still be indexed
                rec.status = "error"
                rec.messages.append(f"parse error: {type(e).__name__}: {e}")
                return rec
        for w in caught:
            rec.messages.append(f"parser: {w.message}")
        if not recs or not len(recs[0].seq):
            rec.status = "error"
            rec.messages.append("no sequence")
            return rec
        if len(recs) > 1:
            rec.messages.append(f"{len(recs)} records in file; using the first")
        r = recs[0]
        rec.seq = str(r.seq).upper()
        rec.name = r.name if fmt == "genbank" else r.id
        rec.description = r.description
        rec.topology = r.annotations.get("topology") if fmt == "genbank" else None
        rec.features = _features(r) if fmt == "genbank" else []
        if caught or len(recs) > 1:
            rec.status = "warn"
    bad = set(rec.seq) - IUPAC_CHARS
    if bad:
        rec.status = "warn"
        rec.messages.append(f"non-IUPAC characters in sequence: {''.join(sorted(bad))}")
    rec.content_hash = canonical_hash(rec.seq, rec.circular)
    return rec
