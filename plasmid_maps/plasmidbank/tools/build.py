"""Build a new map from a vector and inserts: Gibson design, or documenting a finished cloning.

The vector (one map) is opened either by restriction enzymes or by PCR between two positions. The
inserts (up to MAX_INSERTS) are each a region of a map in the bank, on either strand, or a pasted
sequence (a synthetic fragment: it gets no primers). The product is the circular joining

    backbone · insert 1 · … · insert N

with the features of every piece carried over, so it can be downloaded, saved to a project, or
uploaded to the bank with its lab ID.

mode "document": no primers, just the map of a cloning that is already done. With enzymes, the
    sticky-end ligation keeps the top strand: backbone from the right site's top cut to the left
    site's top cut, and the insert is exactly the bases that end up between them.
mode "gibson": overlaps at every junction (OVERLAP nt by default). The exonuclease chews 5' ends,
    so an enzyme-opened backbone keeps the 3' ends: from the right site's bottom cut to the left
    site's top cut (5' overhangs are lost, 3' overhangs kept), which is what the overlaps match.
    At each junction the overlap goes on the PCR pieces' primer tails: split in two when both
    sides are PCR products, whole on one side otherwise; a synthetic fragment next to an
    enzyme-opened vector is extended instead ("order this sequence").

Positions are 1-based map coordinates; strands '+' / '-'.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from ..formats import MapRecord
from ..sequtil import revcomp
from . import thermo
from .cloning import Product, _carry_features, find_cuts, product_genbank
from .design import _primer, _scope
from .mutagenesis import _anneal

MAX_INSERTS = 5
OVERLAP = 25
OVERLAP_RANGE = (15, 40)
CROSS_K = 12  # two overlaps sharing a 12-mer can pair with the wrong partner


@dataclass
class Piece:
    name: str
    seq: str
    pcr: bool  # amplified, so its primers can carry tails
    features: list[dict] = field(default_factory=list)
    # where its primers bind, for the specificity check: (ref, fwd (start, end, strand, 3'), rev (...))
    ref: MapRecord | None = None
    fwd_site: tuple | None = None
    rev_site: tuple | None = None
    synthetic: bool = False


def _clean_seq(text: str) -> str:
    lines = [ln for ln in text.splitlines() if not ln.startswith(">")]
    s = re.sub(r"[\s\d/]", "", "".join(lines)).upper()
    bad = sorted(set(s) - set("ACGTRYSWKMBDHVN"))
    if bad:
        raise ValueError(
            f"the sequence has characters that are not DNA: {''.join(bad)}"
        )
    if not s:
        raise ValueError("empty sequence")
    return s


def _flip(feats: list[dict], length: int) -> list[dict]:
    return [
        {
            **f,
            "start": length - f["end"],
            "end": length - f["start"],
            "strand": -f["strand"] if f["strand"] else 0,
        }
        for f in feats
    ]


def _span(L: int, start: int, end: int) -> int:
    """Length of the 1-based closed stretch start..end on a circle (end < start wraps)."""
    return end - start + 1 if end >= start else end + L - start + 1


# ── the pieces ────────────────────────────────────────────────────────────────


def open_vector(
    ref: MapRecord,
    mode: str,
    enzymes: list[str] | None = None,
    start: int = 0,
    end: int = 0,
) -> tuple[Piece, dict]:
    """The backbone, read from the end of the removed stretch round to its start."""
    if not ref.seq or not ref.circular:
        raise ValueError("the vector must be a circular map with a sequence")
    L, s2 = ref.length, ref.seq.upper() * 2
    if enzymes:
        cuts = find_cuts(ref.seq.upper(), True, enzymes)
        per = {e: [c for c in cuts if c["enzyme"] == e] for e in enzymes}
        for e, cs in per.items():
            if len(cs) != 1:
                raise ValueError(
                    f"{e} cuts {ref.plasmid_id or 'the vector'} {len(cs)} times; pick enzymes that cut once"
                )
        sites = sorted((c for cs in per.values() for c in cs), key=lambda c: c["start"])

        def top0(c):
            return c["start"] + (c["len"] if c["kind"] == "3'" else 0)

        def bottom0(c):
            return c["start"] + (c["len"] if c["kind"] == "5'" else 0)

        if (
            len(sites) == 1
        ):  # linearised: backbone is the whole circle from the cut back to it
            a = b = sites[0]
            options = [(a, b)]
        else:  # two sites: keep the longer of the two stretches, drop the other
            a, b = sites
            options = [(a, b), (b, a)]

        def keep(left, right):  # backbone runs from `right` round to `left`
            lo = bottom0(right) if mode == "gibson" else top0(right)
            hi = top0(left)
            n = (hi - lo) % L or L
            return lo % L, n  # a cut near the origin can put lo past L

        left, right = max(options, key=lambda o: keep(*o)[1])
        lo, n = keep(left, right)
        seq = s2[lo : lo + n]
        info = {
            "opened_by": " + ".join(enzymes),
            "dropped_bp": L - n,
            "kept_from": lo % L + 1,
            "kept_to": (lo + n - 1) % L + 1,
        }
        return Piece(
            ref.plasmid_id or ref.name or "vector",
            seq,
            False,
            _carry_features(ref, lo, n),
        ), info
    # PCR-opened: remove start..end (end = start - 1: nothing removed, the cut is before start)
    if not (1 <= start <= L + 1 and start - 1 <= end <= L):
        raise ValueError(
            f"vector positions must be within 1..{L} (end = start - 1 opens it before start)"
        )
    lo, n = end % L, L - (end - start + 1)
    if n < 40:
        raise ValueError("almost nothing of the vector would be left")
    seq = s2[lo : lo + n]
    p = Piece(
        ref.plasmid_id or ref.name or "vector",
        seq,
        True,
        _carry_features(ref, lo, n),
        ref,
    )
    p.fwd_site = ("+", lo % L + 1)  # forward primer starts at the first kept base
    p.rev_site = (
        "-",
        (lo + n - 1) % L + 1,
    )  # reverse primer's 5' end at the last kept base
    return p, {
        "opened_by": "PCR",
        "dropped_bp": L - n,
        "kept_from": lo % L + 1,
        "kept_to": (lo + n - 1) % L + 1,
    }


def region_insert(
    ref: MapRecord, start: int, end: int, strand: str = "+", name: str = ""
) -> Piece:
    if not ref.seq:
        raise ValueError(f"{ref.plasmid_id or 'that map'} has no sequence")
    L = ref.length
    if not (1 <= start <= L and 1 <= end <= L) or (not ref.circular and end < start):
        raise ValueError(
            f"insert positions must be within 1..{L} of {ref.plasmid_id or 'its map'}"
        )
    n = _span(L, start, end)
    s2 = ref.seq.upper() * 2 if ref.circular else ref.seq.upper()
    seq = s2[start - 1 : start - 1 + n]
    feats = _carry_features(ref, start - 1, n)
    label = name or f"{ref.plasmid_id or 'insert'} {start}-{end}"
    p = Piece(label, seq, True, feats, ref)
    first, last = start, (start - 1 + n - 1) % L + 1
    if strand == "-":
        p.seq, p.features = revcomp(seq), _flip(feats, n)
        p.fwd_site, p.rev_site = ("-", last), ("+", first)
    else:
        p.fwd_site, p.rev_site = ("+", first), ("-", last)
    return p


def pasted_insert(text: str, name: str = "") -> Piece:
    return Piece(name or "synthetic fragment", _clean_seq(text), False, synthetic=True)


# ── primers and the product ───────────────────────────────────────────────────


def _bind(piece: Piece, end: str, anneal: str, tail: str, hits) -> dict:
    """A primer at the piece's start ('fwd') or end ('rev'), with its specificity check."""
    ref = piece.ref
    strand, pos = piece.fwd_site if end == "fwd" else piece.rev_site
    L, n = ref.length, len(anneal)
    if strand == "+":  # 5' end at pos, reading right
        s, e = pos, (pos - 1 + n - 1) % L + 1
        three = e
    else:  # 5' end at pos, reading left on the bottom strand
        e, s = pos, (pos - n) % L + 1
        three = s
    return _primer(anneal, strand, s, e, three, tail, ref, hits(ref))


def _cross_talk(overlaps: list[str]) -> list[str]:
    out = []
    kmers = [
        {o[i : i + CROSS_K] for i in range(len(o) - CROSS_K + 1)}
        | {revcomp(o[i : i + CROSS_K]) for i in range(len(o) - CROSS_K + 1)}
        for o in overlaps
    ]
    for i in range(len(overlaps)):
        for j in range(i + 1, len(overlaps)):
            if {
                overlaps[j][k : k + CROSS_K]
                for k in range(len(overlaps[j]) - CROSS_K + 1)
            } & kmers[i]:
                out.append(
                    f"overlaps {i + 1} and {j + 1} share a {CROSS_K}-nt stretch: the fragments may join in the wrong order"
                )
    return out


def build(
    vector: MapRecord,
    vec_open: dict,
    inserts: list[Piece],
    mode: str = "gibson",
    overlap: int = OVERLAP,
    opt_tm: float = 60.0,
    size: tuple[int, int] | None = None,
    hits=lambda ref: [],
    name: str = "",
) -> dict:
    if mode not in ("gibson", "document"):
        raise ValueError("mode must be gibson or document")
    if not 1 <= len(inserts) <= MAX_INSERTS:
        raise ValueError(f"give 1 to {MAX_INSERTS} inserts")
    if mode == "gibson" and not OVERLAP_RANGE[0] <= overlap <= OVERLAP_RANGE[1]:
        raise ValueError(f"overlap must be {OVERLAP_RANGE[0]}-{OVERLAP_RANGE[1]} nt")
    backbone, vinfo = open_vector(
        vector,
        mode,
        vec_open.get("enzymes"),
        vec_open.get("start", 0),
        vec_open.get("end", 0),
    )
    pieces = [backbone, *inserts]
    size = size or (18, 35)
    warnings: list[str] = []

    # the product: pieces joined in order, features shifted, one feature per insert
    seq, feats, parts, offsets = "", [], [], []
    for p in pieces:
        offsets.append(len(seq))
        feats += [
            {**f, "start": f["start"] + len(seq), "end": f["end"] + len(seq)}
            for f in p.features
        ]
        parts.append(
            {
                "source": p.name,
                "orientation": "+",
                "start": len(seq),
                "length": len(p.seq),
            }
        )
        seq += p.seq
    P = len(seq)
    for p, off in zip(inserts, offsets[1:]):
        feats.append(
            {
                "type": "misc_feature",
                "label": f"insert: {p.name}",
                "start": off,
                "end": off + len(p.seq),
                "strand": 1,
                "source": "assembly",
            }
        )

    junctions, primers, order = [], [], []
    tails = {i: {"fwd": "", "rev": ""} for i in range(len(pieces))}
    extend = {
        i: ["", ""] for i in range(len(pieces))
    }  # synthetic fragments: bases added left, right
    if mode == "gibson":
        for i in range(len(pieces)):
            lp, rp = pieces[i], pieces[(i + 1) % len(pieces)]
            li, ri = i, (i + 1) % len(pieces)
            if lp.pcr and rp.pcr:
                h = overlap // 2
                tails[li]["rev"] = revcomp(rp.seq[:h])
                tails[ri]["fwd"] = lp.seq[-(overlap - h) :]
                ov = lp.seq[-(overlap - h) :] + rp.seq[:h]
            elif rp.pcr:
                tails[ri]["fwd"] = lp.seq[-overlap:]
                ov = lp.seq[-overlap:]
            elif lp.pcr:
                tails[li]["rev"] = revcomp(rp.seq[:overlap])
                ov = rp.seq[:overlap]
            elif rp.synthetic:
                extend[ri][0] = lp.seq[-overlap:]
                ov = lp.seq[-overlap:]
            elif lp.synthetic:
                extend[li][1] = rp.seq[:overlap]
                ov = rp.seq[:overlap]
            else:  # cannot happen: the backbone is the only non-PCR, non-synthetic piece
                raise ValueError(
                    "two pieces that cannot carry an overlap are next to each other"
                )
            junctions.append(
                {
                    "overhang": ov,
                    "kind": "Gibson overlap",
                    "between": f"{lp.name} | {rp.name}",
                    "tm": thermo.tm(ov),
                    "gc": round(100 * sum(c in "GC" for c in ov) / len(ov), 1),
                }
            )
            if (thermo.tm(ov) or 0) < 48:
                warnings.append(
                    f"overlap {i + 1} ({lp.name} | {rp.name}) has Tm {thermo.tm(ov)} °C; lengthen the overlap"
                )
        warnings += _cross_talk([j["overhang"] for j in junctions])

        for i, p in enumerate(pieces):
            if not p.pcr:
                continue
            fa = _anneal(p.seq, opt_tm, size)
            ra = _anneal(revcomp(p.seq), opt_tm, size)
            f = _bind(p, "fwd", fa, tails[i]["fwd"], hits)
            r = _bind(p, "rev", ra, tails[i]["rev"], hits)
            for pr, d in ((f, "F"), (r, "R")):
                pr["piece"], pr["dir"] = p.name, d
                if len(pr["full"]) > 60:
                    warnings.append(
                        f"{p.name} {d} primer is {len(pr['full'])} nt: order HPLC-purified"
                    )
                if pr["other_sites"]:
                    warnings.append(
                        f"{p.name} {d} primer also binds elsewhere on {pr.get('ref', p.ref.plasmid_id)}"
                    )
            primers += [f, r]
            # where the whole primer sits on the product
            off = offsets[i]
            feats.append(
                {
                    "type": "primer_bind",
                    "label": f"{p.name} F",
                    "start": (off - len(f["tail"])) % P,
                    "end": (off + len(fa) - 1) % P + 1,
                    "strand": 1,
                    "source": "assembly",
                }
            )
            e0 = off + len(p.seq)
            feats.append(
                {
                    "type": "primer_bind",
                    "label": f"{p.name} R",
                    "start": (e0 - len(ra)) % P,
                    "end": (e0 + len(r["tail"]) - 1) % P + 1,
                    "strand": -1,
                    "source": "assembly",
                }
            )
        for i, p in enumerate(pieces):
            if p.synthetic:
                order.append(
                    {
                        "name": p.name,
                        "sequence": extend[i][0] + p.seq + extend[i][1],
                        "added_left": len(extend[i][0]),
                        "added_right": len(extend[i][1]),
                    }
                )
    else:
        for i in range(len(pieces)):
            lp, rp = pieces[i], pieces[(i + 1) % len(pieces)]
            junctions.append(
                {"overhang": "", "kind": "join", "between": f"{lp.name} | {rp.name}"}
            )

    product = Product(seq, True, parts, junctions, feats, warnings)
    label = (
        name
        or f"{vector.plasmid_id or 'vector'}_{'_'.join(p.name.split()[0] for p in inserts)}"
    )
    return {
        "mode": mode,
        "vector": vinfo,
        "length": P,
        "seq": seq,
        "pieces": [
            {
                "name": p.name,
                "length": len(p.seq),
                "pcr": p.pcr,
                "synthetic": p.synthetic,
            }
            for p in pieces
        ],
        "junctions": junctions,
        "primers": primers,
        "order_fragments": order,
        "warnings": warnings,
        "genbank": product_genbank(
            product, re.sub(r"[^\w.-]+", "_", label)[:16] or "assembly", label
        ),
        "name": label,
        "checked_against": _scope(vector),
        "tm_model": thermo.MODEL,
    }


# ── a map from a pasted sequence ──────────────────────────────────────────────


def map_from_text(text: str, name: str, circular: bool = True) -> dict:
    """GenBank text (features kept), FASTA or bare sequence -> a GenBank file."""
    from .cloning import Product as _P

    t = text.strip()
    if t.startswith("LOCUS"):
        import io

        from Bio import SeqIO

        rec = SeqIO.read(io.StringIO(t), "genbank")
        seq = _clean_seq(str(rec.seq))
        feats = []
        for f in rec.features:
            if f.type == "source":
                continue
            lab = (
                f.qualifiers.get("label")
                or f.qualifiers.get("gene")
                or f.qualifiers.get("product")
                or [f.type]
            )[0]
            feats.append(
                {
                    "type": f.type,
                    "label": lab,
                    "start": int(f.location.start),
                    "end": int(f.location.end),
                    "strand": f.location.strand or 0,
                    "source": "pasted GenBank",
                }
            )
        circular = (
            rec.annotations.get("topology", "circular" if circular else "linear")
            == "circular"
        )
        name = name or rec.name
    else:
        if t.startswith(">") and not name:
            name = (
                t.splitlines()[0][1:].strip().split()[0]
                if t.splitlines()[0][1:].strip()
                else ""
            )
        seq, feats = _clean_seq(t), []
    if len(seq) < 20:
        raise ValueError("that is too short for a map")
    name = name or "pasted"
    p = _P(
        seq,
        circular,
        [
            {
                "source": "pasted sequence",
                "orientation": "+",
                "start": 0,
                "length": len(seq),
            }
        ],
        [],
        feats,
    )
    gb = product_genbank(p, re.sub(r"[^\w.-]+", "_", name)[:16], name)
    return {
        "name": name,
        "length": len(seq),
        "circular": circular,
        "features": len(feats),
        "genbank": gb,
    }
