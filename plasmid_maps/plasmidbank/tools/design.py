"""Primer design on a map region (primer3), circular-aware, checked against the lab collection.

Tasks:
- cloning:    primers exactly at the region ends; optional Golden Gate tails (enzyme + fusion sites)
              with overlap trimming (e.g. AATG fusion site + a CDS starting ATG).
- pcr:        primer3 pairs flanking a region.
- sequencing: collection primers that already read into the region (filtered, scored), plus new
              forward/reverse primers when no collection primer qualifies.

Every primer shown is checked for other binding sites on this plasmid (3'-anchored matcher, the
same one used for the collection), and new ones for a collection primer already at that spot.
The simulated product is returned so the user can check it, and it is scanned for internal
Golden Gate sites. Positions are 1-based map coordinates.
Defaults follow a postdoc review of common bench practice.
"""

from __future__ import annotations

from dataclasses import dataclass

import primer3

from ..formats import MapRecord
from ..oligos import Oligo
from ..sequtil import revcomp
from . import thermo
from .primers import OligoGroup, PrimerHit, match_map

READ_MIN, READ_MAX, READ_AMBER, READ_TARGET = (
    50,
    800,
    700,
    150,
)  # sequencing primer 3' end -> region
SCAN_ENZYMES = {
    "BsaI": "GGTCTC",
    "BpiI": "GAAGAC",
    "Esp3I": "CGTCTC",
    "SapI": "GCTCTTC",
}

# UNVERIFIED (from a review agent, not re-checked against the papers): Golden Gate tail prefixes,
# 5' padding + recognition site + spacer; the 4-nt fusion site follows.
GG_PREFIX = {"BpiI": "TTGAAGACAA", "BsaI": "TTGGTCTCA"}
# PhytoBrick common-syntax core fusion sites (5' site, 3' site), Patron et al. 2015,
# New Phytol. doi:10.1111/nph.13532. Extend from the paper's full table when needed.
FUSION_SITES = {
    "promoter + 5'UTR": ("GGAG", "AATG"),
    "promoter": ("GGAG", "TACT"),
    "5'UTR": ("TACT", "AATG"),
    "CDS (with stop)": ("AATG", "GCTT"),
    "terminator": ("GCTT", "CGCT"),
}


@dataclass
class Window:
    """A rotated stretch of the map: window index i -> map position (offset + i) % L (0-based)."""

    seq: str
    offset: int
    L: int

    def to_map(self, i: int) -> int:
        return (self.offset + i) % self.L + 1


def _window(ref: MapRecord, start: int, end: int, pad: int) -> tuple[Window, int, int]:
    L = ref.length
    length = end - start + 1 if end >= start else end + L - start + 1
    if ref.circular:
        pad = max(0, min(pad, (L - length) // 2))
        off = (start - 1 - pad) % L
        s2 = ref.seq + ref.seq
        return Window(s2[off : off + length + 2 * pad], off, L), pad, length
    lo = max(0, start - 1 - pad)
    hi = min(L, end + pad)
    return Window(ref.seq[lo:hi], lo, L), start - 1 - lo, length


def _settings(
    opt_tm: float, product: tuple[int, int] | None = None, n: int = 5, size: tuple[int, int] | None = None
) -> dict:
    """primer3 settings. `size` = (min, max) length of the annealing part; default 18-27, optimum 21."""
    lo, hi = size or (18, 27)
    g = {
        **thermo.DESIGN_ARGS,
        "PRIMER_OPT_TM": opt_tm,
        "PRIMER_MIN_TM": opt_tm - 3,
        "PRIMER_MAX_TM": opt_tm + 3,
        "PRIMER_PAIR_MAX_DIFF_TM": 3.0,
        "PRIMER_OPT_SIZE": min(max(21, lo), hi),
        "PRIMER_MIN_SIZE": lo,
        "PRIMER_MAX_SIZE": hi,
        # GC 35-65 % is a warning in the UI, not a filter: AT-rich plant UTRs would return nothing.
        "PRIMER_MIN_GC": 20.0,
        "PRIMER_MAX_GC": 80.0,
        "PRIMER_GC_CLAMP": 1,
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_MAX_END_GC": 3,
        "PRIMER_MAX_SELF_ANY_TH": 45.0,
        "PRIMER_MAX_SELF_END_TH": 35.0,
        "PRIMER_MAX_HAIRPIN_TH": 47.0,
        "PRIMER_NUM_RETURN": n,
    }
    if product:
        g["PRIMER_PRODUCT_SIZE_RANGE"] = [list(product)]
    return g


def _oligo(seq: str, oid: str) -> Oligo:
    return Oligo(oid, oid, "", "", "", seq, seq, "", "")


def _check(
    ref: MapRecord, seq: str, strand: str, three: int, existing: list[PrimerHit]
) -> dict:
    hits = match_map(
        ref, [OligoGroup(seq, [_oligo(seq, "q")])], max_mismatches=2, min_binding=15
    )
    other = [h for h in hits if not (h.strand == strand and h.three_prime == three)]
    same = [
        h for h in existing if h.strand == strand and abs(h.three_prime - three) <= 3
    ]
    return {
        "other_sites": [
            {
                "strand": h.strand,
                "start": h.start,
                "end": h.end,
                "mismatches": h.mismatches,
            }
            for h in other
        ],
        "in_collection": [
            {"ids": h.oligo_ids, "names": h.names, "three_prime": h.three_prime}
            for h in same
        ],
    }


def _primer(
    seq: str,
    strand: str,
    start: int,
    end: int,
    three: int,
    tail: str,
    ref,
    existing,
    penalty: float | None = None,
) -> dict:
    tail = tail.upper()
    full = tail + seq
    gc = round(100 * sum(c in "GC" for c in seq) / len(seq), 1)
    return {
        "sequence": seq,
        "full": full,
        "tail": tail,
        "strand": strand,
        "start": start,
        "end": end,
        "three_prime": three,
        "length": len(seq),
        "tm": thermo.tm(seq),
        "tm_full": thermo.tm(full) if tail else None,
        "gc": gc,
        "gc_warning": not 35 <= gc <= 65,
        "hairpin_tm": thermo.hairpin_tm(full),
        "self_dimer_tm": thermo.self_dimer_tm(full),
        "penalty": None if penalty is None else round(penalty, 2),
        **_check(ref, seq, strand, three, existing),
    }


def _place(
    L: int, strand: str, start: int, end: int, trim: int, length: int
) -> tuple[int, int, int]:
    """Map coordinates (start, end, 3' end) of a binding part anchored at a region end."""
    if strand == "+":
        s1 = (start - 1 + trim) % L + 1
        e1 = (s1 - 1 + length - 1) % L + 1
        return s1, e1, e1
    e1 = (end - 1 - trim) % L + 1
    s1 = (e1 - length) % L + 1
    return s1, e1, s1


def _from_p3(p: dict, win: Window, side: str, tail: str, ref, existing) -> dict:
    pos, ln = p["COORDS"]
    if side == "left":
        s, e = win.to_map(pos), win.to_map(pos + ln - 1)
        return _primer(p["SEQUENCE"], "+", s, e, e, tail, ref, existing, p["PENALTY"])
    # primer3 gives the right primer's 5' end (rightmost base) and its length
    s, e = win.to_map(pos - ln + 1), win.to_map(pos)
    return _primer(p["SEQUENCE"], "-", s, e, s, tail, ref, existing, p["PENALTY"])


def _product(ref: MapRecord, fwd: dict, rev: dict) -> dict:
    """Simulated product, 5'->3' along the forward primer, and the Golden Gate sites in it."""
    L = ref.length
    s2 = ref.seq + ref.seq
    lo, hi = (fwd, rev) if fwd["strand"] == "+" else (rev, fwd)
    a, b = lo["start"] - 1, hi["end"]
    span = s2[a:b] if b > a else s2[a : b + L]
    core = span if fwd["strand"] == "+" else revcomp(span)
    seq = fwd["tail"] + core + revcomp(rev["tail"])
    sites = []
    for enz, site in SCAN_ENZYMES.items():
        for strand, motif in (("+", site), ("-", revcomp(site))):
            i = seq.find(motif)
            while i != -1:
                in_tail = i < len(fwd["tail"]) or i + len(motif) > len(seq) - len(
                    rev["tail"]
                )
                sites.append(
                    {
                        "enzyme": enz,
                        "position": i + 1,
                        "strand": strand,
                        "from_tail": in_tail,
                    }
                )
                i = seq.find(motif, i + 1)
    sites.sort(key=lambda d: d["position"])
    return {
        "product": len(seq),
        "product_seq": seq,
        "product_head": seq[:30],
        "product_tail": seq[-30:],
        "internal_sites": [d for d in sites if not d["from_tail"]],
        "tail_sites": [d for d in sites if d["from_tail"]],
    }


def _scope(ref: MapRecord) -> str:
    return f"{ref.plasmid_id or 'this map'} ({ref.length:,} bp) only"


def gg_tails(enzyme: str, oh5: str, oh3: str) -> tuple[str, str]:
    """Forward and reverse 5' tails for a Golden Gate part. oh3 is written on the part's sense strand."""
    if enzyme not in GG_PREFIX:
        raise ValueError(f"enzyme must be one of {', '.join(GG_PREFIX)}")
    oh5, oh3 = oh5.upper(), oh3.upper()
    if len(oh5) != 4 or len(oh3) != 4 or set(oh5 + oh3) - set("ACGT"):
        raise ValueError("fusion sites must be 4 nt of A/C/G/T")
    return GG_PREFIX[enzyme] + oh5, GG_PREFIX[enzyme] + revcomp(oh3)


def _overlap(oh: str, seq: str) -> int:
    """Leading bases of seq already supplied by the end of the overhang (AATG + ATG... -> 3).

    Only 3 or 4 shared bases count: a 1-2 nt chance overlap (AATG + G...) would silently move a
    CDS start into the fusion site and shift its frame.
    """
    for k in range(min(4, len(seq)), 2, -1):
        if oh[-k:] == seq[:k]:
            return k
    return 0


def design_cloning(
    ref: MapRecord,
    start: int,
    end: int,
    existing: list[PrimerHit],
    opt_tm: float = 60.0,
    tail_f: str = "",
    tail_r: str = "",
    sense: str = "+",
    enzyme: str = "",
    oh5: str = "",
    oh3: str = "",
    size: tuple[int, int] | None = None,
) -> dict:
    """Primers that start exactly at the region ends (clone a CDS/part as-is; tails go 5').

    `sense` is the strand the part reads on: for a minus-strand CDS the forward primer sits at
    the region's right end. With `enzyme` + `oh5`/`oh3`, Golden Gate tails are built here and
    region bases already supplied by a fusion site are trimmed from the binding part.
    """
    win, ri, length = _window(ref, start, end, pad=0)
    if length < 40:
        raise ValueError("region must be at least 40 bp")
    top = win.seq[ri : ri + length]
    part = top if sense == "+" else revcomp(top)  # the part on its own sense strand
    trim5 = trim3 = 0
    if enzyme:
        tail_f, tail_r = gg_tails(enzyme, oh5, oh3)
        trim5 = _overlap(oh5.upper(), part)
        trim3 = _overlap(revcomp(oh3.upper()), revcomp(part))
    other = "-" if sense == "+" else "+"

    def make(src: str, strand: str, tail: str, trim: int) -> dict:
        src = src[trim:]
        lo, hi = size or (18, 35)
        options = [src[:n] for n in range(lo, min(hi, len(src)) + 1)] or [src]
        b = next((o for o in options if (thermo.tm(o) or 0) >= opt_tm), None)
        if b is None:
            b = min(options, key=lambda o: abs((thermo.tm(o) or 0) - opt_tm))
        s1, e1, three = _place(ref.length, strand, start, end, trim, len(b))
        return _primer(b, strand, s1, e1, three, tail, ref, existing)

    fwd = make(part, sense, tail_f, trim5)
    rev = make(revcomp(part), other, tail_r, trim3)
    pair = {
        "forward": fwd,
        "reverse": rev,
        "pair_dimer_tm": thermo.hetero_dimer_tm(fwd["full"], rev["full"]),
        "tm_difference": round(abs((fwd["tm"] or 0) - (rev["tm"] or 0)), 1),
        **_product(ref, fwd, rev),
    }
    warnings = []
    if enzyme and oh5.upper() == "AATG" and not part.startswith("ATG"):
        warnings.append(
            "part does not start with ATG: the AATG fusion site adds an ATG in front of it; "
            "check the reading frame, or use a different 5' fusion site"
        )
    if pair["internal_sites"]:
        enz = ", ".join(sorted({d["enzyme"] for d in pair["internal_sites"]}))
        warnings.append(
            f"internal {enz} site(s) in the product: domesticate before Golden Gate assembly"
        )
    return {
        "task": "cloning",
        "region": [start, end],
        "sense": sense,
        "pairs": [pair],
        "warnings": warnings,
        "trimmed": {"5prime": trim5, "3prime": trim3},
        "tm_model": thermo.MODEL,
        "checked_against": _scope(ref),
    }


def design_pcr(
    ref: MapRecord,
    start: int,
    end: int,
    existing: list[PrimerHit],
    opt_tm: float = 60.0,
    max_extra: int = 300,
    tail_f: str = "",
    tail_r: str = "",
    n: int = 5,
    size: tuple[int, int] | None = None,
) -> dict:
    """Pairs whose product covers [start, end] with at most max_extra bases of flank on each side."""
    win, ri, length = _window(ref, start, end, pad=max_extra + 30)
    r = primer3.design_primers(
        {
            "SEQUENCE_ID": ref.plasmid_id or "map",
            "SEQUENCE_TEMPLATE": win.seq,
            "SEQUENCE_TARGET": [ri, length],
        },
        _settings(opt_tm, (length, length + 2 * max_extra), n, size),
    )
    pairs = []
    for i, pair in enumerate(r.get("PRIMER_PAIR", [])):
        f = _from_p3(r["PRIMER_LEFT"][i], win, "left", tail_f, ref, existing)
        rv = _from_p3(r["PRIMER_RIGHT"][i], win, "right", tail_r, ref, existing)
        pairs.append(
            {
                "forward": f,
                "reverse": rv,
                "penalty": round(pair["PENALTY"], 2),
                "pair_dimer_tm": thermo.hetero_dimer_tm(f["full"], rv["full"]),
                "tm_difference": round(abs((f["tm"] or 0) - (rv["tm"] or 0)), 1),
                **_product(ref, f, rv),
            }
        )
    return {
        "task": "pcr",
        "region": [start, end],
        "pairs": pairs,
        "tm_model": thermo.MODEL,
        "checked_against": _scope(ref),
        "explain": {
            k: r.get(k, "")
            for k in (
                "PRIMER_LEFT_EXPLAIN",
                "PRIMER_RIGHT_EXPLAIN",
                "PRIMER_PAIR_EXPLAIN",
            )
        },
    }


def design_sequencing(
    ref: MapRecord,
    start: int,
    end: int,
    existing: list[PrimerHit],
    opt_tm: float = 58.0,
    n: int = 3,
    size: tuple[int, int] | None = None,
) -> dict:
    """Collection primers reading into the region (filtered, scored), plus new ones where none qualify."""
    L = ref.length
    sites_per_oligo: dict[str, int] = {}
    for h in existing:
        sites_per_oligo[h.oligo_ids[0]] = sites_per_oligo.get(h.oligo_ids[0], 0) + 1
    fwd, rev = [], []
    for h in existing:
        if h.strand == "+":
            d = (start - h.three_prime) % L if ref.circular else start - h.three_prime
        else:
            d = (h.three_prime - end) % L if ref.circular else h.three_prime - end
        if not READ_MIN <= d <= READ_MAX:
            continue
        if sites_per_oligo[h.oligo_ids[0]] != 1 or h.tm is None or h.tm < 50:
            continue
        # one transparent score, lower is better: distance from ideal + mismatches + late-start penalty
        score = round(
            abs(d - READ_TARGET) / 50 + 2 * h.mismatches + (2 if d > READ_AMBER else 0),
            1,
        )
        entry = {
            "ids": h.oligo_ids,
            "names": h.names,
            "three_prime": h.three_prime,
            "distance": d,
            "tm": h.tm,
            "mismatches": h.mismatches,
            "score": score,
            "late": d > READ_AMBER,
        }
        (fwd if h.strand == "+" else rev).append(entry)
    fwd.sort(key=lambda x: x["score"])
    rev.sort(key=lambda x: x["score"])

    new = {"forward": [], "reverse": []}
    win, ri, length = _window(ref, start, end, pad=READ_MAX + 30)
    if not fwd:
        a = max(0, ri - 300)
        r = primer3.design_primers(
            {
                "SEQUENCE_TEMPLATE": win.seq,
                "SEQUENCE_INCLUDED_REGION": [a, ri - READ_MIN - a],
            },
            {**_settings(opt_tm, n=n, size=size), "PRIMER_PICK_RIGHT_PRIMER": 0},
        )
        new["forward"] = [
            _from_p3(p, win, "left", "", ref, existing)
            for p in r.get("PRIMER_LEFT", [])
        ]
    if not rev:
        b = ri + length + READ_MIN
        r = primer3.design_primers(
            {
                "SEQUENCE_TEMPLATE": win.seq,
                "SEQUENCE_INCLUDED_REGION": [b, min(250, len(win.seq) - b)],
            },
            {**_settings(opt_tm, n=n, size=size), "PRIMER_PICK_LEFT_PRIMER": 0},
        )
        new["reverse"] = [
            _from_p3(p, win, "right", "", ref, existing)
            for p in r.get("PRIMER_RIGHT", [])
        ]
    note = None
    if length > READ_AMBER:
        note = (
            f"region is {length:,} bp; a Sanger read is good for ~50–700 bp past the primer, so this "
            f"needs ~{-(-length // 550)} reads (tile every 500–600 bp). Above ~1.5 kb a whole-plasmid run may be cheaper."
        )
    return {
        "task": "sequencing",
        "region": [start, end],
        "existing_forward": fwd[:8],
        "existing_reverse": rev[:8],
        "new": new,
        "note": note,
        "tm_model": thermo.MODEL,
        "checked_against": _scope(ref),
        "filters": "collection primers shown only if: single site on this map, Tm ≥ 50 °C, 3′ end 50–800 bp from the region",
    }


def design_walk(ref: MapRecord, start: int, end: int, existing: list[PrimerHit], step: int = 550,
                opt_tm: float = 58.0, both_strands: bool = False, size: tuple[int, int] | None = None) -> dict:
    """Primer walk for a long region: forward tiles every `step` bp (collection primers first),
    plus one reverse primer anchored past the end so the last stretch is read both ways.
    With both_strands, a reverse tile is placed for every forward one as well."""
    L = ref.length
    length = end - start + 1 if end >= start else end + L - start + 1
    if step < 300 or step > 800:
        raise ValueError("step must be 300-800 bp")
    sites_per_oligo: dict[str, int] = {}
    for h in existing:
        sites_per_oligo[h.oligo_ids[0]] = sites_per_oligo.get(h.oligo_ids[0], 0) + 1

    def usable(h: PrimerHit) -> bool:
        return sites_per_oligo[h.oligo_ids[0]] == 1 and h.tm is not None and h.tm >= 50

    def pos(offset: int) -> int:  # 1-based map position `offset` bp after region start (may be negative)
        return (start - 1 + offset) % L + 1 if ref.circular else start + offset

    tiles = []
    targets = list(range(-READ_TARGET, length - 200, step))  # 3' ends: first one ~150 bp before the region
    for t in targets:
        tiles.append(_tile(ref, existing, usable, "+", pos(t), opt_tm, t, size))
        if both_strands:
            tiles.append(_tile(ref, existing, usable, "-", pos(min(t + step + READ_TARGET, length + READ_TARGET)), opt_tm, t, size))
    if not both_strands:
        tiles.append(_tile(ref, existing, usable, "-", pos(length + READ_TARGET - 1), opt_tm, length + READ_TARGET, size))
    reads = len(tiles)
    covered = sorted({t["covers"] for t in tiles})
    return {"task": "walk", "region": [start, end], "length": length, "step": step, "tiles": tiles,
            "reads": reads, "tm_model": thermo.MODEL, "checked_against": _scope(ref),
            "note": (f"{reads} Sanger reads for {length:,} bp. " + (
                "A whole-plasmid run may be cheaper." if length > 1500 else "")).strip(), "covered": covered}


def _tile(ref, existing, usable, strand: str, target: int, opt_tm: float, offset: int,
          size: tuple[int, int] | None = None) -> dict:
    """One walking primer whose 3' end sits near `target` (1-based); collection primer if one fits."""
    L = ref.length

    def dist(a: int, b: int) -> int:
        d = abs(a - b)
        return min(d, L - d) if ref.circular else d

    near = [h for h in existing if h.strand == strand and usable(h) and dist(h.three_prime, target) <= 80]
    near.sort(key=lambda h: (h.mismatches, dist(h.three_prime, target)))
    base = {"strand": strand, "target": target, "covers": offset}
    if near:
        h = near[0]
        return {**base, "source": "collection", "ids": h.oligo_ids, "names": h.names, "three_prime": h.three_prime,
                "tm": h.tm, "mismatches": h.mismatches, "offset_from_target": h.three_prime - target}
    # design: a window of ~160 bp ending at the target (forward) or starting at it (reverse)
    w_start = target - 160 if strand == "+" else target
    win, ri, wlen = _window(ref, max(1, w_start) if not ref.circular else (w_start - 1) % L + 1,
                            min(L, w_start + 159) if not ref.circular else (w_start + 158) % L + 1, pad=0)
    side = "LEFT" if strand == "+" else "RIGHT"
    r = primer3.design_primers({"SEQUENCE_TEMPLATE": win.seq, "SEQUENCE_INCLUDED_REGION": [0, len(win.seq)]},
                               {**_settings(opt_tm, n=1, size=size), f"PRIMER_PICK_{'RIGHT' if strand == '+' else 'LEFT'}_PRIMER": 0})
    picks = r.get(f"PRIMER_{side}", [])
    if not picks:
        return {**base, "source": "none", "explain": r.get(f"PRIMER_{side}_EXPLAIN", "")}
    p = _from_p3(picks[0], win, side.lower(), "", ref, existing)
    return {**base, "source": "new", **p}
