"""Which oligos in the lab collection bind a map, where, and how well.

Matching is anchored on the oligo's 3' end, because that is what primes:
- the 3'-terminal SEED nt must match exactly (IUPAC-compatible) - a 3' mismatch kills priming;
- binding then extends towards the 5' end while mismatches stay <= max_mismatches;
- anything 5' of the binding part is a tail (Golden Gate sites, overhangs, tags), reported, not penalised.
Both strands are searched, and across the origin for circular maps.
Tm is nearest-neighbour (SantaLucia) via primer3 for the binding part only, under the reaction
conditions in thermo.py.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from ..formats import MapRecord
from ..oligos import Oligo
from ..sequtil import compatible, revcomp
from .thermo import tm

SEED = 12


@dataclass
class OligoGroup:
    """Oligos with an identical normalised sequence, matched once."""

    seq: str
    oligos: list[Oligo] = field(default_factory=list)


@dataclass
class PrimerHit:
    oligo_ids: list[str]
    names: list[str]
    strand: str  # "+" binds/extends along the top strand (forward primer), "-" reverse
    start: int  # 1-based first map base of the binding part (top-strand coordinates)
    end: int  # 1-based last map base of the binding part
    three_prime: int  # 1-based map position of the primer's 3' end
    binding_length: int
    mismatches: int
    tail_length: int
    tm: float | None
    oligo_length: int


def group_oligos(oligos: list[Oligo]) -> tuple[list[OligoGroup], dict[str, int]]:
    """Group by sequence. Returns the groups plus counts of oligos that cannot be matched."""
    groups: dict[str, OligoGroup] = {}
    skipped: dict[str, int] = defaultdict(int)
    for o in oligos:
        if not o.seq:
            skipped["no sequence"] += 1
            continue
        if len(o.seq) < SEED + 3:
            skipped["shorter than 15 nt"] += 1
            continue
        if set(o.seq[-SEED:]) - set("ACGT"):
            skipped["degenerate 3' end"] += 1
            continue
        groups.setdefault(o.seq, OligoGroup(o.seq)).oligos.append(o)
    return list(groups.values()), dict(skipped)


def _seed_index(s: str) -> dict[str, list[int]]:
    idx: dict[str, list[int]] = defaultdict(list)
    for i in range(len(s) - SEED + 1):
        idx[s[i : i + SEED]].append(i)
    return idx


def _tm(seq: str) -> float | None:
    return tm(seq)


def match_map(
    ref: MapRecord,
    groups: list[OligoGroup],
    max_mismatches: int = 2,
    min_binding: int = 15,
) -> list[PrimerHit]:
    L = ref.length
    if L < SEED:
        return []
    top = ref.seq + (ref.seq[:200] if ref.circular else "")
    bottom = revcomp(ref.seq)
    bottom = bottom + (bottom[:200] if ref.circular else "")
    idx = {"+": _seed_index(top), "-": _seed_index(bottom)}
    strings = {"+": top, "-": bottom}
    hits: dict[tuple, PrimerHit] = {}
    for g in groups:
        o = g.seq
        seed = o[-SEED:]
        for strand in ("+", "-"):
            s = strings[strand]
            for p in idx[strand].get(seed, ()):
                # extend 5'-wards from the seed
                mm, k, best_len, best_mm = 0, 0, SEED, 0
                while k < len(o) - SEED and p - 1 - k >= 0:
                    if not compatible(o[-SEED - 1 - k], s[p - 1 - k]):
                        mm += 1
                        if mm > max_mismatches:
                            break
                    k += 1
                    if compatible(o[-SEED - k], s[p - k]):
                        best_len, best_mm = SEED + k, mm
                if best_len < min_binding:
                    continue
                b_start = p + SEED - best_len  # on this strand's string, 0-based
                b_end = p + SEED  # exclusive
                if strand == "+":
                    start0, end0 = b_start % L, (b_end - 1) % L
                    three = end0
                else:
                    start0, end0 = (L - b_end) % L, (L - 1 - b_start) % L
                    three = start0
                key = (strand, three, o)
                if key in hits:
                    continue
                hits[key] = PrimerHit(
                    [x.oligo_id for x in g.oligos],
                    [x.name for x in g.oligos],
                    strand,
                    start0 + 1,
                    end0 + 1,
                    three + 1,
                    best_len,
                    best_mm,
                    len(o) - best_len,
                    _tm(s[b_start:b_end]),
                    len(o),
                )
    return sorted(hits.values(), key=lambda h: (h.three_prime, h.strand))
