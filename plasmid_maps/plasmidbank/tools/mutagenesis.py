"""Site-directed mutagenesis primers, Q5 SDM style (back-to-back, non-overlapping).

The whole plasmid is amplified by two primers whose 5' ends meet at the change, and the linear
product is phosphorylated, ligated and template-digested (KLD). Here:
- the forward primer anneals from the first base after the change; the new bases are its 5' tail;
- the reverse primer anneals from the last base before the change, pointing the other way;
- an insertion longer than SPLIT_ABOVE nt is split: the second half on the forward primer, the
  first half (reverse complemented) on the reverse primer;
- a deletion needs no tail: the primers simply flank the deleted bases.

A mutation is given either as DNA (replace map bases start..end with new; end = start - 1 means a
pure insertion before start) or as a protein change in a CDS feature: K48R, Y112*, K48del,
K48_L50del, K48_L49insGS. New codons follow this CDS's own codon usage (the most used codon for
that amino acid, and for a substitution the fewest base changes first). Positions are 1-based.
"""

from __future__ import annotations

import re
from collections import Counter

from Bio.Seq import Seq

from ..formats import Feature, MapRecord
from ..sequtil import revcomp
from . import thermo
from .design import _primer, _scope

# UNVERIFIED: NEB's guidance puts short insertions wholly in the forward primer and splits longer
# ones; the exact cut-off was not re-checked. Change here if the lab prefers another.
SPLIT_ABOVE = 6
MAX_CHANGE = 120  # inserted bases; longer belongs in a Golden Gate or Gibson design

_CODONS: dict[str, list[str]] = {}
for _c in (a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"):
    _CODONS.setdefault(str(Seq(_c).translate()), []).append(_c)


# ── protein changes -> DNA changes ────────────────────────────────────────────


def cds_positions(ref: MapRecord, f: Feature) -> list[int]:
    """0-based map positions of the CDS bases, in coding order (5'->3' on its own strand)."""
    L = ref.length
    n = (f.end - f.start) % L or (L if f.end == f.start else 0)
    pos = [(f.start + i) % L for i in range(n)]
    return pos[::-1] if f.strand < 0 else pos


def _coding(ref: MapRecord, pos: list[int], minus: bool) -> str:
    s = "".join(ref.seq[p] for p in pos).upper()
    return s.translate(str.maketrans("ACGT", "TGCA")) if minus else s


def _usage(cds: str) -> Counter:
    return Counter(cds[i : i + 3] for i in range(0, len(cds) - 2, 3))


def _codon_for(aa: str, usage: Counter, old: str | None = None) -> str:
    opts = _CODONS.get(aa)
    if not opts:
        raise ValueError(f"unknown amino acid {aa!r}")
    diff = (lambda c: sum(x != y for x, y in zip(c, old))) if old else (lambda c: 0)
    return min(opts, key=lambda c: (diff(c), -usage[c], c))


_MUT = re.compile(
    r"^(?P<a1>[A-Z*])(?P<n1>\d+)(?:_(?P<a2>[A-Z*])(?P<n2>\d+))?(?P<op>DEL|INS(?P<ins>[A-Z*]+)|(?P<to>[A-Z*]))$"
)


def protein_change(ref: MapRecord, f: Feature, mutation: str) -> dict:
    """Turn K48R-style notation into a DNA change on the map (start, end, new; 1-based)."""
    m = _MUT.match("".join(mutation.split()).upper())
    if not m:
        raise ValueError(
            "write the change as K48R, Y112*, K48del, K48_L50del or K48_L49insGS"
        )
    pos = cds_positions(ref, f)
    minus = f.strand < 0
    cds = _coding(ref, pos, minus)
    prot = str(Seq(cds[: len(cds) - len(cds) % 3]).translate())
    usage = _usage(cds)
    checks = []
    if f.type != "CDS":
        checks.append(f"{f.label} is annotated as {f.type}, not CDS: check that it starts at the start codon")
    if len(cds) % 3:
        checks.append(f"{f.label} is {len(cds)} bp, not a multiple of 3")
    if "*" in prot[:-1]:
        checks.append(f"{f.label} has a stop codon inside it in this frame")
    if not prot.startswith("M"):
        checks.append(f"{f.label} does not start with M (it starts {prot[:6]}…)")
    other = str(Seq(revcomp(cds)[: len(cds) - len(cds) % 3]).translate())
    if "*" not in other[:-1] and ("*" in prot[:-1] or not prot.startswith("M")):
        checks.append(f"the other strand reads through without a stop ({other[:10]}…): the feature's strand may be annotated wrongly")
    n1 = int(m["n1"])
    n2 = int(m["n2"]) if m["n2"] else n1
    for a, n in ((m["a1"], n1), (m["a2"], n2)):
        if a is None:
            continue
        if not 1 <= n <= len(prot):
            raise ValueError(f"{f.label} has {len(prot)} codons; there is no codon {n}")
        if prot[n - 1] != a:
            raise ValueError(f"codon {n} of {f.label} is {prot[n - 1]}, not {a}")
    if n2 < n1:
        raise ValueError("the range must run from the lower to the higher codon")

    def span(
        c1: int, c2: int
    ) -> list[int]:  # map positions of codons c1..c2, ascending
        return sorted(pos[3 * (c1 - 1) : 3 * c2])

    op = m["op"]
    if op == "DEL":
        ps = span(n1, n2)
        if max(ps) - min(ps) != len(ps) - 1:
            raise ValueError("the deletion crosses the origin; rotate the map first")
        old, new_dna, new_prot = (
            "".join(ref.seq[p] for p in ps),
            "",
            prot[: n1 - 1] + prot[n2:],
        )
        start, end = min(ps) + 1, max(ps) + 1
        codons = []
    elif m["ins"]:
        if n2 != n1 + 1:
            raise ValueError(
                "an insertion goes between two neighbouring codons: K48_L49insGS"
            )
        codons = [_codon_for(a, usage) for a in m["ins"]]
        coding = "".join(codons)
        last = pos[3 * n1 - 1]  # last base of codon n1
        if (
            minus
        ):  # the coding strand runs right to left: insert before that base on the map
            start, end, new_dna = last + 1, last, revcomp(coding)
        else:
            start, end, new_dna = last + 2, last + 1, coding
        old, new_prot = "", prot[:n1] + m["ins"] + prot[n1:]
    else:
        if m["n2"]:
            raise ValueError("a substitution names one codon: K48R")
        to = m["to"]
        old_codon = cds[3 * (n1 - 1) : 3 * n1]
        codon = _codon_for(to, usage, old_codon)
        codons = [codon]
        ps = span(n1, n1)
        new_plus = revcomp(codon) if minus else codon
        old_plus = "".join(ref.seq[p] for p in ps).upper()
        # only the bases that differ need to be in the primer
        i = next((k for k in range(3) if old_plus[k] != new_plus[k]), None)
        if i is None:
            raise ValueError(f"codon {n1} already encodes {to}")
        j = max(k for k in range(3) if old_plus[k] != new_plus[k])
        if ps[j] - ps[i] != j - i:
            raise ValueError("that codon crosses the origin; rotate the map first")
        start, end = ps[i] + 1, ps[j] + 1
        old, new_dna = old_plus[i : j + 1], new_plus[i : j + 1]
        new_prot = prot[: n1 - 1] + to + prot[n1:]
    return {
        "start": start,
        "end": end,
        "new": new_dna,
        "old": old,
        "feature": f.label,
        "notation": mutation.strip(),
        "codons": codons,
        "protein_before": prot,
        "protein_after": new_prot,
        "checks": checks,
        "codon_note": "codons are numbered from the feature's first base; new codons are this CDS's most used codon for the amino acid (for a substitution, fewest base changes first)",
    }


# ── primers ───────────────────────────────────────────────────────────────────


def _anneal(seq: str, opt_tm: float, size: tuple[int, int]) -> str:
    """Shortest prefix of seq (lo..hi nt) reaching opt_tm; else the one closest to it."""
    lo, hi = size
    options = [seq[:n] for n in range(lo, min(hi, len(seq)) + 1)]
    if not options:
        raise ValueError("not enough sequence next to the change for a primer")
    hit = next((o for o in options if (thermo.tm(o) or 0) >= opt_tm), None)
    return hit or min(options, key=lambda o: abs((thermo.tm(o) or 0) - opt_tm))


def design_sdm(
    ref: MapRecord,
    start: int,
    end: int,
    new: str,
    existing: list,
    opt_tm: float = 60.0,
    size: tuple[int, int] | None = None,
    protein: dict | None = None,
) -> dict:
    L, seq = ref.length, ref.seq.upper()
    new = "".join(new.upper().split())
    if set(new) - set("ACGT"):
        raise ValueError("the new sequence may only contain A, C, G, T")
    if not (1 <= start <= L + 1 and start - 1 <= end <= L):
        raise ValueError(
            f"positions must be within 1..{L}, with end ≥ start - 1 (end = start - 1: insertion)"
        )
    if len(new) > MAX_CHANGE:
        raise ValueError(
            f"insertions above {MAX_CHANGE} nt belong in an assembly design (Golden Gate or Gibson)"
        )
    old = seq[start - 1 : end]
    if old == new:
        raise ValueError("the new sequence is the same as the map")
    size = size or (18, 35)
    hi = size[1]
    if ref.circular:
        s2 = seq + seq
        f_from = end % L  # 0-based index of the first base after the change
        down = s2[f_from : f_from + hi]
        up = s2[start - 1 + L - hi : start - 1 + L]  # the bases ending just before the change
    else:
        f_from = end
        down = seq[end : end + hi]
        up = seq[max(0, start - 1 - hi) : start - 1]
    if len(down) < size[0] or len(up) < size[0]:
        raise ValueError("the change is too close to the end of this linear sequence for back-to-back primers")
    fa = _anneal(down, opt_tm, size)  # forward: anneals from the first base after the change, on +
    ra = _anneal(revcomp(up), opt_tm, size)  # reverse: from the last base before the change, on -
    if len(new) - len(old) > SPLIT_ABOVE:  # a long insertion: half on each primer's 5' end
        h = len(new) // 2
        f_tail, r_tail = new[h:], revcomp(new[:h])
    else:
        f_tail, r_tail = new, ""
    fs = f_from % L + 1
    fe = (fs - 1 + len(fa) - 1) % L + 1
    re_ = (start - 2) % L + 1 if ref.circular else start - 1
    rs = (re_ - len(ra)) % L + 1
    fwd = _primer(fa, "+", fs, fe, fe, f_tail, ref, existing)
    rev = _primer(ra, "-", rs, re_, rs, r_tail, ref, existing)
    mutant = seq[: start - 1] + new + seq[end:]
    k = start - 1
    kind = (
        "deletion"
        if not new
        else "insertion"
        if not old
        else "substitution"
        if len(old) == len(new)
        else "replacement"
    )
    warnings = []
    for p, n in ((fwd, "forward"), (rev, "reverse")):
        if len(p["full"]) > 60:
            warnings.append(f"{n} primer is {len(p['full'])} nt: order HPLC-purified")
        if p["other_sites"]:
            warnings.append(f"{n} primer also binds elsewhere on this plasmid")
    tms = [t for t in (fwd["tm"], rev["tm"]) if t]
    return {
        "task": "sdm",
        "kind": kind,
        "region": [start, max(start, end)],
        "change": {"start": start, "end": end, "old": old, "new": new},
        "protein": protein,
        "forward": fwd,
        "reverse": rev,
        "tm_difference": round(abs(tms[0] - tms[1]), 1) if len(tms) == 2 else None,
        "pair_dimer_tm": thermo.hetero_dimer_tm(fwd["full"], rev["full"]),
        "context_before": seq[max(0, k - 15) : end + 15],
        "context_after": mutant[max(0, k - 15) : k + len(new) + 15],
        "mutant_length": len(mutant),
        "mutant_seq": mutant,
        "warnings": warnings,
        "checked_against": _scope(ref),
        "tm_model": thermo.MODEL,
    }
