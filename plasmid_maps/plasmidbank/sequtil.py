"""Sequence helpers shared by the tools: IUPAC sets, reverse complement, canonical hash."""

from __future__ import annotations

import hashlib

IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}
IUPAC_CHARS = frozenset(IUPAC)
_COMP = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn"
)


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def compatible(a: str, b: str) -> bool:
    """True when two IUPAC codes can denote the same base (N matches anything)."""
    sa, sb = IUPAC.get(a.upper()), IUPAC.get(b.upper())
    if sa is None or sb is None:
        return False
    return bool(set(sa) & set(sb))


def _least_rotation(s: str) -> int:
    """Booth's algorithm: start index of the lexicographically least rotation."""
    s2 = s + s
    f = [-1] * len(s2)
    k = 0
    for j in range(1, len(s2)):
        i = f[j - k - 1]
        while i != -1 and s2[j] != s2[k + i + 1]:
            if s2[j] < s2[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if s2[j] != s2[k + i + 1]:
            if s2[j] < s2[k]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return k


def canonical_seq(seq: str, circular: bool) -> str:
    """Same string for the same molecule regardless of start position or strand."""
    seq = seq.upper()
    rc = revcomp(seq)
    if not circular:
        return min(seq, rc)
    rots = []
    for x in (seq, rc):
        k = _least_rotation(x)
        rots.append(x[k:] + x[:k])
    return min(rots)


def canonical_hash(seq: str, circular: bool) -> str:
    return hashlib.sha1(canonical_seq(seq, circular).encode()).hexdigest()
