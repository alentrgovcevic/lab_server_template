"""One place for Tm conditions, so every number the app shows uses the same model.

Q5/Phusion-like reaction: 50 mM monovalent, 2.0 mM Mg2+, 0.2 mM dNTP,
500 nM primer. primer3's own default is 1.5 mM Mg2+, 0.6 mM dNTP, 50 nM primer.
"""

from __future__ import annotations

import primer3

CONDITIONS = {"mv_conc": 50.0, "dv_conc": 2.0, "dntp_conc": 0.2, "dna_conc": 500.0}
MODEL = "primer3 nearest-neighbour (SantaLucia 1998), 50 mM Na⁺, 2.0 mM Mg²⁺, 0.2 mM dNTP, 500 nM primer"

# The same conditions under primer3's design_primers() argument names.
DESIGN_ARGS = {
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 2.0,
    "PRIMER_DNTP_CONC": 0.2,
    "PRIMER_DNA_CONC": 500.0,
}


def tm(seq: str) -> float | None:
    if not seq or set(seq.upper()) - set("ACGT"):
        return None
    return round(primer3.calc_tm(seq.upper(), **CONDITIONS), 1)


def hairpin_tm(seq: str) -> float:
    r = primer3.calc_hairpin(seq.upper(), **CONDITIONS)
    return max(0.0, round(r.tm, 1)) if r.structure_found else 0.0


def self_dimer_tm(seq: str) -> float:
    r = primer3.calc_homodimer(seq.upper(), **CONDITIONS)
    return max(0.0, round(r.tm, 1)) if r.structure_found else 0.0


def hetero_dimer_tm(a: str, b: str) -> float:
    r = primer3.calc_heterodimer(a.upper(), b.upper(), **CONDITIONS)
    return max(0.0, round(r.tm, 1)) if r.structure_found else 0.0
