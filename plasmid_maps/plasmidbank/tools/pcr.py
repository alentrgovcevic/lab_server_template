"""A PCR from two primers on a map: the product, and a cycling program for a chosen polymerase.

"select two primers and get the PCR product using selected polymerases, instead of
going to Thermo Fisher and calculating yourself".

The rules below are the vendors' own, retrieved 2026-10-27 by a web-research sub-agent. Every value
carries its source and how sure we are of it (`confidence`): NEB's site refused direct fetches, so
the NEB values come from NEB PDFs hosted by a reseller (MEDIUM) or from search summaries (LOW).

🛑 The Tm here is this app's model (`thermo.MODEL`), not the vendor's calculator: NEB and Thermo each
use their own method and salt, so their numbers can differ by a few °C. The app says so next to every
annealing temperature; the vendor calculator stays the reference for a critical PCR.
"""

from __future__ import annotations

import math

from ..formats import MapRecord
from ..oligos import Oligo
from . import thermo
from ..sequtil import revcomp
from .primers import group_oligos, match_map

NEB_Q5 = "https://www.neb.com/protocols/2013/12/13/pcr-using-q5-high-fidelity-dna-polymerase-m0491"
NEB_Q5_COPY = "https://genexpress.cl/wp-content/uploads/2025/01/pcr-using-q5-high-fidelity-dna-polymerase-m0491.pdf"
THERMO_PHUSION = "https://documents.thermofisher.com/TFS-Assets/BID/Reference-Materials/phusion-dna-polymerases-labaid.pdf"
THERMO_TA = (
    "https://www.thermofisher.com/global/en/home/brands/thermo-scientific/molecular-biology/molecular-biology-"
    "learning-center/molecular-biology-resource-library/optimizing-tm-and-annealing.html"
)
THERMO_SUPERFI2 = "https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0018859_Platinum_SuperFi_II_DNA_Pol_UG.pdf"
NEB_TAQ = (
    "https://www.neb.com/protocols/taq-dna-polymerase-with-standard-taq-buffer-m0273"
)
NEB_ONETAQ_COPY = (
    "https://mejoras.genexpress.cl/wp-content/uploads/2023/01/NEB.M0481-OneTaq.pdf"
)

# s per kb of product at the extension step, by template; "min" is the shortest step worth programming.
POLYMERASES = {
    "Q5": {
        "name": "NEB Q5 High-Fidelity",
        "ta": "lower_tm_plus_3",
        "ta_max": 72,
        "two_step_at": 72,
        "denature": (98, 30),
        "cycle_denature": (98, "5–10 s"),
        "anneal_s": "10–30 s",
        "extend_c": 72,
        "s_per_kb": {"plasmid": 10, "genomic": 30, "cdna": 40},
        "min_s": 10,
        "cycles": {"plasmid": "25–35", "genomic": "30–35", "cdna": "25–35"},
        "final": (72, "2 min"),
        "confidence": "MEDIUM (NEB's protocol PDF via a reseller copy; NEB's site refused direct reads)",
        "source": [NEB_Q5, NEB_Q5_COPY],
        "notes": [
            "NEB: anneal at 3 °C above the lower primer's Tm from the NEB Tm Calculator; "
            "with an annealing temperature ≥ 72 °C a 2-step program is possible.",
            "Extension 10 s/kb for plasmid (or any product < 1 kb), 20–30 s/kb genomic, 40 s/kb cDNA or long.",
        ],
    },
    "Phusion": {
        "name": "Thermo Phusion High-Fidelity",
        "ta": "phusion",
        "ta_max": 72,
        "two_step_at": 72,
        "denature": (98, 30),
        "cycle_denature": (98, "5–10 s"),
        "anneal_s": "10–30 s",
        "extend_c": 72,
        "s_per_kb": {"plasmid": 15, "genomic": 30, "cdna": 30},
        "min_s": 10,
        "cycles": {"plasmid": "25–35", "genomic": "25–35", "cdna": "25–35"},
        "final": (72, "5–10 min"),
        "confidence": "HIGH for cycling (Thermo LabAid 2018); MEDIUM for the annealing rule (Thermo web page, undated)",
        "source": [THERMO_PHUSION, THERMO_TA],
        "notes": [
            "Thermo: primers ≤ 20 nt anneal at the lower Tm; longer primers at the lower Tm + 3 °C "
            "(Tm from Thermo's calculator, which uses its own method).",
            "Extension 15–30 s/kb, never more than 1 min/kb.",
        ],
    },
    "SuperFi II": {
        "name": "Thermo Platinum SuperFi II",
        "ta": "fixed_60",
        "ta_max": 60,
        "two_step_at": None,
        "denature": (98, 30),
        "cycle_denature": (98, "5–10 s"),
        "anneal_s": "10 s",
        "extend_c": 72,
        "s_per_kb": {"plasmid": 15, "genomic": 30, "cdna": 30},
        "min_s": 15,
        "cycles": {"plasmid": "25–30", "genomic": "30–35", "cdna": "30–35"},
        "final": (72, "5 min"),
        "confidence": "HIGH (Thermo user guide MAN0018859 Rev. B.0, 2022)",
        "source": [THERMO_SUPERFI2],
        "notes": [
            "Thermo: universal annealing at 60 °C, no Tm calculation needed; a 2-step program at 72 °C "
            "is an option for primers over 30 nt without non-binding tails."
        ],
    },
    "OneTaq": {
        "name": "NEB OneTaq",
        "ta": "lower_tm_minus_5",
        "ta_max": 68,
        "two_step_at": 68,
        "denature": (94, 30),
        "cycle_denature": (94, "15–30 s"),
        "anneal_s": "15–60 s",
        "extend_c": 68,
        "s_per_kb": {"plasmid": 60, "genomic": 60, "cdna": 60},
        "min_s": 30,
        "cycles": {"plasmid": "30", "genomic": "30", "cdna": "30"},
        "final": (68, "5 min"),
        "confidence": "MEDIUM (NEB's OneTaq Hot Start M0481 protocol via a reseller copy; M0480 may differ)",
        "source": [NEB_ONETAQ_COPY],
        "notes": [
            "NEB: annealing typically 45–68 °C; optimise with a gradient starting 5 °C below the calculated Tm. "
            "Shown here: the lower Tm − 5 °C as that starting point.",
            "Extension 1 min/kb at 68 °C.",
        ],
    },
    "Taq": {
        "name": "NEB Taq (standard buffer)",
        "ta": "lower_tm_minus_5",
        "ta_max": 68,
        "two_step_at": None,
        "denature": (95, 30),
        "cycle_denature": (95, "15–30 s"),
        "anneal_s": "15–60 s",
        "extend_c": 68,
        "s_per_kb": {"plasmid": 60, "genomic": 60, "cdna": 60},
        "min_s": 30,
        "cycles": {"plasmid": "25–30", "genomic": "25–30", "cdna": "25–30"},
        "final": (68, "5 min"),
        "confidence": "LOW (search summaries of NEB's M0273 page; not read directly). UNVERIFIED",
        "source": [NEB_TAQ],
        "notes": [
            "UNVERIFIED: NEB recommends its Tm Calculator; 'Tm − 5 °C' is the common Taq starting point."
        ],
    },
}
TEMPLATES = ("plasmid", "genomic", "cdna")


def _anneal(rule: dict, tms: list[float], lengths: list[int]) -> tuple[float, str]:
    low = min(tms)
    kind = rule["ta"]
    if kind == "fixed_60":
        return 60.0, "fixed 60 °C (universal annealing)"
    if kind == "lower_tm_plus_3":
        ta, why = low + 3, f"lower Tm {low} + 3 °C"
    elif kind == "phusion":
        ta, why = (
            (low, f"lower Tm {low} (both primers ≤ 20 nt)")
            if max(lengths) <= 20
            else (low + 3, f"lower Tm {low} + 3 °C (a primer > 20 nt)")
        )
    else:
        ta, why = low - 5, f"lower Tm {low} − 5 °C (gradient start)"
    if ta > rule["ta_max"]:
        ta, why = rule["ta_max"], why + f", capped at {rule['ta_max']} °C"
    return round(ta, 1), why


def pcr(
    ref: MapRecord,
    primers: list[dict],
    polymerase: str = "Q5",
    template: str = "plasmid",
) -> dict:
    """primers: two of {"name", "seq"} (5'->3', tails allowed). The pair that faces each other with the
    shortest product is used; every other binding site is listed as a warning."""
    if polymerase not in POLYMERASES:
        raise ValueError(
            f"unknown polymerase {polymerase}; choose one of {', '.join(POLYMERASES)}"
        )
    if template not in TEMPLATES:
        raise ValueError(f"template must be one of {', '.join(TEMPLATES)}")
    if len(primers) != 2:
        raise ValueError("pick exactly two primers")
    if not ref.seq:
        raise ValueError("map has no sequence")
    rule = POLYMERASES[polymerase]
    olig = [
        Oligo(
            p.get("name") or f"primer {i + 1}",
            p.get("name") or f"primer {i + 1}",
            "",
            "",
            "",
            p["seq"],
            "".join(c for c in p["seq"].upper() if c.isalpha()),
            "",
            "",
        )
        for i, p in enumerate(primers)
    ]
    groups, _ = group_oligos(olig)
    hits = match_map(ref, groups)
    by = lambda h, o: o.oligo_id in h.oligo_ids
    L, circ = ref.length, ref.circular
    pairs = []
    for fi, ri in ((0, 1), (1, 0)):
        for f in (h for h in hits if by(h, olig[fi]) and h.strand == "+"):
            for r in (h for h in hits if by(h, olig[ri]) and h.strand == "-"):
                span = r.end - f.start + 1
                if span <= 0:
                    if not circ:
                        continue
                    span += L  # around the origin of a circular map
                pairs.append((span, f, r, fi, ri))
    if not pairs:
        raise ValueError(
            "these two primers do not face each other on this map (one must bind each strand)"
        )
    span, f, r, fi, ri = min(pairs, key=lambda x: x[0])
    fseq, rseq = olig[fi].seq, olig[ri].seq
    # The product carries the primers' own bases (5' mismatches and tails included), and the map only
    # between the two binding parts: the map's bases under a mismatched primer end are not in it.
    s2 = ref.seq + ref.seq
    gap = (r.start - 1) - f.end  # map bases between the two binding parts
    if circ and r.end < f.start:  # the product runs across the origin
        gap += L
    rc = revcomp(rseq)
    middle = s2[f.end : f.end + gap] if gap > 0 else ""
    product_seq = fseq + middle + (rc[-gap:] if gap < 0 else rc)  # gap < 0: the binding parts overlap
    prod = {"product": len(product_seq), "product_seq": product_seq}
    warnings = []
    extra = [h for h in hits if h is not f and h is not r]
    if extra:
        warnings.append(
            f"{len(extra)} more binding site(s) on this map: "
            + ", ".join(f"{h.names[0]} {h.strand} {h.start}–{h.end}" for h in extra[:6])
        )
    if len(pairs) > 1:
        warnings.append(f"{len(pairs)} possible products; the shortest is shown")
    tms = [f.tm, r.tm]
    if None in tms:
        raise ValueError(
            "a primer's binding part has bases other than A, C, G, T: no Tm"
        )
    if abs(tms[0] - tms[1]) > 5:
        warnings.append(
            f"the binding parts' Tms differ by {abs(tms[0] - tms[1]):.1f} °C (more than 5)"
        )
    ta, why = _anneal(rule, tms, [f.binding_length, r.binding_length])
    kb = prod["product"] / 1000
    ext = max(rule["min_s"], math.ceil(rule["s_per_kb"][template] * kb))  # short products: the minimum step
    two_step = bool(rule["two_step_at"] and ta >= rule["two_step_at"])
    steps = [
        {
            "step": "initial denaturation",
            "c": rule["denature"][0],
            "time": f"{rule['denature'][1]} s",
        },
        {
            "step": "denaturation",
            "c": rule["cycle_denature"][0],
            "time": rule["cycle_denature"][1],
            "cycles": rule["cycles"][template],
        },
    ]
    if two_step:
        steps.append(
            {
                "step": "anneal + extend",
                "c": rule["extend_c"],
                "time": f"{ext} s",
                "cycles": rule["cycles"][template],
            }
        )
    else:
        steps += [
            {
                "step": "annealing",
                "c": ta,
                "time": rule["anneal_s"],
                "cycles": rule["cycles"][template],
            },
            {
                "step": "extension",
                "c": rule["extend_c"],
                "time": f"{ext} s",
                "cycles": rule["cycles"][template],
            },
        ]
    steps.append(
        {"step": "final extension", "c": rule["final"][0], "time": rule["final"][1]}
    )
    return {
        "product": prod["product"],
        "product_seq": prod["product_seq"],
        "forward": {
            "name": olig[fi].name,
            "start": f.start,
            "end": f.end,
            "tm": f.tm,
            "binding": f.binding_length,
            "tail": f.tail_length,
            "mismatches": f.mismatches,
        },
        "reverse": {
            "name": olig[ri].name,
            "start": r.start,
            "end": r.end,
            "tm": r.tm,
            "binding": r.binding_length,
            "tail": r.tail_length,
            "mismatches": r.mismatches,
        },
        "map_start": f.start,
        "map_end": r.end,
        "polymerase": rule["name"],
        "template": template,
        "annealing": ta,
        "annealing_rule": why,
        "two_step": two_step,
        "extension_s": ext,
        "program": steps,
        "notes": rule["notes"],
        "confidence": rule["confidence"],
        "sources": rule["source"],
        "tm_model": thermo.MODEL,
        "warnings": warnings,
    }
