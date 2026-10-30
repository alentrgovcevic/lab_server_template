"""Restriction sites via Biopython's Bio.Restriction (REBASE), circular-aware.

Biopython reports each cut as the 1-based position of the first base after the cut on the
top strand. Curated sets keep the list readable: the full commercial set is ~600 enzymes,
most of them isoschizomers.
"""

from __future__ import annotations

import re

import Bio.Restriction as R
from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Seq import Seq

from ..formats import MapRecord

COMMON = (
    "AatII AflII AgeI ApaI ApaLI AscI AseI AvrII BamHI BclI BglII BsiWI BspEI BsrGI BstBI "
    "BstEII ClaI DraI EagI EcoRI EcoRV FseI HindIII HpaI KpnI MfeI MluI NcoI NdeI NheI NotI "
    "NruI NsiI PacI PmeI PstI PvuI PvuII SacI SacII SalI SbfI ScaI SfiI SmaI SnaBI SpeI SphI "
    "StuI SwaI XbaI XhoI XmaI"
).split()
# Type IIS enzymes used in Golden Gate / MoClo work; one name per isoschizomer group.
GOLDEN_GATE = "BsaI BpiI Esp3I SapI AarI BtgZI".split()  # BpiI = BbsI, Esp3I = BsmBI

SETS = {"common": COMMON, "unique": COMMON, "goldengate": GOLDEN_GATE}


def enzymes_for(set_name: str, custom: str = "") -> list[str]:
    if custom.strip():
        names = [n.strip() for n in custom.replace(",", " ").split()]
        bad = [
            n
            for n in names
            if not isinstance(getattr(R, n, None), R.Restriction.RestrictionType)
        ]
        if bad:
            raise ValueError(f"unknown enzyme(s): {', '.join(bad)}")
        return names
    if set_name not in SETS:
        raise ValueError(f"unknown set {set_name!r}; use one of {', '.join(SETS)}")
    return SETS[set_name]


def overhang(enz) -> str:
    """'5′ AATT', '3′ TGCA', 'blunt' (a clicked site shows its overhang)."""
    if enz.is_blunt():
        return "blunt"
    if enz.is_5overhang():
        return f"5′ {enz.ovhgseq}"
    if enz.is_3overhang():
        return f"3′ {enz.ovhgseq}"
    return "not known"


def find_sites(
    ref: MapRecord, set_name: str = "unique", custom: str = ""
) -> list[dict]:
    names = enzymes_for(set_name, custom)
    batch = RestrictionBatch([getattr(R, n) for n in names])
    hits = Analysis(batch, Seq(ref.seq), linear=not ref.circular).full()
    out = []
    for enz, cuts in hits.items():
        if not cuts:
            continue
        if set_name == "unique" and not custom.strip() and len(cuts) != 1:
            continue
        out.append(
            {
                "enzyme": str(enz),
                "site": enz.elucidate(),
                "recognition": str(enz.site),
                "overhang": overhang(enz),
                "cuts": sorted(cuts),
                "n": len(cuts),
                # for the sequence view: where the recognition sequence sits (0-based starts,
                # either strand) and the bottom-strand cut, which is the top cut minus ovhg (5' overhangs are
                # negative in Biopython: EcoRI G^AATT_C has ovhg -4); orientation does not change that
                "spans": _site_spans(ref.seq, ref.circular, str(enz.site)),
                "size": len(enz.site),
                "ovhg": int(enz.ovhg or 0),
                "methylation": methylation_context(
                    ref.seq, ref.circular, str(enz.site)
                ),
            }
        )
    return sorted(out, key=lambda d: (d["cuts"][0], d["enzyme"]))


# ── methylation context ─────────────────────────────────────────
# Dam methylates the A of GATC; Dcm the second C of CCWGG (Addgene, "Plasmids 101: Methylation and
# Restriction Enzymes", 2016-06-30). A recognition site can be blocked by a motif that the
# FLANKING bases create: "Cleavage by XbaI may be blocked due to methylation if the enzyme's
# recognition site (TCTAGA) is preceded by GA or followed by TC" (same source). Only that
# context-made overlap is flagged here. Whether an enzyme whose own site always contains GATC or
# CCWGG is blocked is a property of the enzyme: that needs the supplier's sensitivity table, which
# we do not have, so those are deliberately not flagged.
_METH = {"Dam": re.compile("(?=(GATC))"), "Dcm": re.compile("(?=(CC[AT]GG))")}
_IUPAC_RE = {
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[CG]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}


def _site_spans(seq: str, circular: bool, site: str) -> list[int]:
    """0-based starts of a recognition site on either strand (wrapping the origin if circular)."""
    n, L = len(site), len(seq)
    s2 = seq + seq[: n - 1] if circular else seq
    rc = site.translate(str.maketrans("ACGTRYKMBDHVN", "TGCAYRMKVHDBN"))[::-1]
    starts = set()
    for pat in {site, rc}:
        rx = re.compile("(?=(" + "".join(_IUPAC_RE.get(c, c) for c in pat) + "))")
        starts.update(m.start() % L for m in rx.finditer(s2))
    return sorted(starts)


def methylation_context(seq: str, circular: bool, site: str) -> list[dict]:
    """Sites where a Dam/Dcm motif overlaps the recognition site but runs past its edge."""
    seq = seq.upper()
    n, L = len(site), len(seq)
    pad = 5
    out = []
    for a in _site_spans(seq, circular, site.upper()):
        if circular:
            window = (seq[-pad:] + seq + seq[:pad])[a : a + n + 2 * pad]
        else:
            lo = max(0, a - pad)
            window = " " * (pad - (a - lo)) + seq[lo : a + n + pad]
        for name, rx in _METH.items():
            for m in rx.finditer(window):
                b, e = m.start(), m.start() + len(m.group(1))
                inside = pad <= b and e <= pad + n
                overlaps = b < pad + n and e > pad
                if overlaps and not inside:
                    out.append(
                        {"position": a + 1, "motif": name, "context": window.strip()}
                    )
    return out


# ── enzyme finder: ApE-style conditions across several maps ──────
# A user on Reddit praised ApE for "sites that appear only once in your vector and not at all in
# your insert", limited "to only the enzymes in your lab". Conditions are ANDed; NOT is max=0.


def commercial() -> list[str]:
    return sorted(str(e) for e in R.CommOnly)


def select_enzymes(
    conditions: list[dict], names: list[str], one_per_site: bool = True, prefer=()
) -> list[dict]:
    """Enzymes meeting every condition.

    Each condition is {"map": MapRecord, "min": int, "max": int | None, "region": (start, end) | None}:
    the number of cuts in the whole map must lie in [min, max], and when a region is given
    (1-based, inclusive, may wrap the origin) at least one of those cuts must fall inside it.
    With one_per_site, isoschizomers collapse to one name, preferring names in ``prefer``.
    """
    batch = RestrictionBatch([getattr(R, n) for n in names])
    found = []
    per_map = [
        Analysis(batch, Seq(c["map"].seq), linear=not c["map"].circular).full()
        for c in conditions
    ]
    for enz in batch:
        counts = []
        ok = True
        for c, hits in zip(conditions, per_map):
            cuts = sorted(hits.get(enz, []))
            hi = c.get("max")
            if len(cuts) < c.get("min", 0) or (hi is not None and len(cuts) > hi):
                ok = False
                break
            reg = c.get("region")
            if reg:
                s, e = reg
                inside = [
                    x for x in cuts if (s <= x <= e if s <= e else x >= s or x <= e)
                ]
                if not inside:
                    ok = False
                    break
            counts.append(cuts)
        if ok:
            found.append(
                {
                    "enzyme": str(enz),
                    "site": str(enz.site),
                    "pattern": enz.elucidate(),
                    "cuts": counts,
                    "suppliers": list(enz.suppl),
                }
            )
    if one_per_site:
        best: dict[str, dict] = {}
        pref = set(prefer)
        for f in sorted(found, key=lambda f: (f["enzyme"] not in pref, f["enzyme"])):
            key = f["pattern"]
            if key in best:
                best[key].setdefault("isoschizomers", []).append(f["enzyme"])
            else:
                best[key] = f
        found = list(best.values())
    return sorted(found, key=lambda f: f["enzyme"])
