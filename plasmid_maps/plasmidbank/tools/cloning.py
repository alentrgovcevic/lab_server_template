"""In-house cloning simulation: restriction digest, ligation, Golden Gate assembly, virtual gel.

Cut geometry (verified against Bio.Restriction.search on both strands, 2026-10-02): for a site
match starting at 0-based i on the top strand, the top strand is cut before index i + fst5 and the
bottom strand before i + size + fst3 (top-strand coordinates). For a match of the reverse
complement, top = i - fst3 and bottom = i + size - fst5. The region between the two cuts is the
overhang: 5' when the top cut comes first, 3' when the bottom cut does, blunt when they coincide.

A fragment's sequence includes the overhang regions at both ends, so ligating A then B is
A.seq + B.seq[len(overhang):] when A's right end and B's left end match (same sequence, same type).
"""

from __future__ import annotations

import io
from datetime import date
import itertools
import re
import warnings
from dataclasses import dataclass, field

import Bio.Restriction as R
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from ..formats import MapRecord
from ..sequtil import IUPAC, revcomp

_IUPAC_RE = {k: "[" + v + "]" if len(v) > 1 else v for k, v in IUPAC.items()}


@dataclass
class End:
    overhang: str  # top-strand sequence of the overhang region ("" for blunt or a natural end)
    kind: str  # "5'", "3'", "blunt" or "natural" (end of a linear molecule)

    def matches(self, other: "End") -> bool:
        return (
            self.kind == other.kind
            and self.kind != "natural"
            and self.overhang == other.overhang
        )

    def flipped(self) -> "End":
        return End(revcomp(self.overhang), self.kind)


@dataclass
class Fragment:
    source: str  # plasmid ID or file name of the parent
    seq: str
    left: End
    right: End
    src_start: int  # 0-based start on the parent (top strand); may be >= L for wrapped fragments
    src_len: int  # parent length
    features: list[dict] = field(
        default_factory=list
    )  # carried features, fragment coordinates
    reversed: bool = False

    @property
    def length(self) -> int:
        return len(self.seq)

    def flip(self) -> "Fragment":
        n = len(self.seq)
        feats = [
            {**f, "start": n - f["end"], "end": n - f["start"], "strand": -f["strand"]}
            for f in self.features
        ]
        return Fragment(
            self.source,
            revcomp(self.seq),
            self.right.flipped(),
            self.left.flipped(),
            self.src_start,
            self.src_len,
            feats,
            not self.reversed,
        )


def _enzyme(name: str):
    enz = getattr(R, name, None)
    if not isinstance(enz, R.Restriction.RestrictionType):
        raise ValueError(f"unknown enzyme {name!r}")
    if enz.fst5 is None or enz.fst3 is None:
        raise ValueError(f"{name} cuts at an unknown position and cannot be simulated")
    if hasattr(enz, "fst5") and str(enz.elucidate()).count("^") > 1:
        raise ValueError(
            f"{name} cuts twice per site (excision enzyme) and is not supported"
        )
    return enz


def _site_regex(site: str) -> re.Pattern:
    return re.compile(
        "(?=(" + "".join(_IUPAC_RE.get(c, c) for c in site.upper()) + "))"
    )


def find_cuts(seq: str, circular: bool, enzyme_names: list[str]) -> list[dict]:
    """All cuts as {enzyme, top, bottom, start (overhang region start), end, kind, overhang, strand}."""
    L = len(seq)
    s2 = (seq + seq[:60]) if circular else seq
    cuts = {}
    for name in enzyme_names:
        enz = _enzyme(name)
        site = str(enz.site)
        size = enz.size
        hits = [(m.start(), "+") for m in _site_regex(site).finditer(s2)]
        rc_site = revcomp(site)
        if rc_site != site:
            hits += [(m.start(), "-") for m in _site_regex(rc_site).finditer(s2)]
        for i, strand in hits:
            if i >= L:
                continue
            if strand == "+":
                top, bottom = i + enz.fst5, i + size + enz.fst3
            else:
                top, bottom = i - enz.fst3, i + size - enz.fst5
            if not circular and (min(top, bottom) < 0 or max(top, bottom) > L):
                continue  # cut would fall off the end of a linear molecule
            a, b = min(top, bottom), max(top, bottom)
            kind = "blunt" if a == b else ("5'" if top < bottom else "3'")
            a_mod = a % L if circular else a
            key = (a_mod, b - a, kind)
            cuts[key] = {
                "enzyme": str(enz),
                "top": top % L + 1 if circular else top + 1,
                "start": a_mod,
                "len": b - a,
                "kind": kind,
                "strand": strand,
                "overhang": (s2 + s2)[a % L : a % L + (b - a)]
                if circular
                else seq[a:b],
            }
    return sorted(cuts.values(), key=lambda c: c["start"])


def _carry_features(ref: MapRecord, start: int, length: int) -> list[dict]:
    """Features of ref lying wholly inside [start, start+length) (start may run past L)."""
    L = ref.length
    out = []
    for f in ref.features:
        fs = f.start
        fe = f.end if f.end >= f.start else f.end + L
        if fe == fs:
            continue
        for k in (0, 1, -1):
            s, e = fs + k * L, fe + k * L
            if start <= s and e <= start + length:
                out.append(
                    {
                        "type": f.type,
                        "label": f.label,
                        "start": s - start,
                        "end": e - start,
                        "strand": f.strand or 0,
                        "source": ref.plasmid_id or "",
                    }
                )
                break
    return out


def digest(
    ref: MapRecord, enzyme_names: list[str]
) -> tuple[list[Fragment], list[dict]]:
    seq, L = ref.seq.upper(), ref.length
    cuts = find_cuts(seq, ref.circular, enzyme_names)
    src = ref.plasmid_id or ref.name
    frags = []
    if not cuts:
        return [], cuts
    if ref.circular:
        for i, c in enumerate(cuts):
            n = cuts[(i + 1) % len(cuts)]
            start = c["start"]
            stop = n["start"] + n["len"]
            if stop <= start or (len(cuts) == 1):
                stop += L
            s2 = seq * 3
            frags.append(
                Fragment(
                    src,
                    s2[start:stop],
                    End(c["overhang"], c["kind"]),
                    End(n["overhang"], n["kind"]),
                    start,
                    L,
                    _carry_features(ref, start, stop - start),
                )
            )
    else:
        bounds = [None] + cuts + [None]
        for a, b in zip(bounds, bounds[1:]):
            start = 0 if a is None else a["start"]
            stop = L if b is None else b["start"] + b["len"]
            left = End("", "natural") if a is None else End(a["overhang"], a["kind"])
            right = End("", "natural") if b is None else End(b["overhang"], b["kind"])
            frags.append(
                Fragment(
                    src,
                    seq[start:stop],
                    left,
                    right,
                    start,
                    L,
                    _carry_features(ref, start, stop - start),
                )
            )
    return frags, cuts


def _has_site(seq: str, enzyme_names: list[str]) -> bool:
    for name in enzyme_names:
        site = str(_enzyme(name).site)
        if _site_regex(site).search(seq) or _site_regex(revcomp(site)).search(seq):
            return True
    return False


@dataclass
class Product:
    seq: str
    circular: bool
    parts: list[dict]  # [{source, orientation, start, length}]
    junctions: list[dict]  # [{overhang, kind, between}]
    features: list[dict]
    warnings: list[str] = field(default_factory=list)


def _join(order: list[Fragment], circular: bool) -> Product:
    seq = order[0].seq
    parts = [
        {
            "source": order[0].source,
            "orientation": "-" if order[0].reversed else "+",
            "start": 0,
            "length": order[0].length,
        }
    ]
    feats = [dict(f) for f in order[0].features]
    junctions = []
    for prev, frag in zip(order, order[1:]):
        k = len(frag.left.overhang)
        offset = len(seq) - k
        junctions.append(
            {
                "overhang": frag.left.overhang,
                "kind": frag.left.kind,
                "between": f"{prev.source} | {frag.source}",
            }
        )
        seq += frag.seq[k:]
        parts.append(
            {
                "source": frag.source,
                "orientation": "-" if frag.reversed else "+",
                "start": offset,
                "length": frag.length,
            }
        )
        feats += [
            {**f, "start": f["start"] + offset, "end": f["end"] + offset}
            for f in frag.features
        ]
    if circular:
        k = len(order[0].left.overhang)
        junctions.append(
            {
                "overhang": order[0].left.overhang,
                "kind": order[0].left.kind,
                "between": f"{order[-1].source} | {order[0].source}",
            }
        )
        seq = seq[: len(seq) - k] if k else seq
        L = len(seq)
        feats = [
            {**f, "start": f["start"] % L, "end": (f["end"] - 1) % L + 1} for f in feats
        ]
    return Product(seq, circular, parts, junctions, feats)


MAX_LIGATION_STEPS = 200_000  # bounds the order/orientation search (security audit 2026-10-08)
MAX_GG_COMBINATIONS = 64  # fragment choices across all inputs


def ligate(fragments: list[Fragment], circular: bool = True) -> list[Product]:
    """Every distinct order/orientation joining all fragments once (first fragment fixed, forward)."""
    if not fragments:
        return []
    first, rest = fragments[0], fragments[1:]
    products = []
    steps = [0]

    def dfs(chain: list[Fragment], remaining: list[Fragment]):
        steps[0] += 1
        if steps[0] > MAX_LIGATION_STEPS:
            raise ValueError("too many possible assemblies to try: use fewer parts, or parts with distinct overhangs")
        if not remaining:
            if not circular or chain[-1].right.matches(chain[0].left):
                products.append(_join(chain, circular))
            return
        for i, f in enumerate(remaining):
            for g in (f, f.flip()):
                if chain[-1].right.matches(g.left):
                    dfs(chain + [g], remaining[:i] + remaining[i + 1 :])

    dfs([first], rest)
    seen, unique = set(), []
    for p in products:
        key = min(p.seq, revcomp(p.seq))
        if key not in seen:
            seen.add(key)
            unique.append(p)
    return unique


def golden_gate(acceptor: MapRecord, parts: list[MapRecord], enzyme: str) -> dict:
    """Type IIS assembly: each input contributes the fragment that no longer carries the site."""
    report = {
        "enzyme": enzyme,
        "inputs": [],
        "products": [],
        "errors": [],
        "warnings": [],
    }
    chosen: list[list[Fragment]] = []
    for ref in [acceptor] + parts:
        frags, cuts = digest(ref, [enzyme])
        name = ref.plasmid_id or ref.name
        clean = [
            f
            for f in frags
            if f.left.kind != "natural"
            and f.right.kind != "natural"
            and not _has_site(f.seq, [enzyme])
        ]
        # a part with no annotation of its own would vanish from the product map
        for f in clean:
            if not any(x["type"] != "source" for x in f.features):
                f.features.append(
                    {"type": "misc_feature", "label": f"cloning module_{name}", "start": 0,
                     "end": f.length, "strand": 1, "source": name}
                )
        report["inputs"].append(
            {
                "source": name,
                "cuts": len(cuts),
                "fragments": [
                    {
                        "length": f.length,
                        "left": f.left.overhang,
                        "right": f.right.overhang,
                        "kept": f in clean,
                    }
                    for f in frags
                ],
            }
        )
        if not cuts:
            report["errors"].append(f"{name}: no {enzyme} site")
        elif not clean:
            report["errors"].append(
                f"{name}: every {enzyme} fragment still contains a {enzyme} site"
            )
        chosen.append(clean)
    if report["errors"]:
        return report
    n_combos = 1
    for c in chosen:
        n_combos *= len(c)
    if n_combos > MAX_GG_COMBINATIONS:
        raise ValueError(f"{n_combos} fragment combinations to try: too many, pick fewer parts")
    for combo in itertools.product(*chosen):
        for p in ligate(list(combo), circular=True):
            ohs = [j["overhang"] for j in p.junctions]
            if len(set(ohs)) != len(ohs):
                p.warnings.append(
                    "the same overhang is used at two junctions: the assembly order is ambiguous"
                )
            for oh in ohs:
                if oh and oh == revcomp(oh):
                    p.warnings.append(f"palindromic overhang {oh} can self-ligate")
            if _has_site(p.seq, [enzyme]):
                p.warnings.append(
                    f"the product still contains a {enzyme} site and would be re-cut"
                )
            report["products"].append(p)
    if not report["products"]:
        ends = [
            f"{f.source}: {f.left.overhang}…{f.right.overhang}"
            for fs in chosen
            for f in fs
        ]
        report["errors"].append(
            "no circular product: the overhangs do not chain (" + "; ".join(ends) + ")"
        )
        report["diagnosis"] = diagnose([fs[0] for fs in chosen if len(fs) == 1])
    return report


FUSION_GAP = 30  # CDS features this close (bp) across a junction are read as one fusion
STOPS = {"TAA", "TAG", "TGA"}


def fusion_frames(p: Product) -> list[dict]:
    """Check CDS parts that are joined end to end in a product (e.g. signal peptide + CDS, or
    CDS + C-terminal tag): the downstream CDS must stay in the upstream reading frame, with no stop
    codon before it. A part that is a base or two short is the classic way a fusion goes wrong.

    Positions in the result are 1-based on the product as given.
    """
    L = len(p.seq)
    out = []
    for strand in (1, -1):
        seq = p.seq if strand == 1 else revcomp(p.seq)
        cds = []
        for f in p.features:
            if f.get("type") != "CDS" or f.get("strand") != strand:
                continue
            s, e = (
                (f["start"], f["end"])
                if strand == 1
                else (L - f["end"], L - f["start"])
            )
            s %= L
            e = e if e > s else e + L  # unwrap a feature over the origin
            cds.append((s, e, f.get("label") or "CDS"))
        cds.sort()
        two = seq + seq
        for (s1, e1, a), (s2, e2, b) in zip(cds, cds[1:]):
            gap = s2 - e1
            if not -10 <= gap <= FUSION_GAP:
                continue
            shift = (s2 - s1) % 3
            pos = (s2 if strand == 1 else L - s2 - 1) % L + 1
            check = {
                "upstream": a,
                "downstream": b,
                "strand": "+" if strand == 1 else "-",
                "junction_at": pos,
                "in_frame": shift == 0,
                "stop_at": None,
            }
            if shift:
                check["message"] = (
                    f"frameshift: {b} is fused after {a} but is {shift} base(s) out of its frame"
                )
            else:
                for i in range(s1, s2, 3):
                    if two[i : i + 3] in STOPS:
                        check["stop_at"] = (i if strand == 1 else L - i - 3) % L + 1
                        last = e1 - 3 <= i < e1  # the upstream CDS's own last codon
                        check["message"] = (
                            f"{a} ends in a stop codon, so {b} is not translated as a fusion"
                            if last
                            else f"stop codon inside the {a}–{b} fusion, before {b}"
                        )
                        break
            if "message" in check:
                out.append(check)
    return out


def flip_product(p: Product) -> Product:
    """The same circular product read on the other strand: parts in reverse order, each flipped."""
    L = len(p.seq)
    n = len(p.parts)
    parts = [
        {
            **x,
            "orientation": "+" if x["orientation"] == "-" else "-",
            "start": (L - x["start"] - x["length"]) % L,
        }
        for x in reversed(p.parts)
    ]
    # junction i joins parts i and i+1; the last one closes the circle
    order = [n - 2 - k for k in range(n - 1)] + [n - 1]
    junctions = [
        {
            **p.junctions[i],
            "overhang": revcomp(p.junctions[i]["overhang"]),
            "between": " | ".join(reversed(p.junctions[i]["between"].split(" | "))),
        }
        for i in order
    ]
    feats = [
        {
            **f,
            "start": (L - f["end"]) % L,
            "end": (L - f["start"] - 1) % L + 1,
            "strand": -f["strand"],
        }
        for f in p.features
    ]
    return Product(
        revcomp(p.seq), p.circular, parts, junctions, feats, list(p.warnings)
    )


def _canon(end: End) -> tuple[str, str]:
    """An overhang and its reverse complement are the same junction, read from either side."""
    return end.kind, min(end.overhang, revcomp(end.overhang))


def diagnose(fragments: list[Fragment]) -> list[dict]:
    """Why a set of fragments does not close into one circle.

    Each fragment joins two overhangs; one circle through every fragment needs every overhang to
    be used by exactly two fragment ends (whatever the part's orientation). Reports overhangs used
    once (nothing to ligate to) or more than twice (ambiguous), and, for each lonely overhang,
    any other lonely overhang one base away: the usual sign of a typo in a primer tail.
    """
    uses: dict[tuple[str, str], list[str]] = {}
    # the overhang as the first part reads it, not the canonical form
    shown: dict[tuple[str, str], str] = {}
    for f in fragments:
        for side, end in (("left", f.left), ("right", f.right)):
            if end.kind in ("natural", "blunt"):
                continue
            shown.setdefault(_canon(end), end.overhang)
            uses.setdefault(_canon(end), []).append(
                f"{f.source} ({side} end, {end.overhang})"
            )
    out = []
    lonely = [k for k, v in uses.items() if len(v) == 1]
    for key, who in uses.items():
        kind, oh = key
        if len(who) == 1:
            near = [
                shown[(k2, o)]
                for k2, o in lonely
                if (k2, o) != key
                and len(o) == len(oh)
                and min(
                    sum(a != b for a, b in zip(oh, o)),
                    sum(a != b for a, b in zip(oh, revcomp(o))),
                )
                == 1
            ]
            out.append(
                {
                    "overhang": shown[key],
                    "problem": "no partner",
                    "ends": who,
                    "one_base_from": near,
                }
            )
        elif len(who) > 2:
            out.append(
                {
                    "overhang": shown[key],
                    "problem": "used more than twice",
                    "ends": who,
                    "one_base_from": [],
                }
            )
    return out


_MONTHS = "JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC".split()


def product_genbank(p: Product, name: str = "assembly", description: str = "") -> str:
    rec = SeqRecord(
        Seq(p.seq),
        id=name[:16],
        name=name[:16],
        description=description or name,
        annotations={
            "molecule_type": "DNA",
            "topology": "circular" if p.circular else "linear",
            # without it Biopython writes 01-JAN-1980 (seen on the share, 2026-10-22); months spelled out
            # so a non-English locale cannot change them
            "date": f"{date.today():%d}-{_MONTHS[date.today().month - 1]}-{date.today():%Y}",
        },
    )
    L = len(p.seq)
    for f in p.features:
        s, e = f["start"], f["end"]
        strand = f["strand"] or None
        if e <= s and p.circular:  # wraps the origin of the product
            loc = CompoundLocation(
                [SimpleLocation(s, L, strand), SimpleLocation(0, e, strand)]
            )
        else:
            loc = SimpleLocation(s, e, strand)
        rec.features.append(
            SeqFeature(
                loc,
                type=f["type"],
                qualifiers={
                    "label": [f["label"]],
                    "note": [f"from {f.get('source', '')}"],
                },
            )
        )
    # How it was made, so the record survives a GenBank export (GenBank has no history field)
    rec.annotations["comment"] = (
        f"simulated with plasmidbank on {date.today().isoformat()}; parts in order: "
        + ", ".join(
            f"{x['source']} ({'reverse complement' if x['orientation'] == '-' else 'as drawn'})"
            for x in p.parts
        )
        + "; junctions: "
        + ", ".join(f"{j['overhang'] or 'blunt'} ({j['between']})" for j in p.junctions)
    )
    out = io.StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SeqIO.write(rec, out, "genbank")
    return out.getvalue()


# ── virtual gel ────────────────────────────────────────────────────────────────

# UNVERIFIED: ladder band sizes from memory of the vendors' product sheets; check before relying on them.
LADDERS = {
    "1 kb (NEB N3232)": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500],
    "GeneRuler 1 kb (Thermo SM0311)": [
        10000,
        8000,
        6000,
        5000,
        4000,
        3500,
        3000,
        2500,
        2000,
        1500,
        1000,
        750,
        500,
        250,
    ],
    "100 bp (NEB N3231)": [
        1517,
        1200,
        1000,
        900,
        800,
        700,
        600,
        500,
        400,
        300,
        200,
        100,
    ],
}


def gel_lane(ref: MapRecord, enzyme_names: list[str]) -> dict:
    if not enzyme_names:
        return {
            "label": f"{ref.plasmid_id} uncut",
            "bands": [],
            "uncut": True,
            "note": "uncut plasmid runs as supercoiled/nicked forms; its position cannot be predicted from size",
        }
    frags, cuts = digest(ref, enzyme_names)
    sizes = (
        sorted(
            (
                f.length - len(f.right.overhang) if ref.circular else f.length
                for f in frags
            ),
            reverse=True,
        )
        if frags
        else [ref.length]
    )
    return {
        "label": f"{ref.plasmid_id or ref.name} {'+'.join(enzyme_names)}",
        "bands": sizes,
        "cuts": len(cuts),
        "uncut": False,
        "note": "" if cuts else "no cut: the plasmid stays uncut",
    }
