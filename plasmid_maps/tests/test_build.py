"""Map builder: Gibson design, documenting a cloning, a map from a pasted sequence.

The Gibson tests do not trust the design's own product: they rebuild it the way the reaction does.
Each PCR piece becomes forward tail + template + reverse tail, a synthetic fragment is what the
design says to order, and neighbouring fragments are fused on the longest end overlap (>= 15 nt).
"""

import random

import pytest

from plasmidbank.formats import Feature, MapRecord
from plasmidbank.sequtil import canonical_hash, revcomp
from plasmidbank.tools import build as B


def _map(n, seed, pid, circular=True, features=()):
    rnd = random.Random(seed)
    s = "".join(rnd.choice("ACGT") for _ in range(n))
    return MapRecord(
        "x.gb",
        pid,
        "PMC",
        "filename",
        [],
        "genbank",
        "ok",
        seq=s,
        topology="circular" if circular else "linear",
        features=list(features),
    )


def _fuse(frags):
    """Join fragments in order on their end overlaps, then close the circle."""

    def ov(a, b):
        for k in range(min(len(a), len(b), 80), 14, -1):
            if a.endswith(b[:k]):
                return k
        raise AssertionError("fragments do not overlap")

    s = frags[0]
    for f in frags[1:]:
        s += f[ov(s, f) :]
    k = ov(s, frags[0])
    return s[:-k]


def _amplicons(r, pieces):
    """What the reaction actually contains, in assembly order."""
    prim = {(p["piece"], p["dir"]): p for p in r["primers"]}
    order = {o["name"]: o["sequence"] for o in r["order_fragments"]}
    out = []
    for p in pieces:
        if p.synthetic:
            out.append(order[p.name])
        elif p.pcr:
            f, rv = prim[(p.name, "F")], prim[(p.name, "R")]
            assert p.seq.startswith(f["sequence"]) and p.seq.endswith(
                revcomp(rv["sequence"])
            )
            out.append(f["tail"] + p.seq + revcomp(rv["tail"]))
        else:
            out.append(p.seq)
    return out


def _vector():
    v = _map(4000, 1, "PMC-00100")
    s = v.seq
    # two unique sites: EcoRI at 1000, XhoI at 1500 (make sure neither occurs elsewhere)
    for site in ("GAATTC", "CTCGAG"):
        s = s.replace(site, "GAAGTC" if site == "GAATTC" else "CTCAAG")
    s = s[:1000] + "GAATTC" + s[1006:1500] + "CTCGAG" + s[1506:]
    v.seq = s
    v.features = [
        Feature("CDS", "AmpR", 2000, 2800, 1),
        Feature("misc_feature", "stuffer", 1100, 1400, 1),
    ]
    return v


@pytest.mark.parametrize("opened", ["enzymes", "pcr"])
@pytest.mark.parametrize("n_inserts", [1, 3])
def test_gibson_fragments_rebuild_the_product(opened, n_inserts):
    v = _vector()
    t = _map(3000, 7, "PMC-00200", features=[Feature("CDS", "GFP", 500, 1220, 1)])
    ins = [
        B.region_insert(t, 501, 1220, "+", "GFP"),
        B.region_insert(t, 2001, 2400, "-", "rev part"),
        B.pasted_insert(_map(300, 9, "x").seq, "gBlock"),
    ][:n_inserts]
    vo = (
        {"enzymes": ["EcoRI", "XhoI"]}
        if opened == "enzymes"
        else {"start": 1007, "end": 1500}
    )
    r = B.build(v, vo, ins, "gibson", overlap=25)
    backbone, _ = B.open_vector(
        v, "gibson", vo.get("enzymes"), vo.get("start", 0), vo.get("end", 0)
    )
    frags = _amplicons(r, [backbone, *ins])
    assert canonical_hash(_fuse(frags), True) == canonical_hash(r["seq"], True)
    assert r["length"] == len(r["seq"]) and len(r["junctions"]) == n_inserts + 1
    # the stuffer between the sites is gone, AmpR and GFP came along
    assert (
        "LOCUS" in r["genbank"]
        and "AmpR" in r["genbank"]
        and "stuffer" not in r["genbank"]
    )
    if n_inserts >= 1:
        assert "GFP" in r["genbank"]
    for j in r["junctions"]:
        assert len(j["overhang"]) == 25


def test_enzyme_backbone_keeps_3prime_ends_only():
    v = _vector()
    bb, info = B.open_vector(v, "gibson", ["EcoRI", "XhoI"])
    # EcoRI and XhoI leave 5' overhangs, which the exonuclease removes
    assert bb.seq.startswith(v.seq[1505:1520]) and bb.seq.endswith(v.seq[985:1001])
    doc, _ = B.open_vector(v, "document", ["EcoRI", "XhoI"])
    assert doc.seq.startswith(v.seq[1501:1520]) and doc.seq.endswith(v.seq[985:1001])
    assert info["dropped_bp"] == 4000 - len(bb.seq)
    no_site = _map(4000, 3, "x")
    no_site.seq = no_site.seq.replace("GAATTC", "GAAGTC")
    with pytest.raises(ValueError, match="cuts .* 0 times"):
        B.open_vector(no_site, "gibson", ["EcoRI"])
    two = _vector()
    two.seq = two.seq[:3000] + "GAATTC" + two.seq[3006:]
    with pytest.raises(ValueError, match="cuts .* 2 times"):
        B.open_vector(two, "gibson", ["EcoRI"])


def test_document_mode_joins_without_primers():
    v = _vector()
    t = _map(3000, 7, "PMC-00200")
    r = B.build(
        v, {"enzymes": ["EcoRI", "XhoI"]}, [B.region_insert(t, 101, 400)], "document"
    )
    assert r["primers"] == [] and r["order_fragments"] == []
    assert r["seq"] == v.seq[1501:] + v.seq[:1001] + t.seq[100:400]


def test_pcr_vector_insertion_without_removing_anything():
    v = _vector()
    r = B.build(
        v,
        {"start": 2000, "end": 1999},
        [B.pasted_insert("CATCACCATCACCATCAC", "His6")],
        "document",
    )
    assert r["length"] == 4018 and "CATCACCATCACCATCAC" in r["seq"] + r["seq"]


def test_gibson_checks():
    v, t = _vector(), _map(3000, 7, "PMC-00200")
    ins = B.region_insert(t, 501, 1220)
    with pytest.raises(ValueError, match="overlap"):
        B.build(v, {"enzymes": ["EcoRI"]}, [ins], "gibson", overlap=10)
    with pytest.raises(ValueError, match="1 to 5"):
        B.build(v, {"enzymes": ["EcoRI"]}, [ins] * 6, "gibson")
    r = B.build(v, {"enzymes": ["EcoRI"]}, [ins], "gibson")
    assert not any(
        "share" in w for w in r["warnings"]
    )  # random overlaps: no false alarm
    # two overlaps sharing a 12-mer (directly or reverse-complemented) can pair with the wrong partner
    a = "ACGTTGCAAGGCTTAACGGATCCAT"
    assert B._cross_talk([a, "TTTTT" + a[5:17] + "GGGGGGGG"]) != []
    assert B._cross_talk([a, "TTTTT" + revcomp(a[5:17]) + "GGGGGGGG"]) != []
    assert B._cross_talk([a, "CCCCCCCCCCCCCCCCCCCCCCCCC"]) == []


def test_map_from_pasted_text():
    raw = B.map_from_text(">pTest my plasmid\nACGTACGTAC GTACGTACGT\nACGTACGTAC\n", "")
    assert raw["name"] == "pTest" and raw["length"] == 30 and "LOCUS" in raw["genbank"]
    from datetime import date
    locus = raw["genbank"].splitlines()[0]  # dated today, not Biopython's 01-JAN-1980
    assert locus.endswith(f"{date.today():%d}-{date.today().strftime('%b').upper()}-{date.today():%Y}"), locus
    gb = B.map_from_text(raw["genbank"].replace("pasted sequence", "x"), "renamed")
    assert gb["length"] == 30 and gb["name"] == "renamed"
    with pytest.raises(ValueError, match="not DNA"):
        B.map_from_text("ACGTXYZ" * 5, "bad")
    with pytest.raises(ValueError, match="too short"):
        B.map_from_text("ACGT", "short")
