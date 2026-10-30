import random
from pathlib import Path

import pytest

from plasmidbank.formats import Feature, MapRecord, read_map
from plasmidbank.sequtil import revcomp

SAMPLE = Path(__file__).resolve().parents[2] / "examples" / "bank"


def _rand(n, seed, avoid=("GGTCTC", "GAGACC", "GAATTC", "GAAGAC", "GTCTTC")):
    rnd = random.Random(seed)
    s = "".join(rnd.choice("ACGT") for _ in range(n))
    for a in avoid:
        s = s.replace(a, a[:-1] + ("A" if a[-1] != "A" else "T"))
    for a in avoid:  # replacements can recreate a site; repeat until clean
        while a in s:
            s = s.replace(a, a[:-1] + ("A" if a[-1] != "A" else "T"))
    return s


def _map(seq, pid="PMC-00001", circular=True, features=()):
    return MapRecord(
        "x.gb",
        pid,
        "PMC",
        "filename",
        [],
        "genbank",
        "ok",
        seq=seq,
        topology="circular" if circular else "linear",
        features=list(features),
    )


# ── cloning: digest, ligation, Golden Gate ───────────────────────────────────


def test_cut_geometry_matches_biopython_both_strands():
    from Bio.Restriction import BsaI, EcoRI, KpnI
    from Bio.Seq import Seq
    from plasmidbank.tools.cloning import find_cuts

    s = _rand(400, 1)
    s = (
        s[:50]
        + "GGTCTCA"
        + "GGAG"
        + s[61:200]
        + "CGCT"
        + "T"
        + "GAGACC"
        + s[211:300]
        + "GAATTC"
        + s[306:350]
        + "GGTACC"
        + s[356:]
    )
    cuts = find_cuts(s, False, ["BsaI", "EcoRI", "KpnI"])
    for enz in (BsaI, EcoRI, KpnI):
        mine = sorted(c["top"] for c in cuts if c["enzyme"] == str(enz))
        assert mine == sorted(enz.search(Seq(s), linear=True))
    oh = {(c["enzyme"], c["overhang"], c["kind"]) for c in cuts}
    assert ("BsaI", "GGAG", "5'") in oh and ("BsaI", "CGCT", "5'") in oh
    assert ("EcoRI", "AATT", "5'") in oh and ("KpnI", "GTAC", "3'") in oh


def test_digest_fragments_sum_to_plasmid_and_gel_lane():
    from plasmidbank.tools.cloning import digest, gel_lane

    s = _rand(3000, 2)
    s = s[:500] + "GAATTC" + s[506:2000] + "GAATTC" + s[2006:]
    ref = _map(s)
    frags, cuts = digest(ref, ["EcoRI"])
    assert len(cuts) == 2 and len(frags) == 2
    lane = gel_lane(ref, ["EcoRI"])
    assert sorted(lane["bands"]) == [1500, 1500] and sum(lane["bands"]) == 3000
    assert gel_lane(_map(_rand(2000, 3)), ["EcoRI"])["bands"] == [2000]


def test_restriction_ligation_recircularises():
    from plasmidbank.tools.cloning import digest, ligate

    s = _rand(3000, 4)
    s = s[:500] + "GAATTC" + s[506:2000] + "GAATTC" + s[2006:]
    frags, _ = digest(_map(s), ["EcoRI"])
    products = ligate(frags)
    seqs = {min(p.seq, revcomp(p.seq)) for p in products}
    # original plasmid and the variant with the small fragment flipped
    assert len(products) == 2 and any(len(p.seq) == 3000 for p in products)
    from plasmidbank.sequtil import canonical_hash

    assert canonical_hash(s, True) in {canonical_hash(p.seq, True) for p in products}


def test_golden_gate_assembles_part_into_acceptor_with_features():
    from plasmidbank.tools.cloning import golden_gate, product_genbank

    body = _rand(600, 5)
    bb = _rand(2500, 6)
    part_seq = (
        bb[:800] + "GGTCTCA" + "GGAG" + body + "CGCT" + "T" + "GAGACC" + bb[800:1500]
    )
    part = _map(
        part_seq, "EC00001", features=[Feature("CDS", "insert", 811, 811 + 600, 1)]
    )
    dropout = _rand(400, 7)
    acc_seq = (
        bb[1500:2400] + "GGAG" + "TGAGACC" + dropout + "GGTCTCA" + "CGCT" + bb[2400:]
    )
    acc = _map(
        acc_seq,
        "PMC-00002",
        features=[Feature("CDS", "lacZ", 911 + 7, 911 + 7 + 400, 1)],
    )
    r = golden_gate(acc, [part], "BsaI")
    assert r["errors"] == [] and len(r["products"]) == 1
    p = r["products"][0]
    expect_len = (
        len(acc_seq) - (4 + 7 + 400 + 7) + 4 + 600
    )  # backbone keeps one copy of each overhang
    assert len(p.seq) == expect_len
    assert body in p.seq + p.seq or revcomp(body) in p.seq + p.seq
    assert {j["overhang"] for j in p.junctions} == {"GGAG", "CGCT"}
    labels = {f["label"] for f in p.features}
    assert "insert" in labels and "lacZ" not in labels  # dropout cassette is gone
    gb = product_genbank(p, "GG_test")
    assert "LOCUS" in gb and "insert" in gb and "circular" in gb.splitlines()[0]


@pytest.mark.parametrize("annotated", [True, False])
def test_golden_gate_unannotated_part_gets_a_cloning_module(annotated):
    from plasmidbank.tools.cloning import golden_gate

    body, bb = _rand(600, 5), _rand(2500, 6)
    part_seq = bb[:800] + "GGTCTCA" + "GGAG" + body + "CGCT" + "T" + "GAGACC" + bb[800:1500]
    part = _map(part_seq, "EC00001", features=[Feature("CDS", "insert", 811, 1411, 1)] if annotated else [])
    acc_seq = bb[1500:2400] + "GGAG" + "TGAGACC" + _rand(400, 7) + "GGTCTCA" + "CGCT" + bb[2400:]
    acc = _map(acc_seq, "PMC-00002", features=[Feature("CDS", "lacZ", 918, 1318, 1)])
    p = golden_gate(acc, [part], "BsaI")["products"][0]
    mods = {f["label"]: f for f in p.features if f["label"].startswith("cloning module_")}
    assert "cloning module_PMC-00002" in mods  # the backbone fragment carries nothing of its own
    if annotated:
        assert "cloning module_EC00001" not in mods
    else:
        f = mods["cloning module_EC00001"]
        span = (p.seq + p.seq)[f["start"]: f["end"] if f["end"] > f["start"] else f["end"] + len(p.seq)]
        assert body in span or revcomp(body) in span


def test_golden_gate_reports_unmatched_overhangs():
    from plasmidbank.tools.cloning import golden_gate

    body = _rand(300, 8)
    bb = _rand(2000, 9)
    part = _map(
        bb[:600] + "GGTCTCA" + "GGAG" + body + "TTTT" + "T" + "GAGACC" + bb[600:1200],
        "EC00002",
    )
    acc = _map(
        bb[1200:1800]
        + "GGAG"
        + "TGAGACC"
        + _rand(200, 10)
        + "GGTCTCA"
        + "CGCT"
        + bb[1800:],
        "PMC-00003",
    )
    r = golden_gate(acc, [part], "BsaI")
    assert r["products"] == [] and any("do not chain" in e for e in r["errors"])


# ── ORFs, translation, MSA ───────────────────────────────────────────────────


def test_orf_finder_both_strands_and_origin():
    from plasmidbank.tools.analysis import find_orfs

    cds = "ATG" + "GCT" * 150 + "TAA"
    filler = "T" * 60 + "C" * 60  # no ATG in any frame around it
    s = cds[200:] + filler * 10 + cds[:200]  # ORF spans the origin
    s2 = filler * 5 + revcomp(cds) + filler * 5
    orfs = find_orfs(_map(s), min_aa=100)
    assert any(
        o["strand"] == "+" and o["length_aa"] == 151 and o["start"] == len(s) - 200 + 1
        for o in orfs
    )
    orfs2 = find_orfs(_map(s2, circular=False), min_aa=100)
    assert any(
        o["strand"] == "-" and o["start"] == 601 and o["end"] == 600 + len(cds)
        for o in orfs2
    )


def test_translate_region_flags_frame():
    from plasmidbank.tools.analysis import translate_region

    ref = _map("ATGGCTTAA" + "C" * 100)
    assert translate_region(ref, 1, 9)["protein"] == "MA*"
    assert translate_region(ref, 1, 10)["warning"]


@pytest.mark.skipif(
    not __import__("shutil").which("mafft"), reason="MAFFT not installed"
)
def test_mafft_alignment_rotates_circular_sequences():
    from plasmidbank.tools.analysis import multiple_alignment

    s = _rand(1500, 11)
    variant = list(s)
    variant[700] = "A" if s[700] != "A" else "C"
    variant = "".join(variant)
    rotated = revcomp(variant[400:] + variant[:400])
    r = multiple_alignment([("ref", s, True), ("clone", rotated, True)])
    assert r["columns"] == 1500 and r["identity"] > 0.99
    assert sum(1 for c in r["consensus"] if c.islower()) == 1


# ── primer walking ───────────────────────────────────────────────────────────


def test_walk_tiles_cover_region_and_reuse_collection():
    from plasmidbank.oligos import Oligo
    from plasmidbank.tools.design import design_walk
    from plasmidbank.tools.primers import group_oligos, match_map

    s = _rand(6000, 12)
    ref = _map(s)
    coll = Oligo(
        "PM-00009", "walk_known", "", "", "", s[1380:1402], s[1380:1402], "", ""
    )
    groups, _ = group_oligos([coll])
    r = design_walk(ref, 1001, 3000, match_map(ref, groups), step=550)
    fwd = [t for t in r["tiles"] if t["strand"] == "+"]
    assert len(fwd) >= 4 and r["tiles"][-1]["strand"] == "-"
    assert any(t["source"] == "collection" and t["ids"] == ["PM-00009"] for t in fwd)
    threes = sorted(t["three_prime"] for t in fwd if t["source"] != "none")
    assert all(b - a <= 700 for a, b in zip(threes, threes[1:]))


# ── edit and save ────────────────────────────────────────────────────────────


def test_apply_edits_never_changes_sequence(tmp_path):
    import shutil
    from Bio import SeqIO
    import io
    from plasmidbank.tools.edit import apply_edits

    src = next((SAMPLE / "pmc_gb").glob("PMC-07001*"))
    p = tmp_path / "x.gb"
    shutil.copy(src, p)
    m = read_map(p)
    text, summary = apply_edits(
        str(p),
        {
            "features": [
                {"index": 0, "action": "delete"},
                {"index": 1, "action": "update", "label": "AmpR (edited)"},
            ],
            "add": [
                {
                    "type": "misc_feature",
                    "label": "my site",
                    "start": 10800,
                    "end": 20,
                    "strand": 1,
                }
            ],
        },
    )
    rec = SeqIO.read(io.StringIO(text), "genbank")
    assert str(rec.seq).upper() == m.seq
    labels = [
        f.qualifiers.get("label", [""])[0] for f in rec.features if f.type != "source"
    ]
    assert (
        "AmpR (edited)" in labels
        and "my site" in labels
        and len(labels) == len(m.features)
    )
    wrapped = next(f for f in rec.features if f.qualifiers.get("label") == ["my site"])
    assert len(wrapped.location.parts) == 2
    assert p.read_bytes() == src.read_bytes()  # original untouched


def test_golden_gate_diagnosis_names_lonely_overhangs_and_one_base_typos():
    from plasmidbank.tools.cloning import golden_gate

    body = _rand(300, 11)
    bb = _rand(2000, 12)
    # the part's 3' fusion site has a typo: CGCA instead of the acceptor's CGCT
    part = _map(
        bb[:600] + "GGTCTCA" + "GGAG" + body + "CGCA" + "T" + "GAGACC" + bb[600:1200],
        "EC00003",
    )
    acc = _map(
        bb[1200:1800] + "GGAG" + "TGAGACC" + _rand(200, 13) + "GGTCTCA" + "CGCT" + bb[1800:],
        "PMC-00004",
    )
    r = golden_gate(acc, [part], "BsaI")
    assert r["products"] == []
    lonely = {d["overhang"]: d for d in r["diagnosis"] if d["problem"] == "no partner"}
    # GGAG pairs up; the two 3' sites are lonely and each is one base from the other
    assert "GGAG" not in lonely and "CTCC" not in lonely
    assert len(lonely) == 2
    for d in lonely.values():
        assert len(d["one_base_from"]) == 1


def test_flip_product_is_the_same_circle_read_on_the_other_strand():
    from plasmidbank.sequtil import canonical_hash, revcomp
    from plasmidbank.tools.cloning import flip_product, golden_gate

    body, bb = _rand(600, 5), _rand(2500, 6)
    part = _map(
        revcomp(bb[:800] + "GGTCTCA" + "GGAG" + body + "CGCT" + "T" + "GAGACC" + bb[800:1500]),
        "EC00001",
        features=[Feature("CDS", "insert", 700, 1300, 1)],
    )
    acc = _map(
        bb[1500:2400] + "GGAG" + "TGAGACC" + _rand(400, 7) + "GGTCTCA" + "CGCT" + bb[2400:],
        "PMC-00002",
    )
    p = golden_gate(acc, [part], "BsaI")["products"][0]
    q = flip_product(p)
    assert canonical_hash(q.seq, True) == canonical_hash(p.seq, True)
    assert q.seq == revcomp(p.seq)
    assert {j["overhang"] for j in q.junctions} == {revcomp(j["overhang"]) for j in p.junctions}
    for f, g in zip(p.features, q.features):  # the insert keeps its sequence, read on the other strand
        L = len(p.seq)
        assert g["strand"] == -f["strand"]
        assert revcomp((p.seq * 2)[f["start"] : f["start"] + (f["end"] - f["start"]) % L or L]) == (
            q.seq * 2
        )[g["start"] : g["start"] + (g["end"] - g["start"]) % L or L]
    assert flip_product(q).seq == p.seq


def _orf(n_codons, seed):
    """n sense codons with no stop codon."""
    import random

    rnd = random.Random(seed)
    codons = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"]
    sense = [c for c in codons if c not in ("TAA", "TAG", "TGA")]
    return "".join(rnd.choice(sense) for _ in range(n_codons))


def test_fusion_frame_check_flags_a_part_one_base_short():
    from plasmidbank.tools.cloning import Product, fusion_frames

    sp, cds = "ATG" + _orf(20, 1), _orf(60, 2) + "TAA"
    scar = "AGGT" + "CA"  # 6 bp between the two CDS features keeps the frame

    def product(link):
        seq = _rand(300, 3) + sp + link + cds + _rand(300, 4)
        a0 = 300
        b0 = a0 + len(sp) + len(link)
        feats = [
            {"type": "CDS", "label": "SP", "start": a0, "end": a0 + len(sp), "strand": 1},
            {"type": "CDS", "label": "GFP", "start": b0, "end": b0 + len(cds), "strand": 1},
        ]
        return Product(seq, True, [], [], feats)

    assert fusion_frames(product(scar)) == []
    bad = fusion_frames(product(scar[:-1]))  # the part lost one base
    assert len(bad) == 1 and not bad[0]["in_frame"] and "out of its frame" in bad[0]["message"]
    stop = fusion_frames(product("TAA" + "CCC"))  # in frame, but a stop codon in the linker
    assert len(stop) == 1 and stop[0]["in_frame"] and stop[0]["stop_at"] == 300 + len(sp) + 1


def test_methylation_flags_only_context_made_overlaps():
    from plasmidbank.tools.sites import methylation_context

    # Addgene's example: XbaI (TCTAGA) preceded by GA makes GATC across the site's edge
    assert [m["motif"] for m in methylation_context("AAAAAGATCTAGAAAAAA", False, "TCTAGA")] == ["Dam"]
    assert [m["motif"] for m in methylation_context("AAAAATCTAGATCAAAAA", False, "TCTAGA")] == ["Dam"]
    assert methylation_context("AAAAACCTCTAGAAAAAA", False, "TCTAGA") == []
    # a motif entirely inside the site is the enzyme's own property: not flagged here
    assert methylation_context("AAAAAGGATCCAAAAAA", False, "GGATCC") == []
    # Dcm CCWGG, and the check wraps the origin of a circular map
    assert [m["motif"] for m in methylation_context("TGGAAAAAAAAAAAAAGAATTCC", True, "GAATTC")] == ["Dcm"]
