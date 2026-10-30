"""Annotation library: common + private lists, exact matching on both strands, seeding from the bank."""

import random

import pytest

from plasmidbank.annotations import Library, annotate, suggest_from_bank
from plasmidbank.formats import Feature, MapRecord
from plasmidbank.sequtil import revcomp

PART = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG"  # a 39-nt "GFP start"


def _seq(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def _map(seq, pid="PMC-00001", features=(), circular=True):
    return MapRecord(f"/b/{pid}.gb", pid, "PMC", "filename", [], "genbank", "ok", seq=seq,
                     topology="circular" if circular else "linear", features=list(features))


def test_common_is_shared_and_private_is_private(tmp_path):
    lib = Library(tmp_path)
    c = lib.add("ann", "common", "GFP start", "CDS", PART.lower())
    p = lib.add("ann", "mine", "my tag", "misc_feature", "CATCACCATCACCATCAC")
    assert [(e["label"], e["scope"]) for e in lib.entries("ann")] == [("GFP start", "common"), ("my tag", "mine")]
    assert [e["label"] for e in lib.entries("bob")] == ["GFP start"]
    lib.update("bob", c, label="GFP N-term")  # anyone edits the common list
    assert lib.entries("ann")[0]["label"] == "GFP N-term"
    with pytest.raises(LookupError):
        lib.delete("bob", p)  # but not someone else's private entry
    with pytest.raises(ValueError, match="already"):
        lib.add("bob", "common", "GFP N-term", "CDS", PART)
    with pytest.raises(ValueError, match="DNA"):
        lib.add("ann", "mine", "x", "CDS", "ACGTJ" * 5)
    with pytest.raises(ValueError, match="long"):
        lib.add("ann", "mine", "x", "CDS", "ACGT")
    lib.delete("ann", p)
    assert len(lib.entries("ann")) == 1
    assert not lib.on_upload("ann")
    lib.set_on_upload("ann", True)
    assert lib.on_upload("ann") and not lib.on_upload("bob")


def test_annotate_finds_both_strands_and_across_the_origin(tmp_path):
    lib = Library(tmp_path)
    lib.add("ann", "common", "GFP start", "CDS", PART)
    s = _seq(3000, 1)
    s = PART[20:] + s[19:1000] + PART + s[1000 + len(PART):2000] + revcomp(PART) + s[2000 + len(PART):]
    s = s[: 3000 - 20] + PART[:20]  # the third copy wraps the origin
    m = _map(s)
    got = [(p["start"], p["end"], p["strand"]) for p in annotate(m, lib.entries("ann"))]
    assert (1001, 1000 + len(PART), 1) in got
    assert (2001, 2000 + len(PART), -1) in got
    assert (2981, len(PART) - 20, 1) in got  # wraps: starts near the end, ends near the start
    # an existing annotation with the same label and span is not proposed again
    m.features = [Feature("CDS", "GFP start", 1000, 1000 + len(PART), 1)]
    assert (1001, 1000 + len(PART), 1) not in [(p["start"], p["end"], p["strand"]) for p in annotate(m, lib.entries("ann"))]


def test_suggestions_rank_parts_by_how_many_maps_carry_them():
    maps = []
    for i in range(4):
        s = _seq(2000, 10 + i)
        s = s[:500] + PART + s[500 + len(PART):]
        feats = [Feature("CDS", "GFP start", 500, 500 + len(PART), 1), Feature("misc_feature", "8938", 100, 140, 1)]
        maps.append(_map(s, f"PMC-0000{i}", feats))
    sug = suggest_from_bank(maps, known=set(), min_maps=3)
    assert [(c["label"], c["maps"], c["seq"]) for c in sug] == [("GFP start", 4, PART)]  # '8938' is not a part name
    assert suggest_from_bank(maps, known={PART}, min_maps=3) == []
    assert suggest_from_bank(maps, known=set(), min_maps=5) == []


def test_suggested_labels_are_unescaped():
    maps = []
    for i in range(3):
        s = _seq(1500, 30 + i)
        s = s[:300] + PART + s[300 + len(PART):]
        maps.append(_map(s, f"PMC-0001{i}", [Feature("promoter", "LjUBI1\\promoter\\(no\\BsaI)", 300, 300 + len(PART), 1)]))
    assert suggest_from_bank(maps, set(), 3)[0]["label"] == "LjUBI1 promoter (no BsaI)"


def test_a_shorter_variant_inside_a_longer_one_is_dropped(tmp_path):
    lib = Library(tmp_path)
    lib.add("ann", "common", "ori", "rep_origin", PART)
    lib.add("ann", "common", "ori", "rep_origin", PART[5:30])  # a short variant of the same part
    s = _seq(2000, 3)
    m = _map(s[:500] + PART + s[500 + len(PART):])
    assert [(p["start"], p["length"]) for p in annotate(m, lib.entries("ann"))] == [(501, len(PART))]


def test_import_third_party_list(tmp_path):
    """pLannotate-style FASTA + CSV: lab names win, junk and repeats are skipped, dry run writes nothing."""
    from plasmidbank.annotation_import import main

    other = _seq(40, 7)
    (tmp_path / "f.fa").write_text(f">GFP_(2)\n{PART[:20]}\n{PART[20:]}\n>CmR\n{other}\n>CmR_copy\n{other}\n"
                                   f">short\nACGT\n>oriX\n{_seq(30, 8)}\n>nometa\n{_seq(25, 9)}\n")
    (tmp_path / "m.csv").write_text("sseqid,Feature,Type,Description\nGFP_(2),GFP,CDS,green\n"
                                    "CmR,CmR,CDS,chloramphenicol\nCmR_copy,CmR,CDS,dup\noriX,ori X,Origin,odd type\n")
    lib = Library(tmp_path / "s")
    lib.add("ann", "common", "GFP start", "CDS", PART)  # already in the lab's list under its own name
    args = [str(tmp_path / "f.fa"), "--meta", str(tmp_path / "m.csv"), "--state", str(tmp_path / "s"), "--source", "pLannotate"]
    main(args)
    assert len(lib.entries("x")) == 1
    main(args + ["--write"])
    got = {e["label"]: (e["type"], e["note"]) for e in lib.entries("x")}
    assert got == {"GFP start": ("CDS", ""), "CmR": ("CDS", "pLannotate; chloramphenicol"),
                   "ori X": ("rep_origin", "pLannotate; odd type; type Origin"), "nometa": ("misc_feature", "pLannotate")}
