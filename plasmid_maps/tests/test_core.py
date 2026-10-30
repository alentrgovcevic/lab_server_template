import random
from pathlib import Path

import pytest

from plasmidbank.bank import Bank
from plasmidbank.formats import Feature, MapRecord, read_map, sniff
from plasmidbank.ids import identify
from plasmidbank.oligos import Oligo, normalise
from plasmidbank.sequtil import _least_rotation, canonical_hash, revcomp
from plasmidbank.tools.align import Query, load_queries, mott_trim, verify
from plasmidbank.tools.primers import group_oligos, match_map

SAMPLE = Path(__file__).resolve().parents[2] / "examples" / "bank"


# ── IDs, from real filenames ──────────────────────────────────────────────────


@pytest.mark.parametrize(
    "name,pid,xrefs",
    [
        ("PMC-5832.gb", "PMC-05832", []),
        ("27ABCDEF_1234567_PMC-06982.gb", "PMC-06982", []),
        ("PMC-03456-pBIN19-p35S(long)-GUS-tNOS.gb", "PMC-03456", []),
        (
            "PMC-07001_pCAMBIA1301_binary_vector.gbk",
            "PMC-07001",
            [],
        ),
        ("PMC-03405GeneX PU in PMC-02483_x.gb", "PMC-03405", ["PMC-02483"]),
        ("PMC-06841_EC55699_x.gb", "PMC-06841", ["EC55699"]),
        ("---EC25800 in EC41308.seq", "EC25800", ["EC41308"]),
        ("EC25111 pL0M-PU-promoterA-25111.gb", "EC25111", []),
        ("EC31476-pl1m-f1-p35s-genex-egfp-t35s.gb", "EC31476", []),
    ],
)
def test_identify(name, pid, xrefs):
    r = identify(name)
    assert (r.plasmid_id, r.source, r.xrefs) == (pid, "filename", xrefs)


def test_identify_locus_fallback_and_none():
    assert identify("pUC19.gb", "PMC-5993").plasmid_id == "PMC-05993"
    assert identify("pUC19.gb", "PMC-5993").source == "locus"
    assert identify("pUC19.gb", "Exported").kind == "NONE"
    assert identify("SPEC12345.gb").kind == "NONE"  # EC needs a word boundary


# ── formats over every sample file ────────────────────────────────────────────


def test_every_sample_file_reads():
    bank = Bank.load(SAMPLE)
    assert len(bank.maps) >= 10
    assert all(m.status != "error" for m in bank.maps), [
        m.path for m in bank.maps if m.status == "error"
    ]
    assert all(m.plasmid_id for m in bank.maps)
    for m in bank.maps:
        if m.format != "clc":
            assert m.seq and m.content_hash


def test_clc_pairs_have_parsable_genbank_twin():
    for clc in (SAMPLE / "clc_gb_pairs").glob("*.clc"):
        assert read_map(clc).status == "clc_binary"
        twin = read_map(clc.with_suffix(".gbk"))
        assert twin.plasmid_id == read_map(clc).plasmid_id and twin.length > 0


def test_sniff():
    assert sniff(b"CLC\x00\x00") == "clc"
    assert sniff(b"  LOCUS x") == "genbank"
    assert sniff(b">a\nACGT") == "fasta"
    assert sniff(b"ACGTN\nACGT\n") == "raw"
    assert sniff(b"hello world") == "unknown"


def test_non_utf8_file_decodes(tmp_path):
    p = tmp_path / "PMC-00001.gb"
    p.write_bytes(
        b"LOCUS       x  20 bp DNA circular UNA 01-JAN-2000\nDEFINITION  \xb4tick.\n"
        b"ORIGIN\n        1 acgtacgtac gtacgtacgt\n//\n"
    )
    m = read_map(p)
    assert m.length == 20 and m.encoding == "cp1252"


# ── sequence helpers ─────────────────────────────────────────────────────────


def test_least_rotation_matches_brute_force():
    rnd = random.Random(0)
    for _ in range(300):
        s = "".join(rnd.choice("AC") for _ in range(rnd.randint(1, 12)))
        brute = min(s[i:] + s[:i] for i in range(len(s)))
        k = _least_rotation(s)
        assert s[k:] + s[:k] == brute


def test_canonical_hash_rotation_and_strand_invariant():
    rnd = random.Random(1)
    s = "".join(rnd.choice("ACGT") for _ in range(500))
    h = canonical_hash(s, True)
    assert canonical_hash(s[123:] + s[:123], True) == h
    assert canonical_hash(revcomp(s[77:] + s[:77]), True) == h
    assert canonical_hash(s, False) != canonical_hash(s[1:] + s[:1], False)


# ── oligos ───────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "raw,seq,flags",
    [
        ("gaTTGTTCAG tacaaa", "GATTGTTCAGTACAAA", set()),
        ("[Phos]ACGTACGTAC", "ACGTACGTAC", {"modified"}),
        ("ACG[I]ACGTACGT", "ACGNACGTACGT", {"inosine", "degenerate"}),
        ("(6-FAM)gcgattccgg(BHQ1)", "GCGATTCCGG", {"modified"}),
        ("use marker 71R", "", {"nonstandard", "not_a_sequence", "empty"}),
        ("", "", {"empty"}),
    ],
)
def test_oligo_normalise(raw, seq, flags):
    assert normalise(raw) == (seq, flags)


# ── clone verification ───────────────────────────────────────────────────────


def _random_map(n=3000, seed=2, circular=True, features=()):
    rnd = random.Random(seed)
    s = "".join(rnd.choice("ACGT") for _ in range(n))
    return MapRecord(
        "x.gb",
        "PMC-00001",
        "PMC",
        "filename",
        [],
        "genbank",
        "ok",
        seq=s,
        topology="circular" if circular else "linear",
        features=list(features),
    )


def _sub(base):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


@pytest.mark.parametrize("reverse", [False, True])
def test_read_across_origin_reports_substitutions(reverse):
    ref = _random_map(features=[Feature("rep_origin", "ori", 2800, 100, 1)])
    s = ref.seq
    read = list(s[-200:] + s[:300])  # spans the origin
    read[50] = _sub(read[50])  # map position 2851 (1-based)
    read[260] = _sub(read[260])  # map position 61
    read = "".join(read)
    if reverse:
        read = revcomp(read)
    r = verify(ref, Query("r", read))
    assert r.mode == "read" and r.strand == ("-" if reverse else "+")
    assert sorted(d.ref_pos for d in r.discrepancies) == [61, 2851]
    assert all(d.kind == "substitution" for d in r.discrepancies)
    assert r.ref_start == 2801 and r.ref_end == 300
    assert any(
        "ori" in f for d in r.discrepancies if d.ref_pos == 2851 for f in d.features
    )


def test_read_with_n_is_not_a_difference():
    ref = _random_map()
    read = ref.seq[1000:1100] + "N" + ref.seq[1101:1500]
    r = verify(ref, Query("r", read))
    assert r.discrepancies == [] and r.verdict == "perfect match"


def test_whole_plasmid_rotated_revcomp_with_deletion():
    ref = _random_map(n=5000, seed=3)
    s = ref.seq
    cons = s[:2000] + s[2001:]  # delete map base 2001
    cons = revcomp(cons[3100:] + cons[:3100])
    r = verify(ref, Query("c", cons))
    assert r.mode == "whole" and r.strand == "-"
    assert [(d.kind, d.ref_pos) for d in r.discrepancies] == [("deletion", 2001)]


def test_forced_whole_mode_on_a_short_read_checks_it_as_a_read():
    # several Sanger reads with "whole plasmid" chosen all failed
    ref = _random_map(n=5000, seed=5)
    read = list(ref.seq[1000:1800])
    read[300] = _sub(read[300])  # map position 1301
    r = verify(ref, Query("r", "".join(read)), "whole")
    assert r.mode == "read" and any("too short" in w for w in r.warnings)
    assert [(d.kind, d.ref_pos) for d in r.discrepancies] == [("substitution", 1301)]


def test_search_puts_the_map_before_maps_that_only_name_it(tmp_path):
    seq = "ACGT" * 300
    for name in [f"EC4{i:04d}_built_from_PMC-01118.gb" for i in range(40)] + ["PMC-01118_pICH41414.gb"]:
        (tmp_path / name).write_text(f">x\n{seq}\n")
    bank = Bank.load(tmp_path)
    hits = bank.search("PMC-01118")
    assert len(hits) == 41 and bank.maps[hits[0]].plasmid_id == "PMC-01118"


def test_whole_plasmid_homopolymer_indel_flagged():
    ref = _random_map(n=4000, seed=4)
    s = ref.seq[:1500] + "AAAAAA" + ref.seq[1506:]
    ref.seq = s
    cons = s[:1500] + "AAAAAAA" + s[1506:]  # one extra A in the run
    r = verify(ref, Query("c", cons))
    assert len(r.discrepancies) == 1
    d = r.discrepancies[0]
    assert d.kind == "insertion" and d.homopolymer
    assert (
        r.verdict == "perfect match"
    )  # homopolymer indels in whole-plasmid mode are not counted


def test_unrelated_read_is_no_match():
    ref = _random_map(seed=5)
    other = _random_map(seed=6).seq[:600]
    assert verify(ref, Query("x", other)).verdict == "no match"


def test_mott_trim_removes_bad_ends():
    qual = [5] * 30 + [40] * 500 + [8] * 40
    assert mott_trim(qual) == (30, 530)


def test_load_queries_fastq_trims():
    fq = b"@r1\n" + b"ACGT" * 50 + b"\n+\n" + b"#" * 20 + b"I" * 160 + b"#" * 20 + b"\n"
    q = load_queries("r.fastq", fq)[0]
    assert len(q.seq) == 160 and q.trimmed == (20, 20)


def test_real_map_against_its_own_sequence():
    m = read_map(next((SAMPLE / "pmc_gb").glob("PMC-07001*")))
    r = verify(m, Query("self", m.seq[100:900]))
    assert r.verdict == "perfect match" and r.ref_start == 101 and r.ref_end == 900


# ── primer matching ──────────────────────────────────────────────────────────


def _oligo(i, seq):
    return Oligo(f"PM-{i:05d}", f"p{i}", "", "", "", seq, seq, "", "")


def test_primer_matching_strands_tails_and_origin():
    ref = _random_map(n=3000, seed=7)
    s = ref.seq
    fwd = "GGTCTCAAGT" + s[100:122]  # 10-nt Golden Gate tail + 22 nt binding
    rev = revcomp(s[2990:] + s[:12])  # binds across the origin on the bottom strand
    mm = list(s[500:520])
    mm[3] = _sub(mm[3])  # one 5'-side mismatch
    groups, skipped = group_oligos(
        [
            _oligo(1, fwd),
            _oligo(2, rev),
            _oligo(3, "".join(mm)),
            _oligo(4, "ACGT"),
            _oligo(5, s[700:720] + "N"),
        ]
    )
    assert skipped == {"shorter than 15 nt": 1, "degenerate 3' end": 1}
    hits = {h.oligo_ids[0]: h for h in match_map(ref, groups)}
    # with mismatches allowed, a tail base can match the template by chance; check tails exactly
    exact = {h.oligo_ids[0]: h for h in match_map(ref, groups, max_mismatches=0)}
    f = exact["PM-00001"]
    assert (f.strand, f.start, f.end, f.three_prime, f.tail_length, f.mismatches) == (
        "+",
        101,
        122,
        122,
        10,
        0,
    )
    r = hits["PM-00002"]
    assert (r.strand, r.start, r.end, r.three_prime) == ("-", 2991, 12, 2991)
    assert hits["PM-00003"].mismatches == 1 and hits["PM-00003"].binding_length == 20
    assert f.tm is not None


def test_primers_on_real_map_find_amp_primer():
    m = read_map(next((SAMPLE / "pmc_gb").glob("PMC-07001*")))
    oligo = _oligo(9, m.seq[1000:1022])
    groups, _ = group_oligos([oligo])
    hits = match_map(m, groups)
    assert [(h.start, h.end) for h in hits] == [(1001, 1022)]


# ── restriction sites ────────────────────────────────────────────────────────

def test_sites_across_origin_and_sets():
    from plasmidbank.tools.sites import find_sites
    ref = _random_map(n=2000, seed=11)
    clean = ref.seq.replace("GAATTC", "GAATTA").replace("GGTCTC", "GGTCTA").replace("GAGACC", "GAGACA")
    ref.seq = "AATTC" + clean[5:-1] + "G"  # EcoRI site G^AATTC spans the origin
    hits = {h["enzyme"]: h for h in find_sites(ref, "common", "EcoRI")}
    assert hits["EcoRI"]["cuts"] == [1]
    ref.topology = "linear"
    assert find_sites(ref, "common", "EcoRI") == []
    with pytest.raises(ValueError):
        find_sites(ref, "common", "NotAnEnzyme")


def test_unique_set_only_single_cutters():
    from plasmidbank.tools.sites import find_sites
    m = read_map(next((SAMPLE / "pmc_gb").glob("PMC-07001*")))
    uniq = find_sites(m, "unique")
    assert uniq and all(h["n"] == 1 for h in uniq)
    bsa = {h["enzyme"]: h for h in find_sites(m, "goldengate")}["BsaI"]
    assert bsa["n"] == m.seq.count("GGTCTC") + m.seq.count("GAGACC")


# ── primer design ────────────────────────────────────────────────────────────

def test_cloning_primers_start_exactly_at_region_ends_and_follow_strand():
    from plasmidbank.tools.design import design_cloning
    ref = _random_map(n=3000, seed=21)
    s = ref.seq
    p = design_cloning(ref, 1001, 1600, [], tail_f="TTGGTCTCAAATG", tail_r="TTGGTCTCAAAGC")["pairs"][0]
    f, r = p["forward"], p["reverse"]
    assert f["start"] == 1001 and s[1000:1000 + f["length"]] == f["sequence"]
    assert r["end"] == 1600 and revcomp(s[1600 - r["length"]:1600]) == r["sequence"]
    assert f["full"].startswith("TTGGTCTCAAATG") and p["product"] == 600 + 26
    q = design_cloning(ref, 1001, 1600, [], tail_f="AAAA", sense="-")["pairs"][0]
    assert q["forward"]["strand"] == "-" and q["forward"]["end"] == 1600 and q["forward"]["tail"] == "AAAA"


def test_primer_length_range_is_respected():
    from plasmidbank.tools.design import design_cloning, design_pcr
    ref = _random_map(n=3000, seed=21)
    for pr in design_pcr(ref, 1200, 1500, [], size=(24, 26))["pairs"]:
        assert 24 <= pr["forward"]["length"] <= 26 and 24 <= pr["reverse"]["length"] <= 26
    p = design_cloning(ref, 1001, 1600, [], size=(30, 32))["pairs"][0]
    assert 30 <= p["forward"]["length"] <= 32 and 30 <= p["reverse"]["length"] <= 32


def test_pcr_pairs_flank_region_across_origin():
    from plasmidbank.tools.design import design_pcr
    ref = _random_map(n=4000, seed=22)
    r = design_pcr(ref, 3900, 150, [], max_extra=200)
    assert r["pairs"]
    for pr in r["pairs"]:
        f, rv = pr["forward"], pr["reverse"]
        assert f["strand"] == "+" and rv["strand"] == "-"
        # forward 3' end before the region start, reverse 3' end after the region end (modulo L)
        assert 0 < (3900 - f["three_prime"]) % 4000 <= 230
        assert 0 < (rv["three_prime"] - 150) % 4000 <= 230
        assert ref.seq[f["start"] - 1:f["start"] - 1 + f["length"]] == f["sequence"]


def test_sequencing_reuses_collection_primer():
    from plasmidbank.tools.design import design_sequencing
    ref = _random_map(n=3000, seed=23)
    s = ref.seq
    groups, _ = group_oligos([_oligo(1, s[880:902]), _oligo(2, revcomp(s[1700:1722]))])
    r = design_sequencing(ref, 1001, 1500, match_map(ref, groups))
    assert [x["ids"] for x in r["existing_forward"]] == [["PM-00001"]]
    assert r["existing_forward"][0]["distance"] == 1001 - 902
    assert [x["ids"] for x in r["existing_reverse"]] == [["PM-00002"]]
    assert r["new"] == {"forward": [], "reverse": []}


def test_golden_gate_cloning_trims_overlap_and_scans_product():
    from plasmidbank.tools.design import design_cloning
    ref = _random_map(n=3000, seed=24)
    body = ref.seq.replace("GAAGAC", "GAAGAA").replace("GTCTTC", "GTCTTA")
    cds = "ATG" + body[1003:1597] + "TAA"  # 600 bp CDS at map 1001..1600, plus strand
    ref.seq = body[:1000] + cds + body[1600:]
    r = design_cloning(ref, 1001, 1600, [], enzyme="BpiI", oh5="AATG", oh3="GCTT")
    assert r["trimmed"] == {"5prime": 3, "3prime": 0}  # ATG supplied by the AATG fusion site
    p = r["pairs"][0]
    assert p["forward"]["full"].startswith("TTGAAGACAAAATG") and p["forward"]["start"] == 1004
    assert p["reverse"]["full"].startswith("TTGAAGACAA" + "AAGC")
    assert p["product_seq"] == "TTGAAGACAA" + "AATG" + cds[3:] + "GCTT" + "TTGTCTTCAA"
    assert p["internal_sites"] == [] and {s["enzyme"] for s in p["tail_sites"]} == {"BpiI"}


def test_minus_strand_cloning_product_reads_on_part_sense():
    from plasmidbank.tools.design import design_cloning
    ref = _random_map(n=3000, seed=25)
    p = design_cloning(ref, 1001, 1600, [], sense="-")["pairs"][0]
    assert p["product_seq"] == revcomp(ref.seq[1000:1600])


# ── GenBank cleaner ──────────────────────────────────────────────────────────

def _gb_record(features):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    rnd = random.Random(31)
    rec = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(3000))), id="PMC-00001", name="PMC-00001",
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    for ftype, s, e, strand, quals in features:
        rec.features.append(SeqFeature(SimpleLocation(s, e, strand), type=ftype, qualifiers=quals))
    return rec


def test_cleaner_merges_stacks_and_keeps_provenance():
    from plasmidbank.tools.clean import ALL_RULES, clean_record
    rec = _gb_record([
        ("misc_feature", 100, 900, -1, {"label": ["Insert"]}),
        ("misc_feature", 100, 900, -1, {"label": ["Insert Sequence"]}),       # generic stack -> merged
        ("misc_feature", 96, 896, 1, {"label": ["Insert"]}),                  # ±5 bp, other strand -> flagged only
        ("CDS", 1000, 1300, 1, {"label": ["Amp R"]}),                         # synonym
        ("CDS", 1000, 1300, 1, {"label": ["Amp R"]}),                         # exact duplicate
        ("misc_feature", 1000, 1300, 1, {}),                                  # unlabelled under a label
        ("Vector", 1400, 2990, 1, {"label": ["pIV10"], "vntifkey": ["21"]}),  # key map, vector
        ("misc_feature", 2000, 2100, 1, {"label": ["Fragment (PCR (a1 - a2))"]}),
        ("misc_feature", 2200, 2300, 1, {"label": ["LjUBI1\\promoter\\(no\\BsaI)"]}),
        ("primer_bind", 2400, 2420, 1, {"label": ["yfp1 (mismatches: 3)"]}),
        ("primer_bind", 2500, 2520, 1, {"label": ["PM-00042 (mismatches: 2)"]}),
    ])
    seq = str(rec.seq)
    r = clean_record(rec, ALL_RULES, lab_oligo_names={"PM-00042"})
    labels = sorted(f.qualifiers.get("label", ["?"])[0] for f in r.record.features)
    assert str(r.record.seq) == seq
    assert "AmpR" in labels and labels.count("AmpR") == 1
    assert "LjUBI1 promoter (no BsaI)" in labels
    assert "pIV10" not in labels and "backbone: pIV10" in r.record.annotations["comment"]
    assert "assembled from: PCR (a1 - a2)" in r.record.annotations["comment"]
    assert "yfp1 (mismatches: 3)" not in labels and "PM-00042 (mismatches: 2)" in labels
    stacked = [f for f in r.record.features if int(f.location.start) == 100]
    assert len(stacked) == 1 and any("Insert" in n for n in stacked[0].qualifiers.get("note", []))
    assert any("near-duplicate" in c.detail for c in r.flags)
    assert any(f.type == "misc_feature" and int(f.location.start) == 96 for f in r.record.features)


def test_cleaner_defaults_are_conservative_and_idempotent():
    import io
    from Bio import SeqIO
    from plasmidbank.tools.clean import clean_record, read_genbank
    for path in sorted((SAMPLE / "pmc_gb").glob("*.gb*"))[:20]:
        first = clean_record(read_genbank(path))
        again = clean_record(SeqIO.read(io.StringIO(first.genbank()), "genbank"))
        assert again.changes == []
        assert first.n_before - first.n_after <= 5  # defaults never strip a map bare
        assert all(c.rule not in ("fragments_to_comment", "vector_to_comment", "clc_primer_hits", "drop_colours")
                   for c in first.changes)


def test_cleaner_flags_cds_problems_without_fixing():
    from plasmidbank.tools.clean import clean_record
    rec = _gb_record([("CDS", 100, 401, 1, {"label": ["broken"]})])  # 301 bp: not a multiple of 3
    r = clean_record(rec)
    assert any("not a multiple of 3" in c.detail for c in r.flags)
    assert len(r.record.features) == 1


@pytest.mark.parametrize("key, expected", [
    ("Misc._feature", "misc_feature"),   # CLC
    ("Misc.", "misc_feature"),           # CLC
    ("Gene", "gene"),                    # CLC, case only
    ("Replication_ori", "rep_origin"),   # CLC
    ("Probe", "misc_feature"),           # CLC
    ("mutation", "variation"),           # Vector NTI
    ("Seq-region", "misc_feature"),      # ApE
    ("Module", "misc_feature"),          # hand-made
    ("promoter", "promoter"),            # legacy key, left alone
])
def test_cleaner_maps_every_unknown_key(key, expected):
    from plasmidbank.tools.clean import KNOWN_KEYS, clean_record
    r = clean_record(_gb_record([(key, 100, 400, 1, {})]))
    f = r.record.features[0]
    assert f.type == expected and f.type in KNOWN_KEYS
    if key != expected:
        assert f"original key: {key}" in f.qualifiers["note"]
        assert f.qualifiers.get("label") or f.qualifiers.get("gene")  # unlabelled keeps its old key as label


def test_cleaner_fixes_clc_quotes_dates_and_unreadable_sites():
    from Bio.SeqFeature import SeqFeature
    from plasmidbank.tools.clean import clean_record
    rec = _gb_record([("misc_feature", 100, 400, 1, {"label": ["35S promoter", '"35S promoter"']})])
    rec.annotations["date"] = "14-okt-2008"
    rec.features.append(SeqFeature(None, type="Site"))  # what Biopython makes of CLC's 'Site -1^1'
    r = clean_record(rec)
    assert r.record.features[0].qualifiers["label"] == ["35S promoter"]
    assert len(r.record.features) == 1 and any("unreadable" in c.detail for c in r.changes)
    assert "14-OCT-2008" in r.genbank().splitlines()[0]


@pytest.mark.parametrize("path", sorted((SAMPLE / "clc_gb_pairs").glob("*.gbk")), ids=lambda p: p.stem)
def test_cleaner_on_clc_exports(path):
    import io
    from Bio import SeqIO
    from plasmidbank.tools.clean import KNOWN_KEYS, clean_record, read_genbank
    before = read_genbank(path)
    first = clean_record(read_genbank(path))
    out = SeqIO.read(io.StringIO(first.genbank()), "genbank")
    assert str(out.seq) == str(before.seq)
    assert all(f.type in KNOWN_KEYS for f in out.features)
    dropped = sum(c.action == "dropped" for c in first.changes)
    assert len(out.features) == len(before.features) - dropped  # every loss is in the report
    assert out.annotations["date"] != "01-JAN-1980"
    assert clean_record(out).changes == []  # idempotent


def test_batch_cleaner_never_writes_in_place(tmp_path):
    import shutil
    from plasmidbank.clean import main
    src = tmp_path / "in"
    src.mkdir()
    one = next((SAMPLE / "pmc_gb").glob("PMC-07001*"))
    shutil.copy(one, src / "PMC-07001.gb")
    before = (src / "PMC-07001.gb").read_bytes()
    with pytest.raises(SystemExit):
        main([str(src), "--out", str(src)])
    assert main([str(src), "--out", str(tmp_path / "out")]) == 0
    assert (src / "PMC-07001.gb").read_bytes() == before
    assert (tmp_path / "out" / "PMC-07001.gb").exists() and (tmp_path / "out" / "clean_report.html").exists()


def _fake_ab1_query(seq):
    """A read with a synthetic chromatogram: the called base peaks high, the others stay low."""
    from plasmidbank.tools.align import Query

    peaks = [10 * i + 5 for i in range(len(seq))]
    trace = {b: [50] * (10 * len(seq)) for b in "ACGT"}
    for i, b in enumerate(seq):
        trace[b][peaks[i]] = 1000
    return Query("r", seq, [40] * len(seq), (0, 0), "abi", trace, peaks)


def test_trace_window_reads_like_the_map_on_both_strands():
    from plasmidbank.tools.align import trace_window

    q = _fake_ab1_query("ACGTTGCAAGCT")
    w = trace_window(q, "+", 5, flank=2)  # 5th base, a T
    assert w["bases"] == "GTTGC" and w["bases"][w["focus"]] == "T"
    top = max("ACGT", key=lambda b: w["channels"][b][w["peaks"][w["focus"]]])
    assert top == "T"
    # aligned on the minus strand: query_pos counts along the reverse complement
    w = trace_window(q, "-", 8, flank=2)  # oriented read pos 8 = original index 4 (T), shown as A
    assert w["bases"][w["focus"]] == "A"
    top = max("ACGT", key=lambda b: w["channels"][b][w["peaks"][w["focus"]]])
    assert top == "A"
    assert trace_window(type(q)("x", "ACGT"), "+", 1) is None  # no peak table: no window


def _blockmap(r, oq, L):
    """Read index -> map position, from VerifyResult.blocks (what the chromatogram viewer does)."""
    pos = {}
    for rp, qp, n in r.blocks:
        for k in range(n):
            pos[(qp - 1 + k) % len(oq)] = (rp - 1 + k) % L + 1
    return pos


@pytest.mark.parametrize("mode,reverse", [("read", False), ("read", True), ("whole", False), ("whole", True)])
def test_blocks_map_every_aligned_read_base_to_the_map(mode, reverse):
    from plasmidbank.tools.align import Query, verify

    m = _random_map(1500, 7)
    s = m.seq
    if mode == "read":
        clone = s[-150:] + s[:200] + s[203:400]  # across the origin, 3 bp deleted
    else:
        clone = s[600:] + s[:600]  # rotated consensus
        clone = clone[:300] + "T" + clone[300:]  # one inserted base
    oq = clone
    q = Query("q", revcomp(clone) if reverse else clone)
    r = verify(m, q, mode)
    assert r.strand == ("-" if reverse else "+")
    pos = _blockmap(r, oq, m.length)
    assert len(pos) >= len(oq) - 1 and all(s[p - 1] == oq[i] for i, p in pos.items())
    if mode == "read":
        assert pos[149] == m.length and pos[150] == 1 and pos[349] == 200 and pos[350] == 204


def test_trace_full_reads_like_the_map_on_the_minus_strand():
    from plasmidbank.tools.align import trace_full

    q = _fake_ab1_query("AACGTTGCAAGCTT")
    t = trace_full(q, "-")
    assert t["bases"] == revcomp(q.seq) and t["qual"] == q.qual[::-1] and len(t["peaks"]) == len(q.seq)
    for i, b in enumerate(t["bases"]):
        assert max("ACGT", key=lambda c: t["channels"][c][t["peaks"][i]]) == b
    assert trace_full(type(q)("x", "ACGT"), "+") is None


def test_read_support_tells_shared_changes_from_single_read_ones():
    from plasmidbank.tools.align import Query, read_support, verify

    m = _random_map(1200, 41)
    ref = m.seq
    mut = ref[:300] + _sub(ref[300]) + ref[301:]
    reads = [Query("r1", mut[100:500]), Query("r2", mut[250:700]), Query("r3", ref[250:600])]
    res = [verify(m, q, "read") for q in reads]
    sup = read_support(res, len(ref))
    d1 = [(d.ref_pos, s) for d, s in zip(res[0].discrepancies, sup[0])]
    assert d1 == [(301, {"covered_by": 3, "seen_in": 2})]
    assert res[2].discrepancies == []


# ── site-directed mutagenesis (Q5 back-to-back) ──────────────────────────────

def _kld(ref, r):
    """Simulate whole-plasmid PCR with the two primers, then ligation: the circular product."""
    L, s2 = ref.length, ref.seq + ref.seq
    f, rv = r["forward"], r["reverse"]
    a = f["start"] - 1 + len(f["sequence"])  # first template base after the forward binding part
    b = rv["start"] - 1  # first base of the reverse binding part
    middle = s2[a : b if b >= a else b + L]
    return f["full"] + middle + revcomp(rv["full"])


@pytest.mark.parametrize("start,end,new", [
    (1500, 1500, "T"),            # point substitution
    (1500, 1502, ""),             # 3-bp deletion
    (1500, 1499, "GGATCC"),       # 6-nt insertion: all on the forward primer
    (1500, 1499, "CATCACCATCACCATCAC"),  # His6 insertion: split across both primers
    (2, 4, "AAAAAA"),             # replacement right at the origin
    (3000, 3000, "G"),            # last base of the map
])
def test_sdm_primers_rebuild_the_mutant_plasmid(start, end, new):
    from plasmidbank.tools.mutagenesis import design_sdm
    ref = _random_map(n=3000, seed=31)
    if ref.seq[start - 1:end] == new:
        new = _sub(new)
    r = design_sdm(ref, start, end, new, [])
    product = _kld(ref, r)
    assert canonical_hash(product, True) == canonical_hash(r["mutant_seq"], True)
    assert len(r["mutant_seq"]) == 3000 - (end - start + 1) + len(new)
    for p in (r["forward"], r["reverse"]):
        assert 18 <= p["length"] <= 35 and not set(p["full"]) - set("ACGT")
    if len(new) > 6:
        assert r["reverse"]["tail"] and r["forward"]["tail"]


def _cds_map(strand):
    from plasmidbank.formats import Feature
    ref = _random_map(n=3000, seed=32)
    # a CDS from codons with no stop: M + 99 codons, then a stop
    rnd = random.Random(5)
    sense = [c for c in (a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT") if c not in ("TAA", "TAG", "TGA")]
    cds = "ATG" + "".join(rnd.choice(sense) for _ in range(99)) + "TAA"
    s = ref.seq[:1000] + (cds if strand > 0 else revcomp(cds)) + ref.seq[1000 + len(cds):]
    ref.seq = s
    ref.features = [Feature("CDS", "gene", 1000, 1000 + len(cds), strand)]
    return ref


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("change", ["sub", "stop", "del", "ins"])
def test_protein_changes_give_the_intended_protein(strand, change):
    from Bio.Seq import Seq
    from plasmidbank.tools.mutagenesis import cds_positions, design_sdm, protein_change
    ref = _cds_map(strand)
    f = ref.features[0]
    cds = "".join(ref.seq[p] for p in cds_positions(ref, f))
    cds = cds if strand > 0 else cds.translate(str.maketrans("ACGT", "TGCA"))
    wt = str(Seq(cds).translate())
    a, b = wt[47], wt[48]
    to = "W" if a != "W" else "R"
    note = {"sub": f"{a}48{to}", "stop": f"{a}48*", "del": f"{a}48_{b}49del", "ins": f"{a}48_{b}49insGSG"}[change]
    pc = protein_change(ref, f, note)
    expect = {"sub": wt[:47] + to + wt[48:], "stop": wt[:47] + "*" + wt[48:], "del": wt[:47] + wt[49:],
              "ins": wt[:48] + "GSG" + wt[48:]}[change]
    assert pc["protein_after"] == expect
    r = design_sdm(ref, pc["start"], pc["end"], pc["new"], [])
    mutant = r["mutant_seq"]
    # translate the mutant CDS: same start, length changed by the indel
    d = len(mutant) - len(ref.seq)
    if strand > 0:
        mcds = mutant[1000:1000 + len(cds) + d]
    else:
        mcds = revcomp(mutant[1000:1000 + len(cds) + d])
    assert str(Seq(mcds).translate()) == expect
    if change == "sub":
        assert len(pc["new"]) <= 3 and pc["new"] != pc["old"]


def test_protein_change_checks_the_wild_type_residue():
    from plasmidbank.tools.mutagenesis import protein_change
    ref = _cds_map(1)
    with pytest.raises(ValueError, match="codon 1 of gene is M"):
        protein_change(ref, ref.features[0], "K1R")
    with pytest.raises(ValueError):
        protein_change(ref, ref.features[0], "M500R")
    with pytest.raises(ValueError):
        protein_change(ref, ref.features[0], "hello")


@pytest.mark.parametrize("start,end,bases,want_seq,want_feats", [
    (5, 4, "GG", "AAAAGGCCCCCCCCCCTTTTT", [(7, 16)]),          # insert before the feature: it shifts
    (8, 9, "", "AAAACCCCCCCCTTTTT", [(5, 12)]),               # delete inside: it shrinks
    (5, 14, "", "AAAATTTTT", []),                              # the whole feature deleted: dropped
    (4, 15, "G", "AAAGTTTT", []),                              # inside a longer replaced stretch: dropped
    (3, 6, "TT", "AATTCCCCCCCCTTTTT", [(5, 12)]),             # overlaps its start: trimmed, follows the new bases
    (5, 14, "GGG", "AAAAGGGTTTTT", [(5, 7)]),                  # exactly its bases replaced: kept on the new ones
])
def test_sequence_edits_move_features(tmp_path, start, end, bases, want_seq, want_feats):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord
    from plasmidbank.tools.edit import apply_edits

    rec = SeqRecord(Seq("AAAA" + "C" * 10 + "TTTTT"), id="x", name="x",
                    annotations={"molecule_type": "DNA", "topology": "linear"},
                    features=[SeqFeature(SimpleLocation(4, 14, 1), type="exon", qualifiers={"label": ["e1"]})])
    f = tmp_path / "x.gb"
    SeqIO.write(rec, f, "genbank")
    text, _ = apply_edits(str(f), {"seq": [{"start": start, "end": end, "bases": bases}]})
    f.write_text(text)
    out = SeqIO.read(f, "genbank")
    assert str(out.seq) == want_seq
    assert [(int(x.location.start) + 1, int(x.location.end)) for x in out.features if x.type == "exon"] == want_feats


def test_sequence_edit_refuses_bad_input(tmp_path):
    from plasmidbank.tools.edit import apply_edits

    m = _random_map(n=200)
    f = tmp_path / "x.gb"
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    SeqIO.write(SeqRecord(Seq(m.seq), id="x", annotations={"molecule_type": "DNA"}), f, "genbank")
    for bad in ([{"start": 1, "end": 1, "bases": "XYZ"}], [{"start": 199, "end": 205, "bases": ""}],
                [{"start": 10, "end": 20, "bases": "A"}, {"start": 15, "end": 25, "bases": "C"}]):
        with pytest.raises(ValueError):
            apply_edits(str(f), {"seq": bad})


def test_bank_and_scan_skip_workspaces(tmp_path):
    (tmp_path / "Workspaces" / "ann").mkdir(parents=True)
    (tmp_path / "Workspaces" / "ann" / "PMC-00001_mine.gb").write_text(">x\nACGT\n")
    (tmp_path / "PMC-00002.gb").write_text(">y\nACGT\n")
    assert [m.plasmid_id for m in Bank.load(tmp_path).maps] == ["PMC-00002"]


def test_sites_carry_recognition_and_overhang():
    from plasmidbank.tools.sites import find_sites

    m = _random_map(n=600, seed=9)
    m.seq = m.seq[:100] + "GAATTC" + m.seq[106:300] + "GATATC" + m.seq[306:400] + "CTGCAG" + m.seq[406:]
    got = {s["enzyme"]: (s["recognition"], s["overhang"]) for s in find_sites(m, "common", "EcoRI, EcoRV, PstI")}
    assert got == {"EcoRI": ("GAATTC", "5′ AATT"), "EcoRV": ("GATATC", "blunt"), "PstI": ("CTGCAG", "3′ TGCA")}
    # the sequence view draws the recognition sequence and the cut on both strands:
    # top cut between bases c-1 and c (0-based boundary c-1), bottom cut at boundary c-1-ovhg
    by = {s["enzyme"]: s for s in find_sites(m, "common", "EcoRI, PstI")}
    eco, pst = by["EcoRI"], by["PstI"]
    assert 100 in eco["spans"] and eco["size"] == 6 and 102 in eco["cuts"]  # G^AATT_C
    assert (102 - 1, 102 - 1 - eco["ovhg"]) == (101, 105)
    assert 400 in pst["spans"] and 406 in pst["cuts"] and (406 - 1, 406 - 1 - pst["ovhg"]) == (405, 401)  # C_TGCA^G
    r = _random_map(n=300, seed=4)
    r.seq = r.seq[:50] + "GGTCTC" + "A" * 20 + "T" * 20 + "GAGACC" + r.seq[96:]  # BsaI both ways round
    bsa = {s["enzyme"]: s for s in find_sites(r, "goldengate", "BsaI")}["BsaI"]
    assert {50, 96} <= set(bsa["spans"])
    tops = {c - 1 for c in bsa["cuts"]}
    assert {57, 96 - 5} <= tops and all(t - bsa["ovhg"] == t + 4 for t in tops)  # 5' NNNN overhang on both


def test_reverse_complement_edit_flips_features(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord
    from plasmidbank.tools.edit import apply_edits

    rec = SeqRecord(Seq("AAAACCCGGGTT"), id="x", name="x", annotations={"molecule_type": "DNA", "topology": "circular"},
                    features=[SeqFeature(SimpleLocation(0, 4, 1), type="CDS", qualifiers={"label": ["a"]})])
    f = tmp_path / "x.gb"
    SeqIO.write(rec, f, "genbank")
    f.write_text(apply_edits(str(f), {"reverse_complement": True})[0])
    out = SeqIO.read(f, "genbank")
    cds = next(x for x in out.features if x.type == "CDS")
    assert str(out.seq) == "AACCCGGGTTTT" and (int(cds.location.start), int(cds.location.end), cds.location.strand) == (8, 12, -1)


@pytest.mark.parametrize("ftype,ok", [("exon", True), ("sgRNA_target", True), ("Ann's_tag", True),
                                      ("bad type", False), ("source", False), ("9lives", False), ("x" * 21, False)])
def test_own_feature_types(tmp_path, ftype, ok):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from plasmidbank.tools.edit import apply_edits

    f = tmp_path / "x.gb"
    SeqIO.write(SeqRecord(Seq("ACGT" * 20), id="x", annotations={"molecule_type": "DNA"}), f, "genbank")
    edits = {"add": [{"type": ftype, "label": "t", "start": 1, "end": 10}]}
    if not ok:
        with pytest.raises(ValueError):
            apply_edits(str(f), edits)
        return
    f.write_text(apply_edits(str(f), edits)[0])
    assert [x.type for x in SeqIO.read(f, "genbank").features] == [ftype]


@pytest.mark.parametrize("start,end", [(1, 1), (40, 40), (79, 80)])
def test_point_feature_of_one_base_can_be_added(tmp_path, start, end):
    # a point mutation (A-->G) is annotated as a 1-bp feature from "save selection"
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from plasmidbank.tools.edit import apply_edits

    f = tmp_path / "x.gb"
    SeqIO.write(SeqRecord(Seq("ACGT" * 20), id="x", annotations={"molecule_type": "DNA"}), f, "genbank")
    f.write_text(apply_edits(str(f), {"add": [{"type": "variation", "label": "A-->G", "start": start, "end": end}]})[0])
    (feat,) = SeqIO.read(f, "genbank").features
    assert (int(feat.location.start) + 1, int(feat.location.end), feat.type) == (start, end, "variation")


@pytest.mark.parametrize("minus", [False, True])
def test_exons_join_into_one_cds_on_either_strand(minus):
    from plasmidbank.tools.analysis import join_features

    cds = "ATGGCT" + "GGTAAA" + "TAA"  # M A G K *, split after 6 nt by an intron
    seq = "CC" + cds[:6] + "GTAAGTTTTAG" + cds[6:] + "CC"
    L, ex = len(seq), [(2, 8), (19, 28)]
    if minus:
        seq, ex = revcomp(seq), [(L - b, L - a) for a, b in ex]
    m = _random_map()
    m.seq, m.topology = seq, "linear"
    m.features = [Feature("exon", f"e{i}", a, b, -1 if minus else 1) for i, (a, b) in enumerate(ex)]
    r = join_features(m, [1, 0])  # the order picked does not matter
    assert (r["protein"], r["multiple_of_3"], r["starts_with_atg"], r["ends_with_stop"]) == ("MAGK*", True, True, True)
    m.features[1] = Feature("exon", "e1", m.features[1].start, m.features[1].end - 1, m.features[1].strand)
    assert not join_features(m, [0, 1])["multiple_of_3"]  # one base short: the frame breaks


def test_joined_cds_is_saved_as_one_feature(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from plasmidbank.tools.edit import apply_edits

    f = tmp_path / "x.gb"
    SeqIO.write(SeqRecord(Seq("ACGT" * 20), id="x", annotations={"molecule_type": "DNA"}), f, "genbank")
    f.write_text(apply_edits(str(f), {"add": [{"type": "CDS", "label": "g", "strand": -1, "parts": [[40, 50], [5, 20]]}]})[0])
    loc = SeqIO.read(f, "genbank").features[0].location
    assert [(int(x.start) + 1, int(x.end)) for x in loc.parts] == [(40, 50), (5, 20)] and loc.strand == -1


@pytest.mark.parametrize("make", [
    lambda s: revcomp(s[2900:] + s[:400]),                # a read on the minus strand, across the origin
    lambda s: s[1000:1700],                               # a plain read
    lambda s: revcomp(s[1200:] + s[:1200]),               # a whole-plasmid consensus, rotated, minus strand
])
def test_reads_land_on_the_map_where_their_blocks_say(make):
    # the stacked view draws each read from `query` + `blocks`: they must agree with the map
    ref = _random_map(n=3000, seed=12)
    r = verify(ref, Query("r", make(ref.seq)))
    L, q = ref.length, r.query
    assert r.blocks and all(ref.seq[(mp - 1 + k) % L] == q[(qp - 1 + k) % len(q)]
                            for mp, qp, n in r.blocks for k in range(n))


def _pcr_map():
    ref = _random_map(n=3000, seed=21)
    fwd = "GGTCTCA" + ref.seq[100:122]  # a 7-nt tail, 22 nt binding
    rev = revcomp(ref.seq[600:622])
    return ref, [{"name": "F1", "seq": fwd}, {"name": "R1", "seq": rev}]


def test_pcr_product_is_the_region_plus_tails():
    from plasmidbank.tools.pcr import pcr

    ref, prs = _pcr_map()
    r = pcr(ref, prs[::-1], "Q5")  # the order picked does not matter
    assert r["product"] == 7 + 522 and r["product_seq"] == "GGTCTCA" + ref.seq[100:622]
    assert (r["forward"]["name"], r["reverse"]["name"], r["map_end"]) == ("F1", "R1", 622)
    assert r["map_start"] <= 101  # the matcher may count mismatched 5' bases as binding, as the Primers tab does


def test_pcr_across_the_origin_of_a_circular_map():
    from plasmidbank.tools.pcr import pcr

    ref = _random_map(n=3000, seed=22)
    r = pcr(ref, [{"name": "F", "seq": ref.seq[2800:2822]}, {"name": "R", "seq": revcomp(ref.seq[200:222])}])
    assert r["product_seq"] == ref.seq[2800:] + ref.seq[:222]


@pytest.mark.parametrize("pol,check", [
    ("Q5", lambda r, low: r["annealing"] == min(round(low + 3, 1), 72)),
    ("Phusion", lambda r, low: r["annealing"] == min(round(low + 3, 1), 72)),  # 22-nt binding parts: > 20 nt
    ("SuperFi II", lambda r, low: r["annealing"] == 60 and not r["two_step"]),
    ("OneTaq", lambda r, low: r["annealing"] == min(round(low - 5, 1), 68)),
])
def test_pcr_annealing_follows_the_vendor_rule(pol, check):
    from plasmidbank.tools.pcr import pcr

    ref, prs = _pcr_map()
    r = pcr(ref, prs, pol)
    low = min(r["forward"]["tm"], r["reverse"]["tm"])
    assert check(r, low) and r["sources"] and r["program"][-1]["step"] == "final extension"


def test_pcr_refuses_primers_that_do_not_face_each_other():
    from plasmidbank.tools.pcr import pcr

    ref, prs = _pcr_map()
    ref.topology = "linear"
    with pytest.raises(ValueError):
        pcr(ref, [prs[0], {"name": "F2", "seq": ref.seq[300:322]}])  # both forward
    with pytest.raises(ValueError):
        pcr(ref, [{"name": "R", "seq": revcomp(ref.seq[50:72])}, {"name": "F", "seq": ref.seq[900:922]}])  # facing away


def test_joined_feature_on_a_linear_record_never_wraps(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    loc = CompoundLocation([SimpleLocation(4, 19, 1), SimpleLocation(35, 50, 1)])
    for topo, want in (("linear", (4, 50)), ("circular", (35, 19))):  # a circle keeps the widest-gap rule
        f = tmp_path / f"{topo}.gb"
        SeqIO.write(SeqRecord(Seq("A" * 54), id="g", annotations={"molecule_type": "DNA", "topology": topo},
                              features=[SeqFeature(loc, type="CDS")]), f, "genbank")
        x = read_map(f).features[0]
        assert (x.start, x.end) == want, topo


def test_ab1_reads_are_named_after_their_file():
    # the stacked reads showed the trace's sample number, not the name she gave the tube
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from plasmidbank.tools import align

    fake = SeqRecord(Seq("ACGT" * 30), id="31234567890", letter_annotations={"phred_quality": [40] * 120})
    fake.annotations["abif_raw"] = {}
    orig = align.SeqIO.read
    align.SeqIO.read = lambda *a, **k: fake
    try:
        q = align.load_queries("reads/ILC178_M13F.ab1", b"ABIF....")
    finally:
        align.SeqIO.read = orig
    assert [x.name for x in q] == ["ILC178_M13F"]
