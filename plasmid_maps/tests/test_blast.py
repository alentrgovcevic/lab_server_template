"""BLAST against the bank. Skipped when BLAST+ is not installed."""

import random

import pytest

from plasmidbank.formats import MapRecord
from plasmidbank.sequtil import revcomp
from plasmidbank.tools import blast

pytestmark = pytest.mark.skipif(not blast.available(), reason="BLAST+ not installed")


def _map(n, seed, pid, circular=True):
    rnd = random.Random(seed)
    return MapRecord(f"/bank/{pid}.gb", pid, "PMC", "filename", [], "genbank", "ok",
                     seq="".join(rnd.choice("ACGT") for _ in range(n)), topology="circular" if circular else "linear")


@pytest.fixture
def bank(tmp_path):
    maps = [_map(3000, i, f"PMC-0010{i}") for i in range(4)]
    maps.append(_map(500, 99, "PMC-00199", circular=False))
    return maps, tmp_path / "blast"


def test_a_region_finds_its_map_on_both_strands(bank):
    maps, folder = bank
    q = maps[2].seq[1000:1400]
    r = blast.search(q, maps, folder)
    top = r["hits"][0]
    assert top["plasmid_id"] == "PMC-00102" and top["identity"] == 100 and top["coverage"] == 100
    assert (top["map_from"], top["map_to"], top["strand"]) == (1001, 1400, "+")
    top = blast.search(revcomp(q), maps, folder)["hits"][0]
    assert top["plasmid_id"] == "PMC-00102" and top["strand"] == "-" and (top["map_from"], top["map_to"]) == (1001, 1400)


def test_a_query_across_the_origin_is_found_and_folded_back(bank):
    maps, folder = bank
    s = maps[1].seq
    top = blast.search(s[-150:] + s[:150], maps, folder)["hits"][0]
    assert top["plasmid_id"] == "PMC-00101" and top["coverage"] == 100
    assert (top["map_from"], top["map_to"]) == (2851, 150)  # folded: wraps the origin


def test_short_query_and_rebuild_on_change(bank):
    maps, folder = bank
    r = blast.search(maps[0].seq[200:225], maps, folder)  # a primer-sized query
    assert r["program"] == "blastn-short" and r["hits"][0]["plasmid_id"] == "PMC-00100"
    meta = blast.ensure_db(maps, folder)
    maps.append(_map(2000, 50, "PMC-00150"))
    assert blast.ensure_db(maps, folder)["fingerprint"] != meta["fingerprint"]
    assert blast.search(maps[-1].seq[100:600], maps, folder)["hits"][0]["plasmid_id"] == "PMC-00150"
    with pytest.raises(ValueError):
        blast.search("ACGT", maps, folder)
