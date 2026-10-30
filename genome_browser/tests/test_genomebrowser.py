"""Build + browse + import, on small synthetic genomes (no real data needed).

Requests go through the real ASGI app (middleware included) with a tiny in-process client.
"""

import asyncio
import gzip
import json
import random
import sqlite3
from pathlib import Path

import pytest

from genomebrowser import build, imports, server, store

random.seed(1)
CHR1 = "".join(random.choice("ACGT") for _ in range(250_000))
CHR2 = "".join(random.choice("ACGT") for _ in range(1_234))

GENES = """##gff-version 3
chr1\tsrc\tgene\t1001\t3000\t.\t-\t.\tID=G1;Name=KIN1;description=TEST kinase
chr1\tsrc\tmRNA\t1001\t3000\t.\t-\t.\tID=G1.1;Parent=G1
chr1\tsrc\texon\t1001\t1500\t.\t-\t.\tID=G1.1.exon.1;Parent=G1.1
chr1\tsrc\texon\t2501\t3000\t.\t-\t.\tID=G1.1.exon.2;Parent=G1.1
chr1\tsrc\tCDS\t1101\t1500\t.\t-\t0\tID=G1.1.cds.1;Parent=G1.1
chr1\tsrc\tCDS\t2501\t2900\t.\t-\t0\tID=G1.1.cds.2;Parent=G1.1
chr1\tsrc\tgene\t200001\t201000\t.\t+\t.\tID=G2_LC
chr1\tsrc\tCDS\t200001\t201000\t.\t+\t0\tID=cds-G2;Parent=G2_LC;product=photosystem%20II%20protein
chr9\tsrc\tgene\t1\t10\t.\t+\t.\tID=NOT_IN_GENOME
"""
INSERTIONS = """chr1\tlab\tinsertion\t1200\t1200\t.\t+\t.\tID=07-030055169
chr1\tlab\tinsertion\t1200\t1200\t.\t+\t.\tID=07-030055169
chr1\tlab\tinsertion\t2000\t2000\t.\t-\t.\tID=07-030055169
chr2\tlab\tinsertion\t50\t50\t.\t+\t.\tID=07-030055169
chr1\tlab\tinsertion\t1050\t1050\t.\t+\t.\tID=30-099999999
"""
BED = "track name=peaks\nchr1\t1199\t1300\tpeak1\t0\t+\nchr1\t5000\t5100\t.\t0\t.\n"


def ragged_fasta(path: Path):
    """Uneven line lengths, CRLF, a description, gzip: everything normalise_fasta must absorb."""
    text = ">chr1 the first one\r\n"
    i = 0
    while i < len(CHR1):
        n = random.choice((50, 61, 80))
        text += CHR1[i : i + n] + "\r\n"
        i += n
    text += ">chr2\n" + CHR2 + "\n"
    path.write_bytes(gzip.compress(text.encode()))


@pytest.fixture
def gdir(tmp_path, monkeypatch):
    d = tmp_path / "genomes"
    monkeypatch.setenv("GENOMES_DIR", str(d))
    monkeypatch.delenv("GENOMES_IMPORT_DIR", raising=False)
    monkeypatch.setattr(server, "AUTH", "none")
    return d


@pytest.fixture
def built(gdir, tmp_path):
    fa = tmp_path / "g.fa.gz"
    ragged_fasta(fa)
    for name, text in (
        ("genes.gff3", GENES),
        ("ins.gff3", INSERTIONS),
        ("peaks.bed", BED),
    ):
        (tmp_path / name).write_text(text)
    return build.new_genome(
        gdir,
        "Test plant v1",
        fa,
        [
            {"name": "genes", "path": str(tmp_path / "genes.gff3"), "kind": "genes"},
            {
                "name": "insertions",
                "path": str(tmp_path / "ins.gff3"),
                "kind": "features",
            },
            {"name": "peaks", "path": str(tmp_path / "peaks.bed"), "kind": "features"},
        ],
        user="anna",
    )


# ── a minimal ASGI client ────────────────────────────────────────────────────


def call(method, path, body=b"", headers=None, query=""):
    headers = dict(headers or {})
    if isinstance(body, (dict, list)):
        body = json.dumps(body).encode()
        headers["content-type"] = "application/json"
    scope = {
        "type": "http",
        "method": method,
        "path": path,
        "raw_path": path.encode(),
        "query_string": query.encode(),
        "root_path": "",
        "scheme": "http",
        "server": ("test", 80),
        "client": ("127.0.0.1", 1),
        "http_version": "1.1",
        "headers": [(b"host", b"test")]
        + [(k.lower().encode(), v.encode()) for k, v in headers.items()],
    }
    sent, out = [False], {"body": b""}

    async def receive():
        if sent[0]:
            return {"type": "http.disconnect"}
        sent[0] = True
        return {"type": "http.request", "body": body, "more_body": False}

    async def send(msg):
        if msg["type"] == "http.response.start":
            out["status"] = msg["status"]
        elif msg["type"] == "http.response.body":
            out["body"] += msg.get("body", b"")

    asyncio.run(server.app(scope, receive, send))
    try:
        out["json"] = json.loads(out["body"])
    except ValueError:
        out["json"] = None
    return out


# ── build ────────────────────────────────────────────────────────────────────


def test_fasta_is_normalised_and_bases_come_back_exactly(built, gdir):
    fa = (gdir / built / "genome.fa").read_text()
    assert fa.startswith(">chr1\n") and "\r" not in fa
    assert {len(x) for x in fa.splitlines()[1:5]} == {60}
    c = store.connect(built)
    assert server.fetch_bases(built, c, "chr1", 1, 70) == CHR1[:70]
    assert server.fetch_bases(built, c, "chr1", 249_990, 250_000) == CHR1[-11:]
    assert (
        server.fetch_bases(built, c, "chr2", 1, 2000) == CHR2
    )  # clipped to the length


def test_tracks_bed_and_unknown_sequences(built):
    c = store.connect(built)
    tr = {
        r["name"]: (r["kind"], r["format"], r["n"])
        for r in c.execute("SELECT * FROM tracks")
    }
    assert tr == {
        "genes": ("genes", "gff3", 8),
        "insertions": ("features", "gff3", 4),
        "peaks": ("features", "bed", 2),
    }
    peak = c.execute(
        "SELECT start, \"end\", strand FROM feat WHERE fid = 'peak1'"
    ).fetchone()
    assert tuple(peak) == (1200, 1300, "+")  # BED 0-based half-open → 1-based closed


def test_bad_inputs_say_what_is_wrong(gdir, tmp_path):
    bad = tmp_path / "x.fa"
    bad.write_text("ACGT\n")
    with pytest.raises(build.BuildError, match="does not start with"):
        build.new_genome(gdir, "Bad one", bad, [])
    fa = tmp_path / "ok.fa"
    fa.write_text(">contig_1\nACGTACGT\n")
    g = tmp_path / "g.gff3"
    g.write_text("Chr1\ts\tgene\t1\t5\t.\t+\t.\tID=a\n")
    with pytest.raises(build.BuildError, match="'Chr1'.*'contig_1'"):
        build.new_genome(
            gdir, "Mismatch", fa, [{"name": "genes", "path": str(g), "kind": "genes"}]
        )
    (tmp_path / "dup.fa").write_text(">a\nAC\n>a\nGT\n")
    with pytest.raises(build.BuildError, match="appears twice"):
        build.new_genome(gdir, "Dup", tmp_path / "dup.fa", [])
    assert not any(gdir.glob("*/genome.sqlite"))  # nothing half-built is left behind


def test_slugs_never_escape(gdir):
    assert build.slugify("Barley Morex V3 (2021)") == "barley-morex-v3-2021"
    assert store.genome_path("../etc") is None and store.genome_path("a/b") is None
    with pytest.raises(build.BuildError):
        build.slugify("!")


# ── reading ──────────────────────────────────────────────────────────────────


def test_region_detail_and_density(built):
    r = call("GET", f"/api/g/{built}/region", query="seq=chr1&start=900&end=3100")[
        "json"
    ]
    gene = [f for f in r["tracks"]["genes"] if f["type"] == "gene"]
    assert [g["fid"] for g in gene] == ["G1"] and gene[0]["defline"] == "TEST kinase"
    assert sorted(f["start"] for f in r["tracks"]["insertions"]) == [1050, 1200, 2000]
    assert [f["start"] for f in r["tracks"]["peaks"]] == [1200]
    assert "bases" not in r
    assert (
        call("GET", f"/api/g/{built}/region", query="seq=chr1&start=1&end=70")["json"][
            "bases"
        ]
        == CHR1[:70]
    )


def test_search_by_id_name_and_description(built):
    s = lambda q: [
        h["fid"] for h in call("GET", f"/api/g/{built}/search", query=f"q={q}")["json"]
    ]  # noqa: E731
    assert s("G1") == ["G1", "G1.1"]
    assert s("KIN1") == ["G1"]
    assert s("kinase") == ["G1"]
    assert s("photosystem") == ["G2_LC"]  # product= on a CDS describes its gene
    assert s("30055169") == ["07-030055169"]  # an ID without its prefix


def test_gene_card_lists_the_copies_of_a_gene(built, gdir):
    # RinRK1 has two copies; the card lists every gene of the PANTHER subfamily
    db = sqlite3.connect(gdir / built / "genome.sqlite")
    with db:
        db.execute("UPDATE gene_info SET defline = '(1 of 2) PTHR48010:SF104 - LRR KINASE' WHERE gene = 'G1'")
        db.execute("INSERT OR REPLACE INTO gene_info (gene, defline) VALUES ('G2_LC', '(1 of 2) PTHR48010:SF104 - LRR KINASE')")
    db.close()
    fam = call("GET", f"/api/g/{built}/feature/G1")["json"]["family"]
    assert fam["id"] == "PTHR48010:SF104"
    assert [(g["fid"], g["seqid"], g["start"]) for g in fam["genes"]] == [("G1", "chr1", 1001), ("G2_LC", "chr1", 200001)]
    with sqlite3.connect(gdir / built / "genome.sqlite") as db:  # a longer subfamily number is another family
        db.execute("UPDATE gene_info SET defline = 'PTHR48010:SF1040' WHERE gene = 'G2_LC'")
    assert [g["fid"] for g in call("GET", f"/api/g/{built}/feature/G1")["json"]["family"]["genes"]] == ["G1"]
    assert call("GET", f"/api/g/{built}/feature/G2_LC")["json"]["family"]["genes"][0]["fid"] == "G2_LC"


def test_gene_card_classifies_overlapping_features(built):
    g = call("GET", f"/api/g/{built}/feature/G1")["json"]
    assert g["defline"] == "TEST kinase" and [t["fid"] for t in g["transcripts"]] == [
        "G1.1"
    ]
    where = {(o["track"], o["start"]): o["where"] for o in g["overlaps"]}
    assert where == {
        ("insertions", 1050): "exon (UTR)",
        ("insertions", 1200): "CDS",
        ("insertions", 2000): "intron",
        ("peaks", 1200): "CDS",
    }


def test_transcript_is_spliced_on_the_minus_strand(built):
    t = call("GET", f"/api/g/{built}/transcript/G1.1")["json"]
    rc = lambda a, b: server.revcomp(CHR1[a - 1 : b])  # noqa: E731
    assert t["mrna"] == rc(2501, 3000) + rc(1001, 1500)
    assert t["cds"] == rc(2501, 2900) + rc(1101, 1500)
    assert t["protein"] == server.translate(t["cds"]).removesuffix("*")
    assert (t["utr5"], t["utr3"]) == (100, 100)
    assert [(e["n"], e["start"], e["utr5"], e["utr3"], e["cds"], e["phase"]) for e in t["exons"]] == [
        (1, 2501, 100, 0, 400, 0),
        (2, 1001, 0, 100, 400, 2),
    ]
    (i,) = t["introns"]
    s = rc(1501, 2500)
    assert (i["start"], i["end"], i["splice"]) == (1501, 2500, f"{s[:2]}…{s[-2:]}")
    g = call("GET", f"/api/g/{built}/transcript/G2_LC")["json"]  # gene → CDS, no mRNA
    assert g["mrna"] == g["cds"] == CHR1[200000:201000] and g["introns"] == []
    assert call("GET", f"/api/g/{built}/transcript/nope")["status"] == 404


def test_translate():
    assert server.translate("ATGTGGTAANNNGGNTTRTA") == "MW*XGL"


def test_occurrences_of_one_id(built):
    d = call("GET", f"/api/g/{built}/occurrences/insertions/07-030055169")["json"]
    assert [
        (o["seqid"], o["start"], o["gene"], o["where"]) for o in d["occurrences"]
    ] == [
        ("chr1", 1200, "G1", "CDS"),
        ("chr1", 2000, "G1", "intron"),
        ("chr2", 50, None, None),
    ]
    assert call("GET", f"/api/g/{built}/occurrences/insertions/nope")["status"] == 404


def test_chromosomes_longer_than_2_31(built, gdir, tmp_path):
    """Faba bean chromosomes pass 2 Gb: the float R*Tree plus the exact filter must still be right."""
    db = sqlite3.connect(gdir / built / "genome.sqlite")
    db.execute(
        "INSERT INTO seq(name, length, offset, linebases, linewidth) VALUES ('chrBig', 3200000000, 0, 60, 61)"
    )
    db.commit()
    db.close()
    big = tmp_path / "big.bed"
    big.write_text(
        "chrBig\t2500000000\t2500000001\tfarA\nchrBig\t2500000300\t2500000301\tfarB\n"
    )
    build.add_tracks(
        gdir / built, [{"name": "far", "path": str(big), "kind": "features"}]
    )

    def q(a, b):
        r = call(
            "GET", f"/api/g/{built}/region", query=f"seq=chrBig&start={a}&end={b}"
        )["json"]
        return [f["fid"] for f in r["tracks"]["far"]]

    assert q(2_500_000_001 - 10, 2_500_000_001 + 10) == ["farA"]
    assert (
        q(2_500_000_002, 2_500_000_299) == []
    )  # within float32 rounding of both, but empty
    assert q(2_500_000_301, 2_500_000_400) == ["farB"]
    dens = call(
        "GET", f"/api/g/{built}/region", query="seq=chrBig&start=1&end=3200000000"
    )["json"]
    assert (
        dens["bin"] > build.BIN and len(dens["density"]["far"]) == 1
    )  # re-binned to ≤ MAX_BINS


def test_fasta_endpoint(built):
    r = call(
        "GET", f"/api/g/{built}/fasta", query="seq=chr1&start=1&end=10&strand=-&fid=x"
    )
    assert r["body"].decode().splitlines() == [
        ">x chr1:1-10(-)",
        server.revcomp(CHR1[:10]),
    ]


# ── importing ────────────────────────────────────────────────────────────────


def test_chunked_upload_resumes_and_refuses_gaps(gdir):
    data = b">a\n" + b"ACGT" * 100 + b"\n"
    up = call("POST", "/api/uploads", {"filename": "a.fa", "size": len(data)})["json"]

    def put(off, chunk):
        return call("PUT", f"/api/uploads/{up['id']}", chunk, query=f"offset={off}")

    assert put(0, data[:100])["json"]["received"] == 100
    assert (
        put(0, data[:100])["status"] == 409
    )  # already have those bytes: resume from 100
    assert put(150, data[150:])["status"] == 409  # a gap
    assert call("GET", f"/api/uploads/{up['id']}")["json"]["received"] == 100
    assert put(100, data[100:])["json"]["received"] == len(data)
    p = imports.upload_file(up["id"], "local")
    assert p.read_bytes() == data and p.name.endswith("__a.fa")
    assert (
        call("POST", "/api/uploads", {"filename": "notes.docx", "size": 10})["status"]
        == 400
    )


def test_upload_refused_when_the_disk_is_full(gdir, monkeypatch):
    import shutil

    monkeypatch.setattr(
        imports.shutil, "disk_usage", lambda p: shutil._ntuple_diskusage(100, 99, 10**6)
    )
    r = call("POST", "/api/uploads", {"filename": "big.fa.gz", "size": 10**6})
    assert r["status"] == 400 and "not enough disk space" in r["json"]["detail"]


def test_upload_to_genome_through_the_queue(gdir):
    fa = b">c1\n" + CHR2.encode() + b"\n"
    gff = b"c1\tx\tgene\t10\t200\t.\t+\t.\tID=g1;Name=abc\n"
    ids = []
    for name, data in (("c.fa", fa), ("c.gff3", gff)):
        u = call("POST", "/api/uploads", {"filename": name, "size": len(data)})["json"]
        call("PUT", f"/api/uploads/{u['id']}", data, query="offset=0")
        ids.append(u["id"])
    job = call(
        "POST",
        "/api/genomes",
        {
            "name": "Tiny",
            "fasta": {"upload": ids[0]},
            "tracks": [
                {"name": "genes", "kind": "genes", "source": {"upload": ids[1]}}
            ],
        },
    )["json"]
    imports._run(
        imports.get_job(job["id"])
    )  # what the worker thread does, synchronously
    j = imports.get_job(job["id"])
    assert j["status"] == "done", j["log"]
    assert [g["slug"] for g in call("GET", "/api/genomes")["json"]] == ["tiny"]
    assert not list((gdir / ".incoming").iterdir())  # inputs are consumed
    again = call("POST", "/api/genomes", {"name": "Tiny", "fasta": {"upload": ids[0]}})
    assert again["status"] == 400 and "exists already" in again["json"]["detail"]


def test_server_import_folder_is_admin_only_and_contained(gdir, tmp_path, monkeypatch):
    root = tmp_path / "share"
    (root / "barley").mkdir(parents=True)
    (root / "barley" / "morex.fa.gz").write_bytes(b"x")
    (tmp_path / "secret.fa").write_text(">s\nA\n")
    monkeypatch.setenv("GENOMES_IMPORT_DIR", str(root))
    assert [f["path"] for f in imports.server_files()] == ["barley/morex.fa.gz"]
    with pytest.raises(imports.ImportError_):
        imports.server_file("../secret.fa")
    monkeypatch.setattr(server, "AUTH", "header")
    h = {"x-auth-user": "ben", "x-auth-role": "user"}
    assert call("GET", "/api/server-files", headers=h)["status"] == 403
    r = call(
        "POST",
        "/api/genomes",
        {"name": "Barley", "fasta": {"path": "barley/morex.fa.gz"}},
        headers=h,
    )
    assert r["status"] == 403  # adding is for admins (GENOMES_UPLOADS=admins, the default)


# ── who may do what ──────────────────────────────────────────────────────────


def test_header_mode_needs_the_proxy_headers(built, monkeypatch):
    monkeypatch.setattr(server, "AUTH", "header")
    assert call("GET", "/healthz")["status"] == 200
    assert call("GET", "/api/genomes")["status"] == 401
    me = call("GET", "/api/me", headers={"x-auth-user": "anna", "x-auth-role": "user"})[
        "json"
    ]
    assert me["user"] == "anna" and me["admin"] is False
    assert (
        call("DELETE", f"/api/genomes/{built}", headers={"x-auth-user": "ben"})[
            "status"
        ]
        == 403
    )
    cross = {"x-auth-user": "anna", "origin": "https://evil.example.org"}
    assert call("DELETE", f"/api/genomes/{built}", headers=cross)["status"] == 403
    assert (
        call("DELETE", f"/api/genomes/{built}", headers={"x-auth-user": "anna"})[
            "status"
        ]
        == 200
    )
    assert store.list_genomes() == []


def test_drop_track(built, gdir):
    build.drop_track(gdir / built, "peaks")
    r = call("GET", f"/api/g/{built}/region", query="seq=chr1&start=1&end=3000")["json"]
    assert "peaks" not in r["tracks"] and r["tracks"]["insertions"]
    assert call("GET", f"/api/g/{built}/search", query="q=peak1")["json"] == []


def test_seed_runs_once(gdir, tmp_path):
    seed = tmp_path / "seed"
    seed.mkdir()
    (seed / "g.fa").write_text(">s1\nACGTACGTAC\n")
    (seed / "seed.json").write_text(json.dumps({"name": "Demo", "fasta": "g.fa"}))
    job = imports.seed(seed)
    assert job and imports.seed(seed) is None  # queued already: not twice
    imports._run(imports.get_job(job["id"]))
    assert store.list_genomes()[0]["name"] == "Demo"


def test_only_admins_add_unless_everyone_may(built, monkeypatch):
    monkeypatch.setattr(server, "AUTH", "header")
    user = {"x-auth-user": "ben", "x-auth-role": "user"}
    assert call("GET", "/api/me", headers=user)["json"]["can_add"] is False
    assert call("POST", "/api/uploads", {"filename": "a.fa", "size": 10}, headers=user)["status"] == 403
    assert call("POST", f"/api/g/{built}/tracks", {"tracks": []}, headers=user)["status"] == 403
    admin = {"x-auth-user": "boss", "x-auth-role": "admin"}
    assert call("POST", "/api/uploads", {"filename": "a.fa", "size": 10}, headers=admin)["status"] == 200
    monkeypatch.setattr(server, "UPLOADS", "everyone")
    assert call("POST", "/api/uploads", {"filename": "a.fa", "size": 10}, headers=user)["status"] == 200
