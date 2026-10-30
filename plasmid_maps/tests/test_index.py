import json
import os
import shutil
from pathlib import Path

import pytest

from plasmidbank.index import Index, bin_of, parse_bin
from plasmidbank.scan import render_report, scan

SAMPLE = Path(__file__).resolve().parents[2] / "examples" / "bank"
CTRL = "PMC-03207 pBR322 CTRL.gbk"
CTRL1 = "PMC-03207 pBR322 CTRL.1.gbk"


@pytest.fixture
def bank(tmp_path):
    root = tmp_path / "bank"
    shutil.copytree(SAMPLE, root)
    return root, tmp_path / "work" / "index.sqlite", tmp_path / "work" / "mirror"


def run(bank):
    root, db, mirror = bank
    return scan(root, db, mirror, log=lambda *_: None)


def names(rows):
    return sorted(Path(r["path"]).name for r in rows)


def test_rescan_reuses_cache(bank):
    first = run(bank)
    n = first["seen"]
    assert n >= 10 and first["parsed"] == n and first["cached"] == 0
    second = run(bank)
    assert (second["seen"], second["parsed"], second["cached"]) == (n, 0, n)
    assert second["mirrored_new"] == 0

    root, db, _ = bank
    p = root / "pmc_gb" / "PMC-03204.gbk"
    p.write_bytes(p.read_bytes() + b"\n")
    third = run(bank)
    assert (third["parsed"], third["cached"]) == (1, n - 1)


def test_missing_since_set_and_row_kept(bank):
    root, db, _ = bank
    run(bank)
    gone = root / "pmc_gb" / "PMC-03204.gbk"
    gone.unlink()
    counts = run(bank)
    assert counts["newly_missing"] == 1
    idx = Index(db)
    row = idx.get(str(gone.resolve()))
    assert row is not None and row["missing_since"] is not None
    assert idx.by_id("PMC-03204") == []
    assert len(idx.by_id("PMC-03204", include_missing=True)) == 1
    assert names(idx.missing()) == ["PMC-03204.gbk"]
    idx.close()


def test_missing_marked_for_back_to_back_runs(bank, monkeypatch):
    import plasmidbank.index as ix

    monkeypatch.setattr(
        ix, "now", lambda: "2026-10-05T09:30:00"
    )  # both runs in one tick
    root, db, mirror = bank
    sub = root / "fasta"
    scan(sub, db, mirror, log=lambda *_: None)
    victim = next(sub.glob("*.fa"))
    victim.unlink()
    assert scan(sub, db, mirror, log=lambda *_: None)["newly_missing"] == 1


def test_missing_scoped_to_root(bank, tmp_path):
    root, db, mirror = bank
    run(bank)
    other = tmp_path / "other"
    other.mkdir()
    shutil.copy(root / "pmc_gb" / "PMC-03204.gbk", other)
    assert scan(other, db, mirror, log=lambda *_: None)["newly_missing"] == 0
    idx = Index(db)
    assert idx.missing() == []
    idx.close()


def test_duplicates_and_ambiguous(bank):
    root, db, _ = bank
    run(bank)
    idx = Index(db)
    dup_groups = [names(g) for g in idx.duplicates()]
    assert sorted([CTRL, CTRL1]) in dup_groups
    amb = idx.ambiguous_ids()
    assert names(amb["PMC-03207"]) == sorted([CTRL, CTRL1])
    assert len(amb["PMC-00147"]) == 3  # clc/ + clc_gb_pairs/ (.clc and .gbk)
    assert "PMC-00029" in idx.clc_only()
    assert "PMC-00147" not in idx.clc_only()  # has a readable .gbk
    xref_hit = [r for r in idx.search("PMC-02483")]
    assert any(r["plasmid_id"] == "PMC-03405" for r in xref_hit)
    idx.close()


def test_mirror_skips_clc_and_keeps_versions(bank):
    root, db, mirror = bank
    run(bank)
    files = [p for p in mirror.rglob("*") if p.is_file()]
    assert files and not any(p.suffix == ".clc" for p in files)
    idx = Index(db)
    assert all(r["mirror_path"] is None for r in idx.files() if r["format"] == "clc")
    src = root / "pmc_gb" / "PMC-03204.gbk"
    old = Path(idx.get(str(src.resolve()))["mirror_path"])
    idx.close()
    assert (
        old.exists() and old.parent.name == "PMC" and old.name.startswith("PMC-03204__")
    )

    src.write_bytes(src.read_bytes().replace(b"LOCUS", b"LOCUS ", 1))
    st = src.stat()
    os.utime(src, ns=(st.st_atime_ns, st.st_mtime_ns + 10**9))
    assert run(bank)["mirrored_new"] == 1
    idx = Index(db)
    new = Path(idx.get(str(src.resolve()))["mirror_path"])
    idx.close()
    assert new != old and new.exists() and old.exists()
    assert sorted(p.name for p in old.parent.glob("PMC-03204__*")) == sorted(
        [old.name, new.name]
    )


def test_mirror_recreated_when_deleted(bank):
    root, db, mirror = bank
    run(bank)
    shutil.rmtree(mirror)
    counts = run(bank)
    assert counts["parsed"] == 0 and counts["mirrored_new"] > 0
    assert any(mirror.rglob("*.gb"))


@pytest.mark.parametrize(
    "folder,rng",
    [
        ("PMC-03001- PMC-03500", (3001, 3500)),
        ("PMC-04001-4500", (4001, 4500)),
        ("PMC-06500-7000", (6500, 7000)),
        ("PMC-00001 - PMC-00500", (1, 500)),
        ("pmc_05501 – pmc_06000", (5501, 6000)),
        ("PMC-03456-pBIN19-p35S(long)-GUS-tNOS", None),
        ("ENSA genbank", None),
        ("PMC-05000-4500", None),
    ],
)
def test_parse_bin(folder, rng):
    assert parse_bin(folder) == rng


def test_out_of_bin(bank):
    root, db, _ = bank
    good = root / "PMC-03001- PMC-03500"
    bad = root / "PMC-04001-4500" / "sub"
    good.mkdir()
    bad.mkdir(parents=True)
    shutil.copy(root / "pmc_gb" / "PMC-03204.gbk", good)
    shutil.copy(root / "pmc_gb" / "PMC-03204.gbk", bad / "PMC-03204 copy.gbk")
    run(bank)
    idx = Index(db)
    hits = idx.out_of_bin(str(root.resolve()))
    assert [(Path(r["path"]).name, f, lo, hi) for r, f, lo, hi in hits] == [
        ("PMC-03204 copy.gbk", "PMC-04001-4500", 4001, 4500)
    ]
    idx.close()
    assert bin_of(str(bad / "x.gb"))[0] == "PMC-04001-4500"


def test_skips_clutter_dotfiles_and_large(bank, monkeypatch):
    root, db, mirror = bank
    (root / "ws.clcinfo").write_text("x")
    (root / ".hidden.gb").write_text("LOCUS x")
    (root / ".dotdir").mkdir()
    (root / ".dotdir" / "PMC-09999.gb").write_text("LOCUS x")
    (root / "segments_1").write_text("x")
    import plasmidbank.scan as s

    monkeypatch.setattr(s, "MAX_BYTES", 20_000)
    counts = run(bank)
    assert counts["skipped"][".clcinfo"] == 1
    assert counts["skipped"]["dotfile"] >= 1
    assert counts["skipped"]["segments"] == 1
    assert counts["skipped"]["over MAX_BYTES"] > 0
    idx = Index(db)
    assert idx.by_id("PMC-09999") == []
    assert all(r["size"] <= 20_000 for r in idx.files())
    idx.close()


def test_refuses_to_write_inside_root(bank):
    from plasmidbank.scan import FatalError

    root, _, _ = bank
    with pytest.raises(FatalError):
        scan(root, root / "index.sqlite", root.parent / "m", log=lambda *_: None)
    assert not (root / "index.sqlite").exists()


def test_report_is_complete_document(bank):
    root, db, _ = bank
    run(bank)
    page = render_report(db, root)
    for s in (
        "<!DOCTYPE html>",
        '<html lang="en">',
        'charset="utf-8"',
        'name="viewport"',
        'name="color-scheme"',
        "<title>",
        "prefers-color-scheme",
        "@media print",
        "data:image/svg+xml",
        "Ambiguous IDs",
        "PMC-03207",
        "Out of bin",
        "CLC-only IDs",
    ):
        assert s in page, s
    assert "http://" not in page.replace("http://www.w3.org/2000/svg", "")
    assert "https://" not in page


def test_run_row(bank):
    root, db, _ = bank
    run(bank)
    idx = Index(db)
    r = idx.last_run(str(root.resolve()))
    assert r["finished"] and json.loads(r["counts"])["complete"] is True
    idx.close()


def test_excluded_folder_is_not_read(bank):
    root, db, _ = bank
    genome = (
        root / "Morex genome v1.0 data files "
    )  # the real folder name ends in a space
    genome.mkdir()
    (genome / "PMC-09999 not a plasmid.gb").write_bytes(
        (root / "pmc_gb" / "PMC-03204.gbk").read_bytes()
    )
    first = run(bank)
    assert first["skipped"].get("excluded folder") == 1
    rows = (
        Index(db)
        .db.execute("SELECT path FROM file WHERE path LIKE '%Morex%'")
        .fetchall()
    )
    assert rows == []


def _row(
    i,
    rel,
    pid,
    fmt="genbank",
    size=5000,
    length=5000,
    status="ok",
    sha="s",
    content=None,
    messages="[]",
):
    return {
        "id": i,
        "path": f"/O/bank/{rel}",
        "plasmid_id": pid,
        "format": fmt,
        "size": size,
        "length": None if fmt == "clc" else length,
        "status": status,
        "sha256": sha + str(i),
        "content_hash": content,
        "messages": messages,
    }


@pytest.mark.parametrize(
    "row, others, category, action",
    [
        (
            _row(1, "PMC0001-1000/PMC-00001.clc", "PMC-00001", "clc", 32),
            [],
            "empty",
            "map lost on the share",
        ),
        (
            _row(1, "PMC0001-1000/PMC-00002.clc", "PMC-00002", "clc", 32),
            [_row(2, "PMC0001-1000/PMC-00002.gbk", "PMC-00002")],
            "empty",
            "delete; PMC0001-1000/PMC-00002.gbk",
        ),
        (
            _row(1, "ENSA/EC1 x copy.gbk", "EC00001", content="h"),
            [_row(2, "ENSA/EC1 x.gbk", "EC00001", content="h")],
            "copy",
            "delete the copy",
        ),
        (
            _row(1, "ENSA/EC1 x (2).gbk", "EC00001", content="a"),
            [_row(2, "ENSA/EC1 x.gbk", "EC00001", content="b")],
            "copy",
            "someone who knows",
        ),
        (_row(1, "ENSA/25305.gbk", None), [], "no_id", "rename to EC25305"),
        (
            _row(1, "PMC5001-6000/PM-5435_CNGC15.gbk", None),
            [],
            "no_id",
            "check, then rename to PMC-05435",
        ),
        (
            _row(
                1, "qPCR amplicons/Lj1g.clc", None, "clc"
            ),
            [],
            "no_id",
            "move the folder",
        ),
        (
            _row(1, "PMC3001-4000/PMC-03777 x.clc.gbk", "PMC-03777"),
            [],
            "odd_name",
            "rename",
        ),
        (
            _row(
                1,
                "x/a.PDW",
                None,
                "unknown",
                status="error",
                messages='["unrecognised format"]',
            ),
            [],
            "unreadable",
            "move out of the bank",
        ),
    ],
)
def test_cleanup_overview_classifies(row, others, category, action):
    from plasmidbank.cleanup_report import classify

    items = classify([row, *others], {1: 3, 2: 3}, set(), "/O/bank")
    mine = [
        it
        for it in items
        if it["path"] == row["path"][len("/O/bank/") :] and it["category"] == category
    ]
    assert mine and mine[0]["action"].startswith(action)


def test_cleanup_overview_leaves_clean_files_out():
    from plasmidbank.cleanup_report import classify

    rows = [
        _row(1, "PMC0001-1000/PMC-00010_pX.gbk", "PMC-00010"),
        _row(2, "PMC0001-1000/PMC-00010.clc", "PMC-00010", "clc"),
    ]
    assert classify(rows, {1: 4}, set(), "/O/bank") == []
