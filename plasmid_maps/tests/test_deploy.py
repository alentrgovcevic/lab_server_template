import importlib
from pathlib import Path

import pytest

SAMPLE = Path(__file__).resolve().parents[2] / "examples" / "bank"


@pytest.fixture
def auth(tmp_path, monkeypatch):
    monkeypatch.setenv("PLASMIDBANK_AUTH", "local")
    monkeypatch.setenv("PLASMIDBANK_USERS", str(tmp_path / "users.json"))
    monkeypatch.setenv("PLASMIDBANK_SECRET", "x" * 40)
    import plasmidbank.app.auth as a

    return importlib.reload(a)


def test_passwords_and_cookies(auth):
    auth.set_password("alice", "a-long-password")
    assert auth.check_password("alice", "a-long-password")
    assert not auth.check_password("alice", "wrong-password")
    assert not auth.check_password("nobody", "a-long-password")
    c = auth.make_cookie("alice")
    assert auth.read_cookie(c) == "alice"
    payload, sig = c.rsplit(".", 1)
    assert auth.read_cookie(payload + "." + "0" * len(sig)) is None  # forged signature
    assert auth.read_cookie("garbage") is None
    assert auth.current_user({}, {auth.COOKIE: c}) == "alice"
    assert auth.current_user({}, {}) is None
    users_file = Path(auth.USERS_FILE)
    assert oct(users_file.stat().st_mode & 0o777) == "0o600"


def test_removed_user_cookie_stops_working(auth):
    auth.set_password("bob", "another-long-pw")
    c = auth.make_cookie("bob")
    users = auth.load_users()
    users.pop("bob")
    auth.save_users(users)
    assert auth.read_cookie(c) is None


def test_bad_names_and_short_passwords_refused(auth):
    with pytest.raises(ValueError):
        auth.set_password("../evil", "a-long-password")
    with pytest.raises(ValueError):
        auth.set_password("carol", "short")


def test_username_for_lab_names(auth):
    assert auth.username_for("Mærta Maria Bjerg Hansen") == "maerta.hansen"
    assert auth.username_for("Ole Mønsted Jensen") == "ole.jensen"
    assert auth.username_for("Ana Maria Lopez Rodríguez") == "ana.rodriguez"
    assert auth.username_for("Carl M. Example") == "carl.example"
    assert auth.username_for("Sara Double-Barrel") == "sara.double-barrel"
    assert auth.username_for("Jose Acenté") == "jose.acente"
    assert auth.username_for("  Eve   Spaced ") == "eve.spaced"


def test_bulk_create_one_time_passwords(auth):
    rows = auth.bulk_create(["Alex Smith", "", "Alex Marie Smith", "Jane Doe"])
    names = [r[1] for r in rows]
    assert names == ["alex.smith", "alex.smith2", "jane.doe"]  # a clash gets a number
    full, user, pw = rows[2]
    assert auth.check_password(user, pw) and auth.must_change(user)
    assert auth.load_users()[user]["full_name"] == "Jane Doe"
    again = auth.bulk_create(["Jane Doe"])  # re-running never resets anyone
    assert again == [("Jane Doe", "jane.doe", "exists")]
    assert auth.check_password(user, pw)
    assert "bulk add" in auth.AUDIT_FILE.read_text()


def test_safe_name_never_escapes():
    from plasmidbank.app.auth import safe_name

    assert safe_name("Ann.T@lab.example.org") == "ann.t_lab.example.org"
    assert "/" not in safe_name("../../etc/passwd") and not safe_name(
        "../x"
    ).startswith(".")


def test_header_mode(monkeypatch):
    monkeypatch.setenv("PLASMIDBANK_AUTH", "header")
    import plasmidbank.app.auth as a

    a = importlib.reload(a)
    assert a.current_user({"x-forwarded-email": "Someone@example.org"}, {}) == "someone_example.org"
    assert a.current_user({}, {}) is None
    monkeypatch.setenv("PLASMIDBANK_AUTH", "none")
    importlib.reload(a)


def test_bank_from_index_uses_mirror_and_keeps_original_paths(tmp_path):
    from plasmidbank.bank import Bank
    from plasmidbank.scan import scan

    db, mirror = tmp_path / "index.sqlite", tmp_path / "mirror"
    scan(str(SAMPLE), db=db, mirror_dir=mirror)
    bank = Bank.from_index(db, str(SAMPLE))
    direct = Bank.load(SAMPLE)
    assert len(bank.maps) == len(direct.maps)
    m = next(x for x in bank.maps if x.plasmid_id == "PMC-07001")
    assert m.path.startswith(str(mirror)) and m.display_path.startswith(str(SAMPLE))
    assert m.seq and m.features
    clc = [x for x in bank.maps if x.format == "clc"]
    assert clc and all(not x.seq and x.status == "clc_binary" for x in clc)
    assert bank.search("pCAMBIA1301")  # matches the original filename, not the mirror name


def test_blast_database_files_are_skipped(tmp_path):
    import shutil
    from plasmidbank.bank import Bank
    db = tmp_path / "Clones_blast_db.db"
    db.mkdir()
    for suffix in (".nhr", ".nin", ".nsq", ".nsd", ".nsi"):
        (db / f"Clones_blast_db{suffix}").write_bytes(b"\x00\x01binary")
    (db / "formatdb.log").write_text("formatdb run")
    shutil.copy(next((SAMPLE / "pmc_gb").glob("PMC-07001*")), tmp_path / "PMC-07001.gbk")
    bank = Bank.load(tmp_path)
    assert [m.plasmid_id for m in bank.maps] == ["PMC-07001"]


def test_app_reloads_the_bank_when_the_index_changes(tmp_path, monkeypatch):
    # 2026-10-27: the scan runs daily at 07:00 and the app picks the new index up without a restart
    import os
    from plasmidbank.scan import scan

    src, db, mirror = tmp_path / "o", tmp_path / "index.sqlite", tmp_path / "mirror"
    src.mkdir()
    (src / "PMC-00001_a.fa").write_text(">a\nACGTACGTAC\n")
    scan(str(src), db=db, mirror_dir=mirror)
    for k, v in {"PLASMIDBANK_INDEX": db, "PLASMIDBANK_ROOT": src, "PLASMIDBANK_STATE": tmp_path / "state",
                 "PLASMIDBANK_UPLOADS": tmp_path / "uploads", "PLASMIDBANK_OUTBOX": tmp_path / "outbox",
                 "PLASMIDBANK_OLIGOS": tmp_path / "none.tab", "PLASMIDBANK_AUTH": "none"}.items():
        monkeypatch.setenv(k, str(v))
    from plasmidbank.app import server
    s = importlib.reload(server)
    maps = s.bank.maps
    assert [m.plasmid_id for m in maps] == ["PMC-00001"] and not s.reload_if_index_changed()
    (src / "PMC-00002_b.fa").write_text(">b\nGGGGCCCCAA\n")
    scan(str(src), db=db, mirror_dir=mirror)
    os.utime(db, (s._index_stamp() + 5,) * 2)
    assert not s.reload_if_index_changed()  # first look: maybe still being written
    assert s.reload_if_index_changed()  # unchanged since: reload
    assert s.bank.maps is maps and sorted(m.plasmid_id for m in maps) == ["PMC-00001", "PMC-00002"]
    assert any(r["plasmid_id"] == "PMC-00002" for r in s.maps(q="PMC-00002"))
    monkeypatch.delenv("PLASMIDBANK_INDEX")
    monkeypatch.delenv("PLASMIDBANK_ROOT")
