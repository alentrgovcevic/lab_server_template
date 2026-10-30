"""Uploads and the add-only sync to the share. Temp folders only — never the real share.

The layout is the one of 2026-10-07: 1000-wide PMC bins named 'PMC7001-8000' and a flat 'ENSA/'.
The fixture also keeps the old, emptied folders, as on the real share, so the tests prove nothing
ever lands in them.
"""

from pathlib import Path

import pytest

from plasmidbank import uploads as up

SAMPLE = (
    Path(__file__).resolve().parents[2] / "examples" / "bank" / "pmc_gb" / "PMC-03204.gbk"
)


@pytest.fixture
def gb():
    if not SAMPLE.exists():
        pytest.skip("examples/bank/ not present")
    return SAMPLE.read_bytes()


@pytest.fixture
def o_drive(tmp_path):
    o = tmp_path / "Plasmid map bank"
    (o / "PMC6001-7000").mkdir(parents=True)
    (o / "PMC7001-8000").mkdir()
    (o / "ENSA").mkdir()
    # the pre-2026-10-07 folders, emptied but still there
    (o / "PMC-07001-7500").mkdir()
    (o / "PMC-06500-7000").mkdir()
    (o / "ENSA genbank" / "EC58401-EC58500").mkdir(parents=True)
    return o


def test_ids_are_validated_not_allocated():
    assert up.normalise_id("pmc 7383") == "PMC-07383"
    assert up.normalise_id("EC-58489") == "EC58489"
    for bad in ("", "7383", "PMC", "EC123", "foo"):
        with pytest.raises(ValueError):
            up.normalise_id(bad)


def test_store_names_file_and_refuses_existing_ids(tmp_path, gb):
    u = tmp_path / "uploads"
    p = up.store_upload(
        u, "anna", gb, "PMC-7383", "pL1 NFR1 / test", existing_ids=set()
    )
    assert p == u / "anna" / "PMC-07383_pL1_NFR1_test.gb"
    with pytest.raises(ValueError, match="already in the map bank"):
        up.store_upload(u, "anna", gb, "PMC-03204", "x", existing_ids={"PMC-03204"})
    with pytest.raises(ValueError, match="already uploaded by anna"):
        up.store_upload(u, "ben", gb, "PMC-07383", "x", existing_ids=set())
    with pytest.raises(ValueError, match="not a readable"):
        up.store_upload(u, "anna", b"hello", "PMC-07384", "x", existing_ids=set())
    assert not list((u / "anna").glob(".upload-*"))


def test_reupload_replaces_own_queued_file(tmp_path, gb):
    u = tmp_path / "uploads"
    up.store_upload(u, "anna", gb, "PMC-07383", "first", set())
    up.store_upload(u, "anna", gb, "PMC-07383", "second", set())
    assert [p.name for p in (u / "anna").iterdir()] == ["PMC-07383_second.gb"]


def test_bin_choice_ignores_the_old_folders(o_drive):
    # PMC-07001-7500 is narrower and covers 7383, but it is an old, dead folder
    assert up.bin_for(o_drive, "PMC-07383").name == "PMC7001-8000"
    assert up.bin_for(o_drive, "PMC-06999").name == "PMC6001-7000"
    assert up.bin_for(o_drive, "PMC-09999") is None
    assert up.bin_for(o_drive, "EC58489") == o_drive / "ENSA"  # EC: always flat


def test_bin_names_follow_the_convention():
    assert up.bin_name(*up.grid_bin(7383)) == "PMC7001-8000"
    assert up.bin_name(*up.grid_bin(7000)) == "PMC6001-7000"
    assert up.bin_name(*up.grid_bin(1)) == "PMC0001-1000"
    assert up.bin_name(*up.grid_bin(1000)) == "PMC0001-1000"
    assert up.bin_name(*up.grid_bin(10001)) == "PMC10001-11000"


def test_only_standard_new_style_bins_are_in_the_series(o_drive):
    (o_drive / "PMC5000-7000").mkdir()  # odd size
    assert [d.name for _, _, d in up.series(o_drive)] == [
        "PMC6001-7000",
        "PMC7001-8000",
    ]


def test_number_without_a_bin_gets_only_its_own_bin(tmp_path, gb, o_drive):
    u = tmp_path / "uploads"
    up.store_upload(u, "anna", gb, "PMC-09123", "a", set())
    r = up.sync(u, o_drive, write=True, alert=lambda m: pytest.fail(m))
    assert [Path(d).name for d in r["folders"]] == ["PMC9001-10000"]
    assert (o_drive / "PMC9001-10000" / "PMC-09123_a.gb").exists()


def test_file_in_second_to_last_bin_creates_five_more(tmp_path, gb, o_drive):
    u = tmp_path / "uploads"
    alerts = []
    up.store_upload(
        u, "anna", gb, "PMC-06999", "a", set()
    )  # PMC6001-7000 is second-to-last
    dry = up.sync(u, o_drive, write=False, alert=alerts.append)
    want = [
        "PMC8001-9000",
        "PMC9001-10000",
        "PMC10001-11000",
        "PMC11001-12000",
        "PMC12001-13000",
    ]
    assert [Path(d).name for d in dry["folders"]] == want and not (
        o_drive / want[0]
    ).exists()
    r = up.sync(u, o_drive, write=True, alert=alerts.append)
    assert [Path(d).name for d in r["folders"]] == want
    assert all((o_drive / n).is_dir() for n in want) and len(alerts) == 1
    up.store_upload(
        u, "ben", gb, "PMC-06998", "b", set()
    )  # that bin is no longer second-to-last
    r = up.sync(u, o_drive, write=True, alert=alerts.append)
    assert r["folders"] == [] and len(alerts) == 1


def test_dry_run_writes_nothing(tmp_path, gb, o_drive):
    u = tmp_path / "uploads"
    up.store_upload(u, "anna", gb, "PMC-07383", "a", set())
    r = up.sync(u, o_drive, write=False, alert=lambda m: pytest.fail(m))
    assert len(r["copied"]) == 1
    assert not list((o_drive / "PMC7001-8000").iterdir())
    assert not (u / "ledger.json").exists()


@pytest.mark.parametrize("o_up", [True, False])
def test_first_run_with_no_uploads_folder(tmp_path, o_drive, o_up):
    # the VM's first hourly run, 2026-10-22: nobody had uploaded, so /state/uploads did not exist
    u = tmp_path / "uploads"
    dest = o_drive if o_up else tmp_path / "unmounted"
    r = up.sync(u, dest, write=True, alert=lambda m: None)
    assert r["queued"] == 0 and r["reachable"] is o_up
    assert (u / "ledger.json").exists()


def test_sync_copies_once_and_never_overwrites(tmp_path, gb, o_drive):
    u = tmp_path / "uploads"
    alerts = []
    up.store_upload(u, "anna", gb, "PMC-07383", "a", set())
    up.store_upload(u, "anna", gb, "EC58489", "b", set())
    r = up.sync(u, o_drive, write=True, alert=alerts.append)
    assert len(r["copied"]) == 2 and not r["refused"] and r["folders"] == []
    assert (o_drive / "PMC7001-8000" / "PMC-07383_a.gb").read_bytes() == gb
    assert (o_drive / "ENSA" / "EC58489_b.gb").exists()
    # nothing went into the old folders
    assert not list((o_drive / "PMC-07001-7500").iterdir())
    assert not list((o_drive / "ENSA genbank" / "EC58401-EC58500").iterdir())
    assert {x["status"] for x in up.status(u)} == {"synced"}
    r = up.sync(u, o_drive, write=True, alert=alerts.append)
    assert r["queued"] == 0 and not r["copied"]
    assert alerts == []


def test_name_or_id_clash_on_o_is_refused_and_alerted_once(tmp_path, gb, o_drive):
    u = tmp_path / "uploads"
    alerts = []
    clash = o_drive / "PMC7001-8000" / "PMC-07383_a.gb"
    clash.write_text("someone else's file")
    (o_drive / "PMC7001-8000" / "PMC-07390 other name.gbk").write_text("x")
    (o_drive / "ENSA" / "EC58489_pL1V-old.clc").write_text(
        "x"
    )  # an EC already in the flat folder
    up.store_upload(u, "anna", gb, "PMC-07383", "a", set())
    up.store_upload(u, "ben", gb, "PMC-07390", "mine", set())
    up.store_upload(u, "ben", gb, "EC58489", "new", set())
    for _ in range(2):
        r = up.sync(u, o_drive, write=True, alert=alerts.append)
        assert len(r["refused"]) == 3 and not r["copied"]
    assert clash.read_text() == "someone else's file"
    assert len(alerts) == 3  # one per refused file, not one per run
    assert {x["status"] for x in up.status(u)} == {"refused"}


def test_identical_file_already_on_o_counts_as_synced(tmp_path, gb, o_drive):
    u = tmp_path / "uploads"
    up.store_upload(u, "anna", gb, "PMC-07383", "a", set())
    (o_drive / "PMC7001-8000" / "PMC-07383_a.gb").write_bytes(gb)
    r = up.sync(u, o_drive, write=True, alert=lambda m: None)
    assert r["already"] and up.status(u)[0]["status"] == "synced"


def test_unreachable_o_queues_and_alerts_on_transitions_only(tmp_path, gb, o_drive):
    u = tmp_path / "uploads"
    alerts = []
    up.store_upload(u, "anna", gb, "PMC-07383", "a", set())
    gone = tmp_path / "not mounted"
    gone.mkdir()  # an empty mount point, as docker sees it when the share is not mounted
    for _ in range(3):
        r = up.sync(u, gone, write=True, alert=alerts.append)
        assert not r["reachable"] and not r["copied"]
    assert len(alerts) == 1 and "not reachable" in alerts[0]
    r = up.sync(u, o_drive, write=True, alert=alerts.append)
    assert len(r["copied"]) == 1
    assert len(alerts) == 2 and "reachable again" in alerts[1]


def test_old_layout_alone_is_not_reachable(tmp_path, gb):
    """A share with only 'ENSA genbank/' (before 2026-10-07) must not receive files."""
    old = tmp_path / "old bank"
    (old / "ENSA genbank").mkdir(parents=True)
    u = tmp_path / "uploads"
    up.store_upload(u, "anna", gb, "EC58489", "b", set())
    r = up.sync(u, old, write=True, alert=lambda m: None)
    assert not r["reachable"] and not r["copied"]


def test_queued_upload_can_be_deleted_but_not_after_syncing(tmp_path, gb, o_drive):
    u = tmp_path / "uploads"
    up.store_upload(u, "anna", gb, "PMC-07383", "a", set())
    up.store_upload(u, "ben", gb, "PMC-07390", "b", set())
    with pytest.raises(PermissionError, match="uploaded by anna"):
        up.delete_queued(u, "ben", "PMC-07383")
    with pytest.raises(ValueError, match="no upload"):
        up.delete_queued(u, "anna", "PMC-09999")
    up.delete_queued(u, "anna", "PMC-07383")
    assert [r["plasmid_id"] for r in up.status(u)] == ["PMC-07390"]
    up.sync(u, o_drive, write=True, alert=lambda m: None)
    with pytest.raises(ValueError, match="already on the share"):
        up.delete_queued(u, "ben", "PMC-07390")
    assert (o_drive / "PMC7001-8000" / "PMC-07390_b.gb").exists()  # the copy on the share is untouched


def _gbtext(seq: str, name: str = "x") -> bytes:
    import io
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    buf = io.StringIO()
    SeqIO.write(SeqRecord(Seq(seq), id=name, name=name, annotations={"molecule_type": "DNA"}), buf, "genbank")
    return buf.getvalue().encode()


LATER = 10**10  # "now" far in the future: every file on the share counts as stable


@pytest.fixture
def ws(tmp_path, o_drive):
    """ann has one draft in My workspace (folder RinRK1), already in sync with the share."""
    import json as _json
    from plasmidbank.projects import Projects

    state = tmp_path / "state"
    pj = Projects(state)
    (state / "users.json").write_text(_json.dumps({"ann": {}, "bob": {}}))
    pid = pj.personal("ann")
    iid = pj.add_draft(pid, "ann", "a.gb", _gbtext("ACGT" * 50, "a"), "exon test", folder="RinRK1")
    alerts = []
    sync = lambda **kw: up.sync_workspaces(state, o_drive, write=True, alert=alerts.append, now=kw.pop("now", LATER), **kw)
    sync()
    o = o_drive / "Workspaces" / "ann" / "RinRK1" / "exon test.gb"
    assert o.read_bytes() == _gbtext("ACGT" * 50, "a")
    return {"state": state, "pj": pj, "pid": pid, "iid": iid, "o": o, "ws": o_drive / "Workspaces", "sync": sync,
            "alerts": alerts, "app": lambda: pj.draft_path(pid, "ann", iid)}


def test_workspace_user_folders_exist_for_every_account(ws):
    assert (ws["ws"] / "bob").is_dir() and (ws["ws"] / "ann").is_dir()


def test_workspace_app_edit_is_pushed(ws):
    ws["pj"].replace_draft(ws["pid"], "ann", ws["iid"], _gbtext("GGGG" * 50))
    assert ws["sync"]()["written"] and ws["o"].read_bytes() == _gbtext("GGGG" * 50)
    assert ws["sync"]()["written"] == []  # nothing changed, nothing written


def test_workspace_edit_on_o_is_pulled_once_stable(ws):
    ws["o"].write_bytes(_gbtext("TTTT" * 50))
    import time
    r = ws["sync"](now=time.time())  # just saved on the share: maybe still being written
    assert r["waiting"] and ws["app"]().read_bytes() == _gbtext("ACGT" * 50, "a")
    assert ws["sync"]()["pulled"] and ws["app"]().read_bytes() == _gbtext("TTTT" * 50)


def test_workspace_both_edited_keeps_both(ws):
    ws["pj"].replace_draft(ws["pid"], "ann", ws["iid"], _gbtext("GGGG" * 50))
    ws["o"].write_bytes(_gbtext("TTTT" * 50))
    r = ws["sync"]()
    assert r["conflicts"] and ws["o"].read_bytes() == _gbtext("GGGG" * 50) and ws["alerts"]
    items = ws["pj"].get(ws["pid"], "ann")["items"]
    copy = next(i for i in items if "changed on the share" in i["label"])
    assert copy["folder"] == "RinRK1" and ws["pj"].draft_path(ws["pid"], "ann", copy["id"]).read_bytes() == _gbtext("TTTT" * 50)
    ws["sync"]()  # the next run puts the share version on the share under its own name
    assert (ws["o"].parent / f"{copy['label']}.gb").read_bytes() == _gbtext("TTTT" * 50)


def test_workspace_removed_in_app_goes_to_deleted_on_o(ws):
    ws["pj"].remove_item(ws["pid"], "ann", ws["iid"])
    ws["sync"]()
    assert not ws["o"].exists() and [p.read_bytes() for p in (ws["ws"] / "ann" / "_deleted").iterdir()] == [_gbtext("ACGT" * 50, "a")]


def test_workspace_removed_on_o_leaves_the_app_with_a_trash_copy(ws):
    ws["o"].unlink()
    assert ws["sync"]()["removed_in_app"] and ws["pj"].get(ws["pid"], "ann")["items"] == []
    assert [p.read_bytes() for p in (ws["state"] / "workspaces_trash").iterdir()] == [_gbtext("ACGT" * 50, "a")]


def test_workspace_removed_on_o_but_edited_in_app_is_pushed_again(ws):
    ws["o"].unlink()
    ws["pj"].replace_draft(ws["pid"], "ann", ws["iid"], _gbtext("GGGG" * 50))
    ws["sync"]()
    assert ws["o"].read_bytes() == _gbtext("GGGG" * 50)


def test_workspace_new_file_on_o_is_imported_into_its_folder(ws):
    new = ws["ws"] / "ann" / "faba" / "pFaba1.fasta"
    new.parent.mkdir()
    new.write_text(">pFaba1\nACGTACGTAA\n")
    assert ws["sync"]()["imported"] == ["ann/faba/pFaba1.fasta"]
    it = next(i for i in ws["pj"].get(ws["pid"], "ann")["items"] if i["label"] == "pFaba1")
    assert it["folder"] == "faba"
    assert ws["sync"]()["imported"] == [] and new.read_text() == ">pFaba1\nACGTACGTAA\n"  # once, untouched


def test_workspace_unreadable_file_on_o_is_reported_once_and_left_alone(ws):
    clc = ws["ws"] / "ann" / "PMC-01118.clc"
    clc.write_bytes(b"CLC\x00\xac\xed\x00\x05 binary")
    assert ws["sync"]()["unreadable"] == ["ann/PMC-01118.clc"] and len(ws["alerts"]) == 1
    assert ws["sync"]()["unreadable"] == [] and len(ws["alerts"]) == 1 and clc.exists()


def test_workspace_hand_placed_file_is_never_overwritten(ws):
    hand = ws["ws"] / "ann" / "RinRK1" / "linear b.gb"
    hand.write_bytes(_gbtext("CCCC" * 10, "hand"))
    iid = ws["pj"].add_draft(ws["pid"], "ann", "b.gb", _gbtext("AAAA" * 10, "b"), "linear b", folder="RinRK1")
    ws["sync"]()
    assert hand.read_bytes() == _gbtext("CCCC" * 10, "hand") and (hand.parent / f"linear b_{iid}.gb").exists()


def test_workspace_move_in_app_moves_on_o(ws):
    ws["pj"].move(ws["pid"], "ann", ws["iid"], folder="RinRK1/done", label="exon v2")
    assert ws["sync"]()["moved"] and (ws["ws"] / "ann" / "RinRK1" / "done" / "exon v2.gb").exists() and not ws["o"].exists()


def test_workspace_sync_waits_for_o_and_for_the_lock(ws, tmp_path):
    import fcntl
    empty = tmp_path / "mnt"
    empty.mkdir()
    assert up.sync_workspaces(ws["state"], empty, write=True)["reachable"] is False and not list(empty.iterdir())
    with open(ws["state"] / "workspaces.lock", "a") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        assert ws["sync"]()["busy"]
