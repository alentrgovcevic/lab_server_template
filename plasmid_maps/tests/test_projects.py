"""Projects: members only, links survive restarts, drafts promote into the bank."""

import importlib
import json
import socket
import threading
import time
import urllib.error
import urllib.request
import uuid

import pytest
import uvicorn

from tests.test_bench import _rand


def _gb(seq: str, name: str) -> str:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    import io

    rec = SeqRecord(Seq(seq), id=name, name=name, description=name,
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    buf = io.StringIO()
    SeqIO.write(rec, buf, "genbank")
    return buf.getvalue()


@pytest.fixture
def app(tmp_path, monkeypatch):
    data = tmp_path / "bank"
    data.mkdir()
    bb, body = _rand(2000, 21), _rand(300, 22)
    (data / "EC00011_part.gb").write_text(
        _gb(bb[:600] + "GGTCTCA" + "GGAG" + body + "CGCT" + "T" + "GAGACC" + bb[600:1200], "part"))
    (data / "PMC-00012_acceptor.gb").write_text(
        _gb(bb[1200:1800] + "GGAG" + "TGAGACC" + _rand(200, 23) + "GGTCTCA" + "CGCT" + bb[1800:], "acc"))
    for k, v in {"PLASMIDBANK_DATA": data, "PLASMIDBANK_STATE": tmp_path / "state",
                 "PLASMIDBANK_UPLOADS": tmp_path / "uploads", "PLASMIDBANK_OUTBOX": tmp_path / "outbox",
                 "PLASMIDBANK_OLIGOS": tmp_path / "none.tab", "PLASMIDBANK_AUTH": "header",
                 "PLASMIDBANK_USERS": tmp_path / "state" / "users.json"}.items():
        monkeypatch.setenv(k, str(v))
    monkeypatch.delenv("PLASMIDBANK_INDEX", raising=False)
    from plasmidbank.app import auth, server

    importlib.reload(auth)
    s = importlib.reload(server)
    sock = socket.socket()
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()
    srv = uvicorn.Server(uvicorn.Config(s.app, host="127.0.0.1", port=port, log_level="warning"))
    th = threading.Thread(target=srv.run, daemon=True)
    th.start()
    while not srv.started:
        time.sleep(0.02)
    yield s, Client(f"http://127.0.0.1:{port}")
    srv.should_exit = True
    th.join(5)
    monkeypatch.setenv("PLASMIDBANK_AUTH", "none")
    importlib.reload(auth)


class Response:
    def __init__(self, status, body):
        self.status_code, self.text = status, body.decode("utf-8", "replace")

    def json(self):
        return json.loads(self.text)


class Client:
    """A few lines of urllib instead of TestClient, which needs a package the venv does not have."""

    def __init__(self, base):
        self.base = base

    def _do(self, method, url, headers=None, json_body=None, files=None):
        headers = dict(headers or {})
        data = None
        if json_body is not None:
            data = json.dumps(json_body).encode()
            headers["Content-Type"] = "application/json"
        elif files:
            b = uuid.uuid4().hex
            parts = []
            for field, items in files.items():
                for fname, content in items if isinstance(items, list) else [items]:
                    body = content if isinstance(content, bytes) else content.encode()
                    if fname is None:  # a plain form field
                        parts.append(f'--{b}\r\nContent-Disposition: form-data; name="{field}"\r\n\r\n'.encode() + body + b"\r\n")
                        continue
                    parts.append(f'--{b}\r\nContent-Disposition: form-data; name="{field}"; filename="{fname}"\r\n'
                                 f"Content-Type: application/octet-stream\r\n\r\n".encode() + body + b"\r\n")
            data = b"".join(parts) + f"--{b}--\r\n".encode()
            headers["Content-Type"] = f"multipart/form-data; boundary={b}"
        req = urllib.request.Request(self.base + url, data=data, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return Response(r.status, r.read())
        except urllib.error.HTTPError as e:
            return Response(e.code, e.read())

    def get(self, url, headers=None):
        return self._do("GET", url, headers)

    def delete(self, url, headers=None):
        return self._do("DELETE", url, headers)

    def post(self, url, json=None, files=None, headers=None):
        return self._do("POST", url, headers, json, files)


def _as(user):
    return {"x-forwarded-user": user}


def _key(client, pid):
    return next(m["key"] for m in client.get("/api/maps", headers=_as("ann")).json() if m["plasmid_id"] == pid)


def test_dashboard_and_bank_routes(app):
    _, c = app
    assert "<h1>Lab server</h1>" in c.get("/", headers=_as("ann")).text
    assert "Plasmid map bank" in c.get("/bank", headers=_as("ann")).text
    assert "<h1>Changelog</h1>" in c.get("/changelog", headers=_as("ann")).text
    assert c.post("/api/feedback", json={"text": "a wish"}, headers=_as("ann")).status_code == 200
    assert c.get("/api/feedback", headers=_as("bob")).json()[0]["text"] == "a wish"


def test_project_is_private_until_the_owner_invites(app):
    _, c = app
    pid = c.post("/api/projects", json={"name": "NFR1"}, headers=_as("ann")).json()["id"]
    assert c.get(f"/api/projects/{pid}", headers=_as("bob")).status_code == 403
    assert c.post(f"/api/projects/{pid}/members", json={"user": "bob"}, headers=_as("bob")).status_code == 403
    assert c.post(f"/api/projects/{pid}/members", json={"user": "Bob"}, headers=_as("ann")).status_code == 200
    assert c.get(f"/api/projects/{pid}", headers=_as("bob")).status_code == 200
    assert [p["id"] for p in c.get("/api/projects", headers=_as("carl")).json()] == []
    assert c.delete(f"/api/projects/{pid}", headers=_as("bob")).status_code == 403


def test_links_are_stored_by_id_and_golden_gate_takes_project_items(app):
    s, c = app
    pid = c.post("/api/projects", json={"name": "GG"}, headers=_as("ann")).json()["id"]
    for plasmid in ("EC00011", "PMC-00012"):
        c.post(f"/api/projects/{pid}/refs", json={"key": _key(c, plasmid)}, headers=_as("ann"))
    items = c.get(f"/api/projects/{pid}", headers=_as("ann")).json()["items"]
    assert {i["plasmid_id"] for i in items} == {"EC00011", "PMC-00012"} and all(i["key"] is not None for i in items)
    s.bank.maps.reverse()  # keys move (as after a restart); the links still resolve by ID
    s._bank_changed()
    items = c.get(f"/api/projects/{pid}", headers=_as("ann")).json()["items"]
    assert all(s.bank.maps[i["key"]].plasmid_id == i["plasmid_id"] for i in items)
    r = c.post("/api/assemble/golden_gate", json={"inputs": [{"project": pid, "item": i["id"]} for i in items],
                                                  "enzyme": "BsaI"}, headers=_as("ann")).json()
    assert len(r["products"]) == 1
    sites = {j["overhang"]: j["site"] for j in r["products"][0]["junctions"]}
    prod = r["products"][0]
    for j in prod["junctions"]:  # the junction zoom slices the product at pos
        assert prod["seq"][j["pos"]:j["pos"] + len(j["overhang"])] == j["overhang"]
    assert "start of promoter" in sites["GGAG"]
    # another user cannot reach the items through the assembler
    bad = c.post("/api/assemble/golden_gate", json={"inputs": [{"project": pid, "item": i["id"]} for i in items]},
                 headers=_as("carl"))
    assert bad.status_code == 403


def test_draft_is_members_only_and_promotion_replaces_it_with_a_link(app):
    s, c = app
    pid = c.post("/api/projects", json={"name": "drafts"}, headers=_as("ann")).json()["id"]
    gb = _gb(_rand(1500, 31), "newconstruct")
    iid = c.post(f"/api/projects/{pid}/drafts", files={"file": ("new.gb", gb)}, headers=_as("ann")).json()["id"]
    assert not any(m["name"] == "newconstruct" for m in c.get("/api/maps", headers=_as("ann")).json())
    assert c.get(f"/api/projects/{pid}/items/{iid}/file", headers=_as("bob")).status_code == 403
    assert c.get(f"/api/projects/{pid}/items/{iid}/file", headers=_as("ann")).text.startswith("LOCUS")
    r = c.post(f"/api/projects/{pid}/items/{iid}/promote", json={"plasmid_id": "PMC-07001", "construct_name": "new"},
               headers=_as("ann"))
    assert r.status_code == 200, r.text
    it = c.get(f"/api/projects/{pid}", headers=_as("ann")).json()["items"][0]
    assert it["kind"] == "ref" and it["plasmid_id"] == "PMC-07001" and it["key"] is not None
    assert not list((s.STATE / "projects" / str(pid)).glob("*"))  # the draft file is gone


def test_upload_without_a_lab_id_is_a_private_draft(app):
    s, c = app
    gb = _gb(_rand(1200, 41), "intermediate")
    drafts = []
    for tag in ("Test", "TMP-3", ""):
        r = c.post("/api/upload", files={"file": ("x.gb", gb), "plasmid_id": (None, tag),
                                         "construct_name": (None, "step1")}, headers=_as("ann"))
        assert r.status_code == 200, r.text
        drafts.append(r.json())
    assert all(d["draft"] for d in drafts) and len({d["project"] for d in drafts}) == 1  # one My drafts
    pid = drafts[0]["project"]
    assert [i["label"] for i in c.get(f"/api/projects/{pid}", headers=_as("ann")).json()["items"]] == [
        "Test step1", "TMP-3 step1", "step1"]
    assert not any(m["name"] == "intermediate" for m in c.get("/api/maps", headers=_as("ann")).json())
    assert not (s.UPLOADS / "ann").exists()  # nothing queued for the sync to the share
    assert c.get(f"/api/projects/{pid}", headers=_as("bob")).status_code == 403
    # a typo of a real number stays an error rather than quietly missing O
    for typo in ("EC123", "PMC-1234567"):
        r = c.post("/api/upload", files={"file": ("x.gb", gb), "plasmid_id": (None, typo)}, headers=_as("ann"))
        assert r.status_code == 400, typo


def test_synced_upload_deleted_on_o_leaves_the_bank(tmp_path, monkeypatch):
    # 2026-10-22: PMC-09999 was synced to the share, then deleted there; the bank kept showing the upload copy
    data, uploads = tmp_path / "bank", tmp_path / "uploads"
    (uploads / "ann").mkdir(parents=True)
    data.mkdir()
    for pid in ("PMC-09999", "PMC-09998"):
        (uploads / "ann" / f"{pid}_t.gb").write_text(_gb(_rand(300, 51), pid))
    (uploads / "ledger.json").write_text(json.dumps({"files": {"ann/PMC-09999_t.gb": {
        "status": "synced", "dest": "/data/bank-rw/PMC9001-10000/PMC-09999_t.gb"}}}))  # not on the share any more
    for k, v in {"PLASMIDBANK_DATA": data, "PLASMIDBANK_STATE": tmp_path / "state", "PLASMIDBANK_UPLOADS": uploads,
                 "PLASMIDBANK_OUTBOX": tmp_path / "outbox", "PLASMIDBANK_OLIGOS": tmp_path / "none.tab",
                 "PLASMIDBANK_AUTH": "none"}.items():
        monkeypatch.setenv(k, str(v))
    monkeypatch.delenv("PLASMIDBANK_INDEX", raising=False)
    from plasmidbank.app import server
    s = importlib.reload(server)
    assert [m.plasmid_id for m in s.bank.maps] == ["PMC-09998"]  # the queued one stays, the deleted one goes
    assert (uploads / "ann" / "PMC-09999_t.gb").exists()  # kept on the server as a backup


def test_personal_project_is_added_to_an_old_store(tmp_path):
    import sqlite3
    from plasmidbank.projects import Projects

    old = sqlite3.connect(tmp_path / "projects.sqlite")  # the schema before 2026-10-22
    old.execute("CREATE TABLE projects (id INTEGER PRIMARY KEY, name TEXT NOT NULL, "
                "description TEXT NOT NULL DEFAULT '', owner TEXT NOT NULL, created TEXT NOT NULL)")
    old.execute("INSERT INTO projects (name, owner, created) VALUES ('NFR1', 'ann', '2026-10-06')")
    old.commit()
    old.close()
    pj = Projects(tmp_path)
    pid = pj.personal("ann")
    assert pid != 1 and pj.personal("ann") == pid and pj.personal("bob") != pid
    assert [p["name"] for p in pj.for_user("ann")] == ["My workspace"]  # ann was never a member of NFR1 here


def test_old_my_drafts_becomes_my_workspace_with_folders(tmp_path):
    import sqlite3
    from plasmidbank.projects import Projects

    pj = Projects(tmp_path)
    pid = pj.personal("ann")
    iid = pj.add_draft(pid, "ann", "a.gb", b"LOCUS x", "a")
    pj.db.execute("UPDATE projects SET name = 'My drafts' WHERE id = ?", (pid,))  # a store from before 2026-10-27
    pj.db.commit()
    pj.db.close()
    pj = Projects(tmp_path)
    assert pj.get(pid, "ann")["name"] == "My workspace"
    pj.move(pid, "ann", iid, folder=" RinRK1 / ../constructs/", label="pRinRK1 v2")
    it = pj.item(pid, "ann", iid)
    assert (it["folder"], it["label"]) == ("RinRK1/constructs", "pRinRK1 v2")
    shared = pj.create("ann", "NFR1")
    with pytest.raises(ValueError):
        pj.add_draft(shared, "ann", "b.gb", b"LOCUS y", "b", folder="sub")  # folders are for My workspace
    with pytest.raises(ValueError):
        pj.move(pid, "ann", iid, folder="a/b/c/d/e")


@pytest.mark.parametrize("old,new,expect", [
    ("Exp12", "Exp 13", {"a": "Exp 13", "b": "Exp 13/L0", "c": "Exp12x", "d": ""}),
    ("Exp12/L0", "L0 modules", {"a": "Exp12", "b": "L0 modules", "c": "Exp12x", "d": ""}),
    ("Exp12", "Done/Exp12", {"a": "Done/Exp12", "b": "Done/Exp12/L0", "c": "Exp12x", "d": ""}),
    ("Exp_2", "x", KeyError),           # '_' is not a LIKE wildcard: 'Exp12/L0' must not match 'Exp_2/%'
    ("Exp12", "Exp12/inside", ValueError),
    ("Exp12", "a/b/c/d", ValueError),   # Exp12/L0 would go 5 deep
])
def test_rename_workspace_folder_moves_its_subfolders(tmp_path, old, new, expect):
    from plasmidbank.projects import Projects

    pj = Projects(tmp_path)
    pid = pj.personal("ann")
    ids = {n: pj.add_draft(pid, "ann", f"{n}.gb", b"LOCUS x", n, folder=f)
           for n, f in {"a": "Exp12", "b": "Exp12/L0", "c": "Exp12x", "d": ""}.items()}
    if isinstance(expect, type):
        with pytest.raises(expect):
            pj.rename_folder(pid, "ann", old, new)
        return
    pj.rename_folder(pid, "ann", old, new)
    assert {n: pj.item(pid, "ann", i)["folder"] for n, i in ids.items()} == expect


def test_rename_folder_only_in_your_own_workspace(tmp_path):
    from plasmidbank.projects import Forbidden, Projects

    pj = Projects(tmp_path)
    pid = pj.personal("ann")
    pj.add_draft(pid, "ann", "a.gb", b"LOCUS x", "a", folder="Exp12")
    with pytest.raises(Forbidden):
        pj.rename_folder(pid, "bob", "Exp12", "Mine")
    shared = pj.create("ann", "NFR1")
    with pytest.raises(ValueError):
        pj.rename_folder(shared, "ann", "Exp12", "Mine")


def test_draft_opens_in_the_viewer_saves_in_place_and_bank_maps_copy_to_workspace(app):
    s, c = app
    gb = _gb(_rand(1200, 61), "linear1")
    r = c.post("/api/upload", files={"file": ("x.gb", gb), "plasmid_id": (None, ""), "construct_name": (None, "exon test"),
                                     "folder": (None, "RinRK1")}, headers=_as("ann")).json()
    key = r["key"]
    assert key < 0
    m = c.get(f"/api/maps/{key}", headers=_as("ann")).json()
    assert m["length"] == 1200 and m["draft"]["folder"] == "RinRK1" and m["file"] == "My workspace/RinRK1/exon test"
    assert c.get(f"/api/maps/{key}", headers=_as("bob")).status_code == 403  # members only
    assert c.get(f"/api/maps/{key}/orfs", headers=_as("bob")).status_code == 403
    assert c.get("/api/maps/-999999", headers=_as("ann")).status_code == 404
    assert c.get(f"/api/maps/{key}/orfs", headers=_as("ann")).status_code == 200  # the tabs work on a draft
    edits = {"add": [{"type": "exon", "label": "exon 1", "start": 10, "end": 90, "strand": 1}],
             "seq": [{"start": 5, "end": 4, "bases": "GAATTC"}]}
    r = c.post(f"/api/maps/{key}/save", json=edits, headers=_as("ann"))
    assert r.status_code == 200 and r.json()["in_place"], r.text
    m = c.get(f"/api/maps/{key}", headers=_as("ann")).json()
    assert m["length"] == 1206 and m["seq"][4:10] == "GAATTC"
    assert [(f["label"], f["start"], f["end"]) for f in m["features"]] == [("exon 1", 15, 96)]  # 0-based start
    assert not list(s.OUTBOX.rglob("*.gb"))  # in place, not to the outbox
    # the share sync pulled an edit in after the file was opened: a stale save must ask first
    import os, time as _t
    opened = c.get(f"/api/maps/{key}", headers=_as("ann")).json()["opened_at"]
    p = s.projects.draft_path(m["draft"]["project"], "ann", -key)
    os.utime(p, (opened + 5, opened + 5))
    stale = {"add": [{"type": "misc_feature", "label": "x", "start": 1, "end": 5}], "opened_at": opened}
    assert c.post(f"/api/maps/{key}/save", json=stale, headers=_as("ann")).status_code == 409
    assert c.post(f"/api/maps/{key}/save", json={**stale, "force": True}, headers=_as("ann")).status_code == 200
    r = c.post(f"/api/projects/{m['draft']['project']}/items/{-key}/move", json={"folder": "RinRK1/done", "label": "exon v2"},
               headers=_as("ann"))
    assert r.status_code == 200, r.text
    assert c.get(f"/api/maps/{key}", headers=_as("ann")).json()["file"] == "My workspace/RinRK1/done/exon v2"
    assert c.post(f"/api/projects/{m['draft']['project']}/items/{-key}/move", json={"folder": "x"},
                  headers=_as("bob")).status_code == 403
    # a bank map: Save a copy to My workspace, the original untouched
    k = _key(c, "EC00011")
    before = c.get(f"/api/maps/{k}", headers=_as("ann")).json()["seq"]
    r = c.post(f"/api/maps/{k}/copy", json={"edits": {"add": [{"type": "primer_bind", "label": "P1", "start": 1, "end": 20}]},
                                           "folder": "primers"}, headers=_as("ann")).json()
    cp = c.get(f"/api/maps/{r['key']}", headers=_as("ann")).json()
    assert cp["draft"]["folder"] == "primers" and cp["seq"] == before and any(f["label"] == "P1" for f in cp["features"])
    assert not any(f["label"] == "P1" for f in c.get(f"/api/maps/{k}", headers=_as("ann")).json()["features"])



def test_enzyme_finder_and_lab_list(app):
    _, c = app
    part, acc = _key(c, "EC00011"), _key(c, "PMC-00012")
    r = c.post("/api/enzymes/select", json={"conditions": [{"key": acc, "min": 1, "max": 1}, {"key": part, "min": 0, "max": 0}],
                                            "set": "commercial"}, headers=_as("ann")).json()
    assert r["maps"] == ["PMC-00012", "EC00011"] and r["searched"] > 500
    assert all(len(e["cuts"][0]) == 1 and e["cuts"][1] == [] for e in r["enzymes"])
    assert c.post("/api/enzymes/select", json={"conditions": [{"key": acc, "min": 1, "max": 1}], "set": "lab"},
                  headers=_as("ann")).status_code == 400  # empty freezer list says so
    assert c._do("PUT", "/api/lab/enzymes", _as("ann"), {"enzymes": "EcoRI, NotI, Bogus1"}).status_code == 400
    saved = c._do("PUT", "/api/lab/enzymes", _as("ann"), {"enzymes": "EcoRI, NotI, BsaI"}).json()
    assert saved["enzymes"] == ["BsaI", "EcoRI", "NotI"] and saved["updated_by"] == "ann"
    r = c.post("/api/enzymes/select", json={"conditions": [{"key": acc, "min": 0, "max": None}], "set": "lab"},
               headers=_as("bob")).json()
    assert {e["enzyme"] for e in r["enzymes"]} <= {"BsaI", "EcoRI", "NotI"} and r["searched"] == 3


def test_several_ab1_reads_get_support_counts_and_traces(app):
    from tests.abif_fixture import make_ab1

    _, c = app
    key = _key(c, "PMC-00012")
    ref = c.get(f"/api/maps/{key}", headers=_as("ann")).json()["seq"]
    sub = {"A": "C", "C": "G", "G": "T", "T": "A"}
    mut = ref[:400] + sub[ref[400]] + ref[401:]

    def ab1(seq):
        peaks = [10 * i + 5 for i in range(len(seq))]
        ch = {b: [50] * (10 * len(seq)) for b in "ACGT"}
        for i, base in enumerate(seq):
            ch[base][peaks[i]] = 1000
        return make_ab1(seq, [40] * len(seq), peaks, ch)

    rc = mut[300:700][::-1].translate(str.maketrans("ACGT", "TGCA"))  # a reverse read
    files = {"files": [("r1.ab1", ab1(mut[200:600])), ("r2.ab1", ab1(rc)), ("r3.ab1", ab1(ref[250:650]))]}
    rows = c.post(f"/api/maps/{key}/verify", files=files, headers=_as("ann")).json()
    assert [r["has_trace"] for r in rows] == [True, True, True]
    d1, d2 = rows[0]["discrepancies"], rows[1]["discrepancies"]
    assert [d["ref_pos"] for d in d1] == [401] and d1[0]["covered_by"] == 3 and d1[0]["seen_in"] == 2
    assert rows[1]["strand"] == "-" and [d["ref_pos"] for d in d2] == [401]
    for d in (d1[0], d2[0]):  # the trace is in map orientation: the focus peak is the changed base
        t = d["trace"]
        assert t["bases"][t["focus"]] == d["query_bases"]
        top = max("ACGT", key=lambda b: t["channels"][b][t["peaks"][t["focus"]]])
        assert top == d["query_bases"]
    assert rows[2]["discrepancies"] == []
    # the whole chromatogram of the reverse read, in map orientation (the viewer re-posts the file)
    t = c.post("/api/trace", files={"file": ("r2.ab1", ab1(rc)), "strand": (None, "-")}, headers=_as("ann")).json()
    assert t["bases"] == mut[300:700] and rows[1]["blocks"] == [[301, 1, 400]]
    assert c.post("/api/trace", files={"file": ("x.fasta", b">x\nACGT\n")}, headers=_as("ann")).status_code == 400


def test_upload_guardrails_count_and_size(app):
    s, c = app
    s.MAX_UPLOAD_MB = 1  # the real default is 10; keep the test payload small
    key = _key(c, "PMC-00012")
    many = {"files": [(f"r{i}.fasta", b">r\nACGT\n") for i in range(s.MAX_UPLOAD_FILES + 1)]}
    r = c.post(f"/api/maps/{key}/verify", files=many, headers=_as("ann"))
    assert r.status_code == 413 and "at most" in r.text
    big = b">r\n" + b"A" * (1024 * 1024 + 10)
    for url, files in ((f"/api/maps/{key}/verify", {"files": [("big.fasta", big)]}),
                       ("/api/trace", {"file": ("big.ab1", big)}),
                       ("/api/align/multiple", {"files": [("big.fasta", big)]})):
        r = c.post(url, files=files, headers=_as("ann"))
        assert r.status_code == 413 and "larger than 1 MB" in r.text, url


# ── admin panel ────────────────────────────────────────────────

def test_admin_routes_need_an_admin_and_guard_the_last_one(app):
    s, c = app
    assert c.get("/api/admin/users", headers=_as("ann")).status_code == 403
    assert c.get("/api/me", headers=_as("ann")).json()["role"] == "user"
    s.auth.save_users({"boss": {"role": "admin"}})
    assert c.get("/api/me", headers=_as("boss")).json()["role"] == "admin"
    assert c.get("/api/admin/sync", headers=_as("ann")).status_code == 403
    (s.STATE / "workspaces_status.json").write_text('{"finished": "2026-10-27T08:00:00", "conflicts": 1}')
    assert c.get("/api/admin/sync", headers=_as("boss")).json()["workspaces"]["conflicts"] == 1
    r = c.post("/api/admin/users", json={"name": "Carl", "role": "student"}, headers=_as("boss")).json()
    assert r == {"name": "carl", "role": "student", "password": ""}  # header mode: the proxy owns passwords
    assert c.post("/api/admin/users", json={"name": "carl"}, headers=_as("boss")).status_code == 400
    assert c.post("/api/admin/users", json={"name": "x y", "role": "user"}, headers=_as("boss")).status_code == 400
    assert c._do("PUT", "/api/admin/users/carl/role", _as("boss"), {"role": "user"}).status_code == 200
    assert c._do("PUT", "/api/admin/users/boss/role", _as("boss"), {"role": "user"}).status_code == 400  # last admin
    assert c._do("PUT", "/api/admin/users/carl/role", _as("boss"), {"role": "king"}).status_code == 400
    assert c.delete("/api/admin/users/boss", headers=_as("boss")).status_code == 400  # not yourself
    assert c.post("/api/admin/users/carl/reset", headers=_as("boss")).status_code == 400  # not in header mode
    assert c.delete("/api/admin/users/carl", headers=_as("boss")).status_code == 200
    users = c.get("/api/admin/users", headers=_as("boss")).json()["users"]
    assert [u["name"] for u in users] == ["boss"]
    audit = c.get("/api/admin/audit", headers=_as("boss")).json()
    assert [a["action"] for a in audit] == ["remove", "role", "add"] and audit[-1]["by"] == "boss"
    assert c.get("/admin", headers=_as("ann")).status_code == 200  # the page itself; its API is what is guarded


_EXPORT = ("Reference\tOligo_type\tOligo id \tOligo name\tOligo target\tOligo sequence\tComment\r\n"
           "AT\tDesalt\tPM-00001\tone\tx\tACGTACGTACGTACGTACGT\tc\t1/2/2020\r\n"
           "AT\tDesalt\tPM-00002\ttwo\tx\tGGGGCCCCAAAATTTTGGCC\tc\t1/2/2020\r\n")


def test_oligo_export_upload_previews_then_merges_without_removing(app):
    s, c = app
    s.auth.save_users({"boss": {"role": "admin"}})
    up = lambda text, apply: c.post("/api/admin/oligos", files={"file": ("export.tab", text), "apply": (None, apply)},
                                    headers=_as("boss"))
    assert c.post("/api/admin/oligos", files={"file": ("e.tab", _EXPORT)}, headers=_as("ann")).status_code == 403
    assert up("just text\n", "false").status_code == 400
    r = up(_EXPORT, "false").json()
    assert [n["id"] for n in r["new"]] == ["PM-00001", "PM-00002"] and not r["applied"] and r["oligos"] == 0
    r = up(_EXPORT, "true").json()
    assert r["applied"] and r["oligos"] == 2 and r["highest_pm"] == "PM-00002" and r["source"] == "uploaded"
    assert c.get("/api/stats", headers=_as("ann")).json()["oligos"] == 2
    # a partial export: one new row, one changed sequence; PM-00001 is absent but must stay
    part = ("Reference\tOligo_type\tOligo id \tOligo name\tOligo target\tOligo sequence\tComment\r\n"
            "AT\tDesalt\tPM-00002\ttwo\tx\tGGGGCCCCAAAATTTTGGCA\tc\t1/2/2020\r\n"
            "AT\tDesalt\tPM-00003\tthree\tx\tttgaagacaaTTTTCCCCGGGGAAAA\tc\t1/3/2020\r\n")
    r = up(part, "true").json()
    assert [n["id"] for n in r["new"]] == ["PM-00003"] and [x["id"] for x in r["changed"]] == ["PM-00002"]
    assert r["oligos"] == 3 and r["unchanged"] == 0
    assert up(part, "true").json()["unchanged"] == 2
    # the scratchpad's lookup: by sequence as ordered, case and whitespace ignored
    found = c.post("/api/oligos/lookup", json={"sequences": ["TTGAAGACAA TTTTCCCCGGGGAAAA", "ACGTACGTACGTACGTACGT", "AAAA"]},
                   headers=_as("ann")).json()
    assert found == {"TTGAAGACAA TTTTCCCCGGGGAAAA": ["PM-00003"], "ACGTACGTACGTACGTACGT": ["PM-00001"], "AAAA": []}
    # labels count: a plain primer never gets the PM number of its labelled twin, and vice versa
    lab = ("Reference\tOligo_type\tOligo id \tOligo name\tOligo target\tOligo sequence\tComment\r\n"
           "AT\tHPLC\tPM-00004\tone FAM\tx\t[6FAM]acgtACGTACGTACGTACGT\tc\t1/4/2020\r\n")
    up(lab, "true")
    found = c.post("/api/oligos/lookup", json={"sequences": ["ACGTACGTACGTACGTACGT", "[6fam] ACGTACGTACGTACGTACGT", "[Btn]ACGTACGTACGTACGTACGT"]},
                   headers=_as("ann")).json()
    assert list(found.values()) == [["PM-00001"], ["PM-00004"], []]
    assert (s.STATE / "oligo_export.prev.tab").exists()
    assert c.get("/api/admin/audit", headers=_as("boss")).json()[0]["action"] == "oligo export"


def test_local_accounts_one_time_passwords_and_cli(tmp_path, monkeypatch):
    for k, v in {"PLASMIDBANK_AUTH": "local", "PLASMIDBANK_USERS": tmp_path / "users.json",
                 "PLASMIDBANK_SECRET": "x" * 40}.items():
        monkeypatch.setenv(k, str(v))
    from plasmidbank.app import auth
    a = importlib.reload(auth)
    try:
        a.set_password("old", "an old password")  # an account from before roles: no role field
        assert a.role_of("old") == "user" and not a.must_change("old")
        pw = a.create_user("boss", "admin")
        assert len(pw) >= 12 and a.check_password("boss", pw) and a.must_change("boss")
        assert pw not in (tmp_path / "users.json").read_text()  # never stored in clear
        a.change_own_password("boss", pw, "my own password")
        assert not a.must_change("boss") and a.check_password("boss", "my own password")
        with pytest.raises(ValueError):
            a.change_own_password("boss", "wrong", "whatever123")
        new = a.reset_password("old")
        assert a.check_password("old", new) and not a.check_password("old", "an old password") and a.must_change("old")
        with pytest.raises(ValueError):
            a.remove_user("boss")  # last admin
        assert a.main(["role", "boss", "user"]) == 1  # the CLI keeps the guard too
        assert a.main(["role", "old", "admin"]) == 0 and a.main(["role", "boss", "user"]) == 0
        # `auth add` on the VM (2026-10-29): a short password is a message, not a traceback; nothing changes
        monkeypatch.setenv("PLASMIDBANK_NEW_PASSWORD", "short")
        assert a.main(["add", "old"]) == 1 and a.check_password("old", new)
        monkeypatch.delenv("PLASMIDBANK_NEW_PASSWORD")
        typed = iter(["a long password", "a different one"])  # typed twice, not the same
        monkeypatch.setattr(a.getpass, "getpass", lambda prompt="": next(typed))
        assert a.main(["add", "old"]) == 1 and a.check_password("old", new)
        # `auth reset NAME`: a one-time password to copy, changed at the next sign-in; creates the account if new
        monkeypatch.setattr(a, "one_time_password", lambda: "copy-me-0123")
        assert a.main(["reset", "old"]) == 0 and a.check_password("old", "copy-me-0123") and a.must_change("old")
        assert a.main(["reset", "newbie"]) == 0 and a.check_password("newbie", "copy-me-0123") and a.role_of("newbie") == "user"
        a.remove_user("newbie")
        a.remove_user("boss")
        assert [u["name"] for u in a.list_users()] == ["old"] and a.list_users()[0]["role"] == "admin"
        assert "role" in a.AUDIT_FILE.read_text()
    finally:
        monkeypatch.setenv("PLASMIDBANK_AUTH", "none")
        importlib.reload(auth)


def test_design_presets_carry_the_order_template_codes(app):
    _, c = app
    p = c.get("/api/design/presets", headers=_as("ann")).json()
    fam = next(m for m in p["oligo_mods"] if m["code"] == "[6FAM]")
    assert set(fam["positions"]) == {"5", "3"} and all(m["code"].startswith("[") for m in p["oligo_mods"])
    assert p["order"]["purification"] == ["Desalt", "Cartridge", "HPLC"] and p["order"]["tubes"][0] == 1


def test_sdm_endpoint(app):
    _, c = app
    key = _key(c, "PMC-00012")
    r = c.get(f"/api/maps/{key}/sdm?start=1000&end=999&new=GGATCC", headers=_as("ann")).json()
    assert r["kind"] == "insertion" and r["forward"]["tail"] == "GGATCC" and r["reverse"]["end"] == 999
    L = c.get(f"/api/maps/{key}", headers=_as("ann")).json()["length"]
    assert r["forward"]["start"] == 1000 and r["mutant_length"] == L + 6
    bad = c.get(f"/api/maps/{key}/sdm?feature=0&mutation=K48R", headers=_as("ann"))
    assert bad.status_code == 400 and "coding feature" in bad.json()["detail"]
    assert c.get(f"/api/maps/{key}/sdm?start=5&end=5&new=XYZ", headers=_as("ann")).status_code == 400


def test_project_list_says_which_item_links_a_map(app):
    _, c = app
    a = c.post("/api/projects", json={"name": "A"}, headers=_as("ann")).json()["id"]
    b = c.post("/api/projects", json={"name": "B"}, headers=_as("ann")).json()["id"]
    item = c.post(f"/api/projects/{a}/refs", json={"key": _key(c, "PMC-00012")}, headers=_as("ann")).json()["id"]
    rows = {p["id"]: p for p in c.get("/api/projects?plasmid_id=PMC-00012", headers=_as("ann")).json()}
    assert rows[a]["item_for"] == item and rows[b]["item_for"] is None
    assert "item_for" not in c.get("/api/projects", headers=_as("ann")).json()[0]
    assert c.get("/api/projects?plasmid_id=PMC-00012", headers=_as("bob")).json() == []  # not a member


def test_build_endpoints(app):
    _, c = app
    acc, part = _key(c, "PMC-00012"), _key(c, "EC00011")
    body = {"vector": acc, "open": {"start": 500, "end": 499}, "mode": "gibson",
            "inserts": [{"key": part, "start": 100, "end": 700, "strand": "+", "name": "part"},
                        {"seq": "CATCACCATCACCATCACCATCACCATCAC", "name": "tag"}]}
    r = c.post("/api/build", json=body, headers=_as("ann")).json()
    assert r["mode"] == "gibson" and len(r["junctions"]) == 3 and "LOCUS" in r["genbank"]
    assert {(p["piece"], p["dir"]) for p in r["primers"]} >= {("part", "F"), ("part", "R")}
    bad = c.post("/api/build", json={**body, "inserts": []}, headers=_as("ann"))
    assert bad.status_code == 400 and "1 to 5" in bad.json()["detail"]
    assert c.post("/api/build", json={"inserts": []}, headers=_as("ann")).status_code == 400
    m = c.post("/api/build/from_sequence", json={"text": ">p1\n" + "ACGT" * 20, "name": ""}, headers=_as("ann")).json()
    assert m["name"] == "p1" and m["length"] == 80


def test_delete_queued_upload(app):
    s, c = app
    assert not any(m["plasmid_id"] == "PMC-07777" for m in c.get("/api/maps", headers=_as("ann")).json())  # list cached
    gb = _gb(_rand(900, 31), "new").encode()
    files = {"file": ("x.gb", gb), "plasmid_id": (None, "PMC-07777"), "construct_name": (None, "mine")}
    r = c.post("/api/upload", files=files, headers=_as("ann")).json()
    assert any(m["plasmid_id"] == "PMC-07777" for m in c.get("/api/maps", headers=_as("ann")).json())
    assert c.delete("/api/uploads/PMC-07777", headers=_as("bob")).status_code == 403
    assert c.delete("/api/uploads/PMC-07777", headers=_as("ann")).status_code == 200
    assert not any(m["plasmid_id"] == "PMC-07777" for m in c.get("/api/maps", headers=_as("ann")).json())
    assert c.get("/api/uploads", headers=_as("ann")).json() == []
    files2 = {"file": ("y.gb", gb), "plasmid_id": (None, "PMC-07778"), "construct_name": (None, "x")}
    c.post("/api/upload", files=files2, headers=_as("ann"))
    rows = c.get("/api/uploads?all=true", headers=_as("bob")).json()
    assert [r["mine"] for r in rows] == [False]  # bob may not delete ann's upload
    assert c.delete("/api/uploads/PMC-07777", headers=_as("ann")).status_code == 400
    assert r["key"] is not None


def test_blast_endpoint(app):
    from plasmidbank.tools import blast
    if not blast.available():
        pytest.skip("BLAST+ not installed")
    _, c = app
    key = _key(c, "PMC-00012")
    seq = c.get(f"/api/maps/{key}", headers=_as("ann")).json()["seq"]
    r = c.post("/api/blast", json={"key": key, "start": 300, "end": 700}, headers=_as("ann")).json()
    top = r["hits"][0]
    assert top["key"] == key and top["self"] and (top["map_from"], top["map_to"]) == (300, 700)
    r = c.post("/api/blast", json={"seq": seq[1000:1300]}, headers=_as("ann")).json()
    assert r["hits"][0]["plasmid_id"] == "PMC-00012" and not r["hits"][0]["self"]
    assert c.post("/api/blast", json={"seq": "ACG"}, headers=_as("ann")).status_code == 400
    assert c.post("/api/blast", json={}, headers=_as("ann")).status_code == 400


def test_annotation_lists_and_annotate_on_upload(app):
    _, c = app
    key = _key(c, "PMC-00012")
    seq = c.get(f"/api/maps/{key}", headers=_as("ann")).json()["seq"]
    part = seq[300:360]
    # from a stretch of a map, into the common list
    r = c.post("/api/annotations", json={"label": "my part", "type": "promoter", "scope": "common",
                                         "key": key, "start": 301, "end": 360}, headers=_as("ann")).json()
    lst = c.get("/api/annotations", headers=_as("bob")).json()
    assert [(e["label"], e["scope"], e["seq"]) for e in lst["entries"]] == [("my part", "common", part)]
    props = c.post(f"/api/maps/{key}/annotate", headers=_as("bob")).json()["proposals"]
    assert [(p["start"], p["end"], p["strand"], p["label"]) for p in props] == [(301, 360, 1, "my part")]
    # bob renames the common entry; a private entry stays private
    assert c._do("PUT", f"/api/annotations/{r['id']}", _as("bob"), {"label": "Pmy"}).status_code == 200
    c.post("/api/annotations", json={"label": "bob only", "seq": seq[1000:1040], "scope": "mine"}, headers=_as("bob"))
    assert [e["label"] for e in c.get("/api/annotations", headers=_as("ann")).json()["entries"]] == ["Pmy"]
    # annotate on upload: off by default, then on
    assert c.get("/api/annotations", headers=_as("ann")).json()["annotate_on_upload"] is False
    c._do("PUT", "/api/annotations/settings", _as("ann"), {"annotate_on_upload": True})
    files = {"file": ("n.gb", _gb(_rand(200, 41) + part + _rand(300, 42), "new")),
             "plasmid_id": (None, "PMC-07801"), "construct_name": (None, "auto")}
    up_ = c.post("/api/upload", files=files, headers=_as("ann")).json()
    assert up_["annotated"] == 1
    feats = c.get(f"/api/maps/{up_['key']}", headers=_as("ann")).json()["features"]
    assert any(f["label"] == "Pmy" and f["start"] == 200 for f in feats)
    assert c.post("/api/annotations", json={"label": "x", "seq": "ACG"}, headers=_as("ann")).status_code == 400


def test_pcr_endpoint_takes_typed_primers_and_lab_oligo_ids(app):
    s, c = app
    k = _key(c, "EC00011")
    seq = c.get(f"/api/maps/{k}", headers=_as("ann")).json()["seq"]
    from plasmidbank.sequtil import revcomp
    body = {"primers": [{"name": "F", "seq": seq[50:72]}, {"name": "R", "seq": revcomp(seq[500:522])}], "polymerase": "Phusion"}
    r = c.post(f"/api/maps/{k}/pcr", json=body, headers=_as("ann"))
    assert r.status_code == 200, r.text
    r = r.json()
    assert r["product"] == 472 and r["genbank"].startswith("LOCUS") and r["polymerase"].startswith("Thermo Phusion")
    assert c.post(f"/api/maps/{k}/pcr", json={**body, "polymerase": "Pfu"}, headers=_as("ann")).status_code == 400
    assert c.post(f"/api/maps/{k}/pcr", json={"primers": [{"oligo_id": "PM-99999"}, {"oligo_id": "PM-1"}]},
                  headers=_as("ann")).status_code == 400


def test_overwrite_on_o_needs_the_id_keeps_the_old_file_and_refuses_a_changed_one(app):
    import json as _json
    from plasmidbank import uploads as up
    s, c = app
    k = _key(c, "EC00011")
    edits = _json.dumps({"add": [{"type": "misc_feature", "label": "new site", "start": 1, "end": 10}]})
    form = lambda cid: {"confirm_id": (None, cid), "edits": (None, edits)}
    assert c.post(f"/api/maps/{k}/overwrite", files=form("EC00012"), headers=_as("ann")).status_code == 400  # wrong ID
    r = c.post(f"/api/maps/{k}/overwrite", files=form("ec00011"), headers=_as("ann"))
    assert r.status_code == 200, r.text
    assert c.post(f"/api/maps/{k}/overwrite", files=form("EC00011"), headers=_as("ann")).status_code == 409  # one at a time
    m = c.get(f"/api/maps/{k}", headers=_as("ann")).json()  # the bank shows the new version at once
    assert any(f["label"] == "new site" for f in m["features"]) and any("waiting for the hourly sync" in x for x in m["messages"])
    bank = s.DATA
    (bank / "ENSA").mkdir()
    old = (bank / "EC00011_part.gb").read_bytes()
    assert up.apply_overwrites(s.STATE, bank)["done"] and (bank / "EC00011_part.gb").read_bytes() == old  # dry run
    up.apply_overwrites(s.STATE, bank, write=True)
    assert b"new site" in (bank / "EC00011_part.gb").read_bytes()
    assert [x.read_bytes() for x in (bank / "_replaced").iterdir()] == [old]
    st = c.get("/api/admin/sync", headers=_as("ann"))
    s.auth.save_users({"boss": {"role": "admin"}})
    assert c.get("/api/admin/sync", headers=_as("boss")).json()["overwrites"][0]["status"] == "done"
    # a second request, then someone edits the file on the share first: refused, nothing lost
    k2 = _key(c, "PMC-00012")
    assert c.post(f"/api/maps/{k2}/overwrite", files={"confirm_id": (None, "PMC-00012"), "edits": (None, edits)},
                  headers=_as("ann")).status_code == 200
    (bank / "PMC-00012_acceptor.gb").write_text((bank / "PMC-00012_acceptor.gb").read_text() + "\n")
    r = up.apply_overwrites(s.STATE, bank, write=True)
    assert r["refused"] and b"new site" not in (bank / "PMC-00012_acceptor.gb").read_bytes()
    # students may not overwrite
    s.auth.save_users({"boss": {"role": "admin"}, "stu": {"role": "student"}})
    assert c.post(f"/api/maps/{k2}/overwrite", files={"confirm_id": (None, "PMC-00012"), "edits": (None, edits)},
                  headers=_as("stu")).status_code == 403
    # the admin's fix button: a queued request is deleted, never written, and the bank shows the index's map again
    assert c.post(f"/api/maps/{k2}/overwrite", files={"confirm_id": (None, "PMC-00012"), "edits": (None, edits)},
                  headers=_as("ann")).status_code == 200
    queued = next(o for o in c.get("/api/admin/sync", headers=_as("boss")).json()["overwrites"] if o["status"] == "queued")
    assert c.delete(f"/api/admin/overwrites/{queued['id']}", headers=_as("ann")).status_code == 403  # admins only
    assert c.delete(f"/api/admin/overwrites/{queued['id']}", headers=_as("boss")).json()["was"] == "queued"
    assert not list(s.OVERWRITES.glob(f"{queued['id']}.*"))
    assert not any(f["label"] == "new site" for f in c.get(f"/api/maps/{_key(c, 'PMC-00012')}", headers=_as("ann")).json()["features"])
    assert up.apply_overwrites(s.STATE, bank, write=True)["done"] == []


def test_overwrite_of_a_clc_becomes_gbk_and_paths_stay_inside_the_bank(tmp_path):
    import hashlib, json as _json
    from pathlib import Path
    from plasmidbank import uploads as up
    state, bank = tmp_path / "state", tmp_path / "bank"
    (bank / "ENSA").mkdir(parents=True)
    (bank / "PMC1001-2000").mkdir()
    clc, twin = bank / "PMC1001-2000" / "PMC-01118.clc", bank / "PMC1001-2000" / "PMC-01118.gbk"
    clc.write_bytes(b"CLC\x00binary")
    twin.write_text("old export")
    q = state / "overwrites"
    q.mkdir(parents=True)
    (q / "a.gb").write_text("LOCUS new")
    reqs = {"a": "PMC1001-2000/PMC-01118.clc", "b": "../outside.gb", "c": "Workspaces/ann/x.gb"}
    for rid, rel in reqs.items():
        (q / f"{rid}.json").write_text(_json.dumps({"id": rid, "rel": rel, "new_name": Path(rel).stem + ".gbk", "file": "a.gb",
                                                   "user": "ann", "status": "queued",
                                                   "sha_before": hashlib.sha256(clc.read_bytes()).hexdigest()}))
    r = up.apply_overwrites(state, bank, write=True)
    assert len(r["done"]) == 1 and len(r["refused"]) == 2
    assert not clc.exists() and twin.read_text() == "LOCUS new"
    assert sorted(x.name.split("_", 2)[-1] for x in (twin.parent / "_replaced").iterdir()) == ["PMC-01118.clc", "PMC-01118.gbk"]
