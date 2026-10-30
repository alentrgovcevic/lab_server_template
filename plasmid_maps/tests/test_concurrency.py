"""Several people at once: the usage log, and the guards against losing someone's write."""

import json
import threading

from tests.test_bench import _rand
from tests.test_projects import _as, _gb, _key, app  # noqa: F401  (app is a fixture)


def test_usage_log_rows(app):
    s, c = app
    key = _key(c, "PMC-00012")
    c.get(f"/api/maps/{key}", headers=_as("ann"))
    c.get("/api/maps", headers=_as("bob"))
    c.get(
        "/auth/check",
        headers={
            **_as("bob"),
            "x-forwarded-uri": "/genomes/api/genomes/3/region?chr=1",
        },
    )
    assert c.get("/auth/check", headers={**_as("bob"), "x-forwarded-uri": "/static/app.js?v=2"}).status_code == 200
    c.get("/healthz")
    rows = [
        json.loads(x)
        for f in (s.USAGE).glob("usage-*.jsonl")
        for x in f.read_text().splitlines()
    ]
    routes = {(r["user"], r["route"]) for r in rows}
    assert ("ann", "/api/maps/{key}") in routes  # the template, not the map number
    assert ("bob", "genomes:/genomes/api/genomes/{n}/region") in routes
    assert not any(r["route"] == "/healthz" or r["route"].endswith(".js") for r in rows)
    assert all(r["inflight"] >= 1 and r["ms"] >= 0 for r in rows)


def test_resource_sample_reads_host_figures(tmp_path):
    from plasmidbank.app.usage import sample

    proc = tmp_path / "proc"
    proc.mkdir()
    (proc / "meminfo").write_text("MemTotal:  4009888 kB\nMemFree: 1 kB\nMemAvailable:  2004944 kB\n")
    (proc / "loadavg").write_text("0.50 0.40 0.30 1/399 1\n")
    (tmp_path / "memory.current").write_text(str(300 * 1048576))
    (proc / "stat").write_text("cpu  100 0 100 700 100 0 0 0 0 0\n")
    _, prev = sample(tmp_path, None, proc, tmp_path)
    (proc / "stat").write_text("cpu  250 0 250 800 100 0 0 0 0 0\n")  # 300 busy of 400 more ticks
    row, _ = sample(tmp_path, prev, proc, tmp_path)
    assert (row["cpu"], row["mem_total"], row["mem_avail"], row["app_mem"], row["load1"]) == (75.0, 3915, 1957, 300, 0.5)
    assert 0 < row["disk_used"] <= row["disk_total"]


def _report():
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[2] / "deploy" / "usage_report.py"
    spec = importlib.util.spec_from_file_location("usage_report", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_usage_report_sessions_tools_and_verdicts():
    r = _report()
    t0 = 1_800_000_000
    rows = [{"t": t0 + dt, "user": u, "route": route, "status": st, "ms": ms, "inflight": 1}
            for dt, u, route, st, ms in [
                (0, "ann", "/api/maps", 200, 300), (600, "ann", "/api/blast", 200, 400),
                (601, "ann", "/api/blast", 429, 1), (3 * 3600, "ann", "/api/maps/{key}/primers", 200, 20),
                (60, "bob", "genomes:/api/genomes/{n}/region", 200, 2)]]
    metrics = [{"t": t0 + i * 86400, "cpu": 90.0, "load1": 2.0, "cpus": 2, "mem_total": 4000, "mem_avail": 400,
                "app_mem": 500, "disk_total": 100.0, "disk_used": 10.0 + i} for i in range(3)]
    a = r.analyse(rows, metrics)
    ann = next(u for u in a["users"] if u["user"] == "ann")
    assert (ann["sessions"], ann["hours"]) == (2, round((11 + 1) / 60, 1))  # 0-601 s is one session, 3 h later another
    assert {t["tool"] for t in a["tools"]} == {"Browsing the bank", "BLAST", "Primers", "Genome browser"}
    assert (a["heavy"], a["busy"], a["peak_users"]) == (2, 1, 2)
    v = {topic: verdict for topic, verdict, _ in r.verdicts(a)}
    assert v == {"Cores": "more would help", "RAM": "add some", "Disk": "plan more", "Bottleneck": "BLAST"}
    assert a["resources"]["days_to_75"] == 63  # 12 GiB used, 1 GiB a day, 75 GiB at 75 %


def test_one_new_id_uploaded_by_many_at_once_is_stored_once(app):
    s, c = app
    gb = _gb(_rand(900, 41), "race").encode()
    codes = []

    def upload(i):
        files = {
            "file": ("r.gb", gb),
            "plasmid_id": (None, "PMC-07001"),
            "construct_name": (None, f"try{i}"),
        }
        codes.append(
            c.post("/api/upload", files=files, headers=_as(f"user{i}")).status_code
        )

    threads = [threading.Thread(target=upload, args=(i,)) for i in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert sorted(codes) == [200] + [400] * 7
    assert len(list((s.UPLOADS).glob("*/PMC-07001*"))) == 1
    assert (
        sum(
            m["plasmid_id"] == "PMC-07001"
            for m in c.get("/api/maps", headers=_as("ann")).json()
        )
        == 1
    )


def test_lab_enzyme_list_refuses_a_stale_version(app):
    _, c = app
    put = lambda user, body: c._do("PUT", "/api/lab/enzymes", _as(user), body)  # noqa: E731
    seen = c.get("/api/lab/enzymes", headers=_as("ann")).json()[
        "version"
    ]  # ann and bob both open v0
    assert (
        put("ann", {"enzymes": "EcoRI", "version": seen}).json()["version"] == seen + 1
    )
    r = put("bob", {"enzymes": "NotI", "version": seen})
    assert r.status_code == 409 and "ann" in r.json()["detail"]
    assert c.get("/api/lab/enzymes", headers=_as("bob")).json()["enzymes"] == [
        "EcoRI"
    ]  # ann's change kept
    assert (
        put("bob", {"enzymes": "EcoRI, NotI", "version": seen + 1}).status_code == 200
    )
    assert (
        put("bob", {"enzymes": "EcoRI"}).status_code == 200
    )  # a page from before the guard still saves


def test_save_as_new_version_flags_newer_versions_and_moved_keys(app):
    s, c = app
    key = _key(c, "PMC-00012")
    m_ann = c.get(f"/api/maps/{key}", headers=_as("ann")).json()
    m_bob = c.get(f"/api/maps/{key}", headers=_as("bob")).json()
    edit = {
        "add": [
            {"label": "x", "type": "misc_feature", "start": 10, "end": 40, "strand": 1}
        ]
    }
    save = lambda user, m, **kw: c.post(
        f"/api/maps/{key}/save",
        headers=_as(user),  # noqa: E731
        json={**edit, "expect_file": m["file"], "opened_at": m["opened_at"], **kw},
    )

    assert save("bob", m_bob).status_code == 200
    r = save("ann", m_ann)
    assert r.status_code == 409 and [n["user"] for n in r.json()["newer"]] == ["bob"]
    first = save("ann", m_ann, force=True).json()
    second = save(
        "ann", {**m_ann, "opened_at": first["saved_at"]}
    )  # nothing newer since her own save
    assert (
        second.status_code == 200 and second.json()["saved_to"] != first["saved_to"]
    )  # same second, two files
    assert len(list((s.OUTBOX / "ann").glob("*.gb"))) == 2

    r = save(
        "ann", {**m_ann, "file": "some/other_map.gb"}, force=True
    )  # the key now points elsewhere
    assert r.status_code == 409 and len(list((s.OUTBOX / "ann").glob("*.gb"))) == 2
