"""Fixes from the 2026-10-08 security audit: server-side must_change, login throttle, Origin check,
the placeholder secret, bounded Golden Gate, and the heavy-job slots."""

import importlib
import socket
import threading
import time
import urllib.parse

import pytest
import uvicorn

from tests.test_bench import _rand
from tests.test_projects import Client, _gb


@pytest.fixture
def local_app(tmp_path, monkeypatch):
    data = tmp_path / "bank"
    data.mkdir()
    (data / "PMC-00012_x.gb").write_text(_gb(_rand(2000, 31), "x"))
    for k, v in {
        "PLASMIDBANK_DATA": data,
        "PLASMIDBANK_STATE": tmp_path / "state",
        "PLASMIDBANK_UPLOADS": tmp_path / "uploads",
        "PLASMIDBANK_OUTBOX": tmp_path / "outbox",
        "PLASMIDBANK_OLIGOS": tmp_path / "none.tab",
        "PLASMIDBANK_AUTH": "local",
        "PLASMIDBANK_SECRET": "s" * 40,
        "PLASMIDBANK_SECURE_COOKIE": "0",
        "PLASMIDBANK_USERS": tmp_path / "state" / "users.json",
    }.items():
        monkeypatch.setenv(k, str(v))
    monkeypatch.delenv("PLASMIDBANK_INDEX", raising=False)
    from plasmidbank.app import auth, server

    a = importlib.reload(auth)
    s = importlib.reload(server)
    sock = socket.socket()
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()
    srv = uvicorn.Server(
        uvicorn.Config(s.app, host="127.0.0.1", port=port, log_level="warning")
    )
    th = threading.Thread(target=srv.run, daemon=True)
    th.start()
    while not srv.started:
        time.sleep(0.02)
    yield s, a, Client(f"http://127.0.0.1:{port}"), port
    srv.should_exit = True
    th.join(5)
    monkeypatch.setenv("PLASMIDBANK_AUTH", "none")
    importlib.reload(auth)


def _login(c, user, pw):
    """POST the login form without following the redirect; returns (status, cookie header)."""
    import http.client

    host, port = c.base.split("//")[1].split(":")
    conn = http.client.HTTPConnection(host, int(port), timeout=30)
    body = urllib.parse.urlencode({"user": user, "password": pw})
    conn.request(
        "POST", "/login", body, {"Content-Type": "application/x-www-form-urlencoded"}
    )
    r = conn.getresponse()
    r.read()
    cookie = (r.getheader("set-cookie") or "").split(";")[0]
    return r.status, cookie


def test_one_time_password_only_reaches_the_password_change(local_app):
    s, a, c, _ = local_app
    pw = a.create_user("eve.spaced")
    status, cookie = _login(c, "eve.spaced", pw)
    assert status == 303 and cookie
    h = {"Cookie": cookie}
    assert c.get("/api/me", headers=h).json()["must_change"] is True
    r = c.get("/api/maps", headers=h)
    assert (
        r.status_code == 403 and "password" in r.text
    )  # the API refuses until the password is changed
    assert (
        c.post(
            "/api/account/password",
            json={"current": pw, "new": "emma's own password"},
            headers=h,
        ).status_code
        == 200
    )
    assert c.get("/api/maps", headers=h).status_code == 200


def test_login_throttle_per_user_name(local_app):
    s, a, c, _ = local_app
    a.set_password("jane.doe", "the right password")
    for _ in range(s._LOGIN_MAX_FAILS):
        assert _login(c, "jane.doe", "wrong password")[0] == 401
    assert (
        _login(c, "jane.doe", "the right password")[0] == 429
    )  # locked, even with the right one
    s._LOGIN_FAILS.clear()
    assert _login(c, "jane.doe", "the right password")[0] == 303


def test_cross_site_writes_refused(local_app):
    s, a, c, port = local_app
    a.set_password("carol", "carols long password")
    _, cookie = _login(c, "carol", "carols long password")
    h = {"Cookie": cookie, "Origin": "https://evil.example.org"}
    assert c.post("/api/blast", json={"seq": "ACGT"}, headers=h).status_code == 403
    same = {"Cookie": cookie, "Origin": f"http://127.0.0.1:{port}"}
    assert c.post("/api/blast", json={"seq": "ACGT"}, headers=same).status_code != 403


def test_placeholder_secret_refused(monkeypatch):
    monkeypatch.setenv("PLASMIDBANK_AUTH", "local")
    monkeypatch.setenv(
        "PLASMIDBANK_SECRET", "change-me-to-a-long-random-string-of-at-least-32-chars"
    )
    from plasmidbank.app import auth

    a = importlib.reload(auth)
    try:
        with pytest.raises(SystemExit):
            a.check_config()
    finally:
        monkeypatch.setenv("PLASMIDBANK_AUTH", "none")
        importlib.reload(auth)


def test_ligation_search_is_bounded(monkeypatch):
    from plasmidbank.tools import cloning

    monkeypatch.setattr(cloning, "MAX_LIGATION_STEPS", 50)
    end = cloning.End("GGAG", "5'")  # identical ends join in every order: n! * 2^n paths
    frags = [cloning.Fragment(f"p{i}", "A" * 100, end, end, 0, 100) for i in range(8)]
    with pytest.raises(ValueError, match="too many possible assemblies"):
        cloning.ligate(frags)


def test_busy_heavy_slots_answer_429(local_app):
    s, a, c, _ = local_app
    a.set_password("dan.lee", "dans long password!")
    _, cookie = _login(c, "dan.lee", "dans long password!")
    held = [s._HEAVY.acquire(blocking=False) for _ in range(2)]
    try:
        assert all(held)
        r = c.post(
            "/api/blast",
            json={"seq": "ACGTACGTACGTACGTACGTACGT"},
            headers={"Cookie": cookie},
        )
        assert r.status_code in (
            429,
            503,
        )  # 503 when BLAST+ is not installed on this machine
    finally:
        for _ in held:
            s._HEAVY.release()
