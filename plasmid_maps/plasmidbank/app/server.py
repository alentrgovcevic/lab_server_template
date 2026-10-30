"""Local web app: browse the map bank, verify clones, match primers.

Run (from plasmid_maps/):  ../.venv/bin/python -m plasmidbank.app.server        (http://127.0.0.1:8765)
Env:  PLASMIDBANK_DATA    folder of maps   (default: examples/bank/)
      PLASMIDBANK_OLIGOS  oligo export     (default: examples/oligos.tab)
"""

from __future__ import annotations

import hmac
import json
import os
import re
import threading
import time
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import asdict
from pathlib import Path

from fastapi import Body, FastAPI, File, Form, HTTPException, Request, UploadFile
from starlette.concurrency import run_in_threadpool
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, PlainTextResponse, RedirectResponse, Response

from . import auth, modules, usage
from .activity import Activity

from ..bank import Bank
from ..formats import read_map
from .. import uploads as up
from .. import labels
from ..projects import Forbidden, Projects
from ..oligos import load_oligos
from ..tools.align import load_queries, read_support, trace_full, trace_window, verify
from ..tools.primers import group_oligos, match_map
from ..tools import thermo
from ..tools.design import FUSION_SITES, GG_PREFIX, design_cloning, design_pcr, design_sequencing, design_walk
from ..tools.analysis import find_orfs, mafft_available, multiple_alignment, translate_region
from ..tools.cloning import LADDERS, digest, flip_product, fusion_frames, gel_lane, golden_gate, product_genbank
from ..formats import MapRecord
from ..sequtil import revcomp
from ..tools.edit import apply_edits
from ..tools.clean import ALL_RULES, DEFAULT_ON, OPT_IN, clean_record, read_genbank
from ..tools.sites import COMMON, SETS, commercial, enzymes_for, find_sites, select_enzymes

REPO = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("PLASMIDBANK_DATA", REPO.parent / "examples" / "bank"))
INDEX_DB = os.environ.get("PLASMIDBANK_INDEX", "")  # v1 index: load maps from it instead of DATA
INDEX_ROOT = os.environ.get("PLASMIDBANK_ROOT", "")  # the scanned root, for display paths
OLIGOS = Path(os.environ.get("PLASMIDBANK_OLIGOS", REPO.parent / "examples" / "oligos.tab"))
UPLOADS = Path(os.environ.get("PLASMIDBANK_UPLOADS", Path.home() / "PlasmidBank" / "uploads"))
STATIC = Path(__file__).parent / "static"
STATE = Path(os.environ.get("PLASMIDBANK_STATE", Path.home() / "PlasmidBank"))  # projects, feedback
GENOMES_URL = os.environ.get("PLASMIDBANK_GENOMES_URL", "http://127.0.0.1:8766")
STATE.mkdir(parents=True, exist_ok=True)
projects = Projects(STATE)
from ..annotations import Library as _Library  # noqa: E402

annotations = _Library(STATE)  # common + private annotation lists

auth.check_config()
app = FastAPI(title="Plasmid map bank")
if INDEX_DB:
    DATA = Path(INDEX_ROOT) if INDEX_ROOT else DATA


def _display(m) -> str:
    p = Path(m.display_path or m.path)
    try:
        return str(p.relative_to(DATA))
    except ValueError:
        return str(p)


def _load_bank() -> Bank:
    """The bank as the index (or folder) holds it, plus queued uploads not yet on the share. Used at start and
    by the in-place reload after a scan (2026-10-27: daily scan, no restart)."""
    global _on_o
    if INDEX_DB:
        # no index yet (first boot on a server): start empty rather than parse the share over SMB
        b = Bank.from_index(INDEX_DB, INDEX_ROOT or None) if Path(INDEX_DB).exists() else Bank(DATA)
    else:
        b = Bank.load(DATA)
    _on_o = set(b.by_id())  # IDs from the share; uploads of these were synced and indexed already
    ledger = up._ledger(UPLOADS)
    for p in up.queued_ids(UPLOADS).values():
        if up.identify(p.name).plasmid_id not in _on_o and not _deleted_from_o(p, ledger):
            b.maps.append(_read_upload(p))
    _overlay_overwrites(b)
    return b


OVERWRITES = STATE / "overwrites"


def _overwrite_requests() -> list[dict]:
    out = []
    for f in sorted(OVERWRITES.glob("*.json")) if OVERWRITES.is_dir() else []:
        try:
            out.append(json.loads(f.read_text()))
        except (OSError, ValueError):
            pass
    return out


def _overlay_overwrites(b: Bank) -> None:
    """Show an overwrite the moment it is asked for: the new version replaces the map in the list until
    the scan has seen the file on the share (queued, or done after the index was last written)."""
    stamp = max((q.stat().st_mtime for q in (Path(INDEX_DB), Path(INDEX_DB + "-wal")) if q.exists()), default=0) if INDEX_DB else None
    from datetime import datetime
    for r in _overwrite_requests():
        if r.get("status") == "refused" or (r.get("status") == "done" and (stamp is None or datetime.fromisoformat(r["finished"]).timestamp() < stamp)):
            continue
        new_rel = str(Path(r["rel"]).with_name(r["new_name"]))
        try:
            m = read_map(OVERWRITES / r["file"])
        except Exception:
            continue
        m.display_path = str(DATA / new_rel)
        ids = up.identify(Path(new_rel).name, m.name)  # the ID comes from the real file name, as in a scan
        m.plasmid_id, m.kind, m.id_source, m.xrefs = ids.plasmid_id, ids.kind, ids.source, ids.xrefs
        m.messages = [*m.messages, f"overwritten on the share by {r['user']} {r['requested'][:16].replace('T', ' ')}"
                      + ("" if r.get("status") == "done" else " — waiting for the hourly sync")]
        old = [i for i, x in enumerate(b.maps) if _display(x) in (r["rel"], new_rel)]
        for i in reversed(old[1:]):
            b.maps.pop(i)
        if old:
            b.maps[old[0]] = m


def _read_upload(p: Path):
    m = read_map(p)
    m.display_path = f"uploads/{p.parent.name}/{p.name}"
    return m


def _deleted_from_o(p: Path, ledger: dict) -> bool:
    """Synced, and the copy on the share is gone: someone deleted it there on purpose, so the bank
    drops it too. The upload stays in the state dir. The ledger's dest is the sync container's path;
    its last two parts (bin folder, file) are the same under this app's read-only root."""
    e = ledger["files"].get(f"{p.parent.name}/{p.name}", {})
    if e.get("status") != "synced" or not e.get("dest"):
        return False
    d = Path(e["dest"])
    return not (DATA / d.parent.name / d.name).exists()


_on_o: set = set()
bank = _load_bank()
def _index_stamp() -> float | None:
    """Newest change to the index, its write-ahead log included."""
    ts = [q.stat().st_mtime for q in (Path(INDEX_DB), Path(INDEX_DB + "-wal")) if q.exists()] if INDEX_DB else []
    return max(ts) if ts else None


_index_mtime = _index_stamp()
_index_seen = _index_mtime


def reload_if_index_changed() -> bool:
    """After the daily scan rewrote the index: rebuild the bank and swap it in place, so nobody waits
    for a restart. Only once the index has stopped changing (two checks in a row see the same stamp),
    so a scan in progress is never read half-way. The old list keeps serving until the new one is
    complete. True when it reloaded."""
    global _index_mtime, _index_seen
    t = _index_stamp()
    if t is None or t == _index_mtime:
        return False
    if t != _index_seen:  # still being written, or just finished: look again next time
        _index_seen = t
        return False
    fresh = _load_bank()
    with _BANK_LOCK:
        bank.maps[:] = fresh.maps
        bank.skipped[:] = fresh.skipped
        _bank_changed()
        _hits_cache.clear()
    _index_mtime = t
    print(f"bank reloaded from the index: {len(bank.maps):,} maps")
    _warm_blast()
    return True


def _watch_index() -> None:
    import time as _t
    while True:
        _t.sleep(60)
        try:
            reload_if_index_changed()
        except Exception as e:  # a half-written index: try again next minute
            print(f"bank reload failed: {e}")
# The admin's uploaded export (state dir) wins over PLASMIDBANK_OLIGOS, which is only the seed
# (on Docker a read-only mount). See POST /api/admin/oligos.
OLIGOS_STATE = STATE / "oligo_export.tab"


def _oligo_source() -> Path:
    return OLIGOS_STATE if OLIGOS_STATE.exists() else OLIGOS


def _load_oligos() -> None:
    """(Re)load the oligo collection and everything derived from it."""
    global oligos, oligo_groups, oligo_skipped, _oligo_names
    src = _oligo_source()
    oligos = load_oligos(src) if src.exists() else []
    oligo_groups, oligo_skipped = group_oligos(oligos)
    _oligo_names = {o.name for o in oligos} | {o.oligo_id for o in oligos}
    _hits_cache.clear()


_hits_cache: dict[int, list] = {}
_load_oligos()
_dups = {i: g for g in bank.duplicate_groups() for i in g}

_OPEN = {"/healthz", "/login", "/logout", "/auth/check"}
# a user still on a one-time password may only reach the dashboard and change it (security audit 2026-10-08)
_MUST_CHANGE_OK = {"/", "/api/me", "/api/account/password", "/api/links"}
# A lab database (e.g. FileMaker's Insert from URL) fetches label images and cannot sign in, so these two routes (and only
# these) also take `Authorization: Bearer <PLASMIDBANK_LABEL_TOKEN>`. No token configured
# = refused (503). Without the header a signed-in browser session works as for any other route.
_LABEL_ROUTE = re.compile(r"^/api/(qr|code128)/[^/]+\.png$")
LABEL_USER = "label-token"


def _label_token_gate(request: Request) -> JSONResponse | None:
    """None when the bearer token is right; otherwise the refusal to send."""
    token = os.environ.get("PLASMIDBANK_LABEL_TOKEN", "")
    if not token:
        return JSONResponse({"detail": "label token not configured"}, status_code=503)
    scheme, _, given = request.headers.get("authorization", "").partition(" ")
    if scheme.lower() != "bearer" or not hmac.compare_digest(given.strip().encode(), token.encode()):
        return JSONResponse({"detail": "wrong label token"}, status_code=401,
                            headers={"WWW-Authenticate": "Bearer"})
    return None


def _cross_site(request: Request) -> bool:
    """A browser write from another site. Requests without an Origin header (curl, tests) pass."""
    origin = request.headers.get("origin")
    if request.method in ("GET", "HEAD", "OPTIONS") or not origin:
        return False
    return origin.split("://", 1)[-1].rstrip("/") != request.headers.get("host", "")


@app.middleware("http")
async def _auth_gate(request: Request, call_next):
    if _cross_site(request):
        return JSONResponse({"detail": "cross-site request refused"}, status_code=403)
    if request.url.path in _OPEN:
        return await call_next(request)
    if "authorization" in request.headers and request.method in ("GET", "HEAD") and _LABEL_ROUTE.match(request.url.path):
        refused = _label_token_gate(request)
        if refused:
            return refused
        request.state.user = LABEL_USER
        return await call_next(request)
    user = auth.current_user(request.headers, request.cookies)
    if not user:
        if request.url.path.startswith("/api/"):
            return JSONResponse({"detail": "not signed in"}, status_code=401)
        return RedirectResponse("/login", status_code=303)
    if request.url.path not in _MUST_CHANGE_OK and auth.must_change(user):
        if request.url.path.startswith("/api/"):
            return JSONResponse({"detail": "choose your own password first"}, status_code=403)
        return RedirectResponse("/#password", status_code=303)
    request.state.user = user
    _USER.set(user)
    activity.touch(user)
    return await call_next(request)


activity = Activity(STATE)  # last active per user, for /admin
modules.setup(STATE)
app.include_router(modules.router)  # /ui.js, the Basic view settings and your own modules (app/modules.py)

USAGE = STATE / "usage"
_inflight = 0  # requests in progress; only touched on the event loop
usage.start_sampler(USAGE, STATE, lambda: _inflight)  # one resource row a minute (usage.py)


@app.middleware("http")
async def _usage_log(request: Request, call_next):
    """One row per request in <state>/usage/ (usage.py). Outermost layer, so it times the whole request."""
    global _inflight
    if request.url.path in usage.SKIP:
        return await call_next(request)
    _inflight += 1
    row = {"t": round(time.time(), 3), "m": request.method, "inflight": _inflight}
    t0, status = time.perf_counter(), 500
    try:
        resp = await call_next(request)
        status = resp.status_code
        return resp
    finally:
        _inflight -= 1
        fwd = request.headers.get("x-forwarded-uri", "").split("?", 1)[0]
        if request.url.path == "/auth/check":  # a genome-browser request, asked about by Caddy
            route = "genomes:" + usage.folded(fwd)
        else:
            route = getattr(request.scope.get("route"), "path", None) or usage.folded(request.url.path)
        row.update(user=getattr(request.state, "user", None), route=route, status=status,
                   ms=round((time.perf_counter() - t0) * 1000, 1))
        # its scripts, styles and icons are skipped: one genome-browser page view would be many rows
        if not (request.url.path == "/auth/check" and fwd.endswith(usage.STATIC_SUFFIXES)):
            try:
                usage.record(USAGE, row)
            except OSError:
                pass  # a full disk must not take the site down


# Login: per user name, a few failures then a pause; and few scrypt checks at once (16 MB each)
_LOGIN_FAILS: dict[str, list[float]] = {}
_LOGIN_MAX_FAILS, _LOGIN_WINDOW = 5, 15 * 60
_SCRYPT_SLOTS = threading.BoundedSemaphore(4)


def _login_blocked(name: str) -> bool:
    now = time.time()
    recent = [t for t in _LOGIN_FAILS.get(name, []) if now - t < _LOGIN_WINDOW]
    _LOGIN_FAILS[name] = recent
    return len(recent) >= _LOGIN_MAX_FAILS


# Heavy jobs (BLAST, MAFFT, read alignment, Golden Gate) share two slots on the 2-CPU VM
_HEAVY = threading.BoundedSemaphore(int(os.environ.get("PLASMIDBANK_HEAVY_JOBS", 2)))
MAX_QUERY_BP = 250_000  # one BLAST query or one map/read in an alignment
MAX_ALIGN_TOTAL_BP = 1_000_000  # all sequences of one multiple alignment together
MAX_READ_BP = 50_000  # one read in clone verification
MAX_GG_INPUTS = 10


@contextmanager
def _heavy_slot():
    if not _HEAVY.acquire(blocking=False):
        raise HTTPException(429, "the server is busy with other analyses; try again in a minute")
    try:
        yield
    finally:
        _HEAVY.release()


@app.get("/healthz")
def healthz():
    return {"ok": True, "maps": len(bank.maps)}


_LOGIN_PAGE = """<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><meta name="color-scheme" content="light dark">
<title>Lab server: sign in</title>
<style>:root{--bg:#f7f7f5;--card:#fff;--ink:#1d2327;--muted:#5f6b73;--line:#dde1e4;--accent:#2b7a78;--bad:#c0392b}
@media (prefers-color-scheme: dark){:root{--bg:#15191c;--card:#1d2327;--ink:#e6e9eb;--muted:#9aa5ad;--line:#333b41;--accent:#5fb3b0;--bad:#f07b6e}}
body{margin:0;min-height:100vh;display:grid;place-items:center;font:15px/1.5 system-ui,sans-serif;background:var(--bg);color:var(--ink)}
form{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:26px 28px;width:min(320px,90vw)}
h1{font-size:18px;margin:0 0 14px}label{display:block;font-size:13px;color:var(--muted);margin-top:10px}
input{width:100%;box-sizing:border-box;padding:8px 10px;margin-top:4px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--ink);font:inherit}
button{margin-top:16px;width:100%;padding:9px;border:0;border-radius:6px;background:var(--accent);color:#fff;font:inherit;cursor:pointer}
.err{color:var(--bad);font-size:13px;margin-top:10px}</style></head><body>
<form method="post" action="/login"><h1>Lab server</h1>
<label>User name<input name="user" autocomplete="username" required autofocus></label>
<label>Password<input name="password" type="password" autocomplete="current-password" required></label>
<button>Sign in</button>__ERR__</form></body></html>"""


@app.get("/login")
def login_page():
    if auth.MODE != "local":
        return RedirectResponse("/", status_code=303)
    return HTMLResponse(_LOGIN_PAGE.replace("__ERR__", ""))


@app.post("/login")
def login(user: str = Form(...), password: str = Form(...)):
    if auth.MODE != "local":
        return RedirectResponse("/", status_code=303)
    name = user.strip().lower()
    if _login_blocked(name):
        return HTMLResponse(_LOGIN_PAGE.replace(
            "__ERR__", '<p class="err">Too many failed attempts. Wait 15 minutes, or ask an admin.</p>'), status_code=429)
    if not _SCRYPT_SLOTS.acquire(timeout=10):
        return HTMLResponse(_LOGIN_PAGE.replace("__ERR__", '<p class="err">Busy, try again.</p>'), status_code=429)
    try:
        ok = auth.check_password(name, password)
    finally:
        _SCRYPT_SLOTS.release()
    if not ok:
        _LOGIN_FAILS.setdefault(name, []).append(time.time())
        return HTMLResponse(_LOGIN_PAGE.replace("__ERR__", '<p class="err">Wrong user name or password.</p>'), status_code=401)
    _LOGIN_FAILS.pop(name, None)
    auth.note_sign_in(name)
    resp = RedirectResponse("/", status_code=303)
    resp.set_cookie(auth.COOKIE, auth.make_cookie(name), max_age=auth.MAX_AGE, httponly=True, samesite="lax",
                    secure=os.environ.get("PLASMIDBANK_SECURE_COOKIE", "1") == "1")
    return resp


@app.get("/auth/check")
def auth_check(request: Request):
    """Sign-in for the other services behind Caddy (forward_auth): 200 with who is asking, or a
    redirect to the login page (401 for API calls, which cannot follow it)."""
    user = auth.current_user(request.headers, request.cookies)
    if user and not auth.must_change(user):
        request.state.user = user  # for the usage log
        return PlainTextResponse("ok", headers={"X-Auth-User": user, "X-Auth-Role": auth.role_of(user)})
    if "/api/" in request.headers.get("x-forwarded-uri", ""):
        return JSONResponse({"detail": "not signed in"}, status_code=401)
    return RedirectResponse("/#password" if user else "/login", status_code=303)


@app.get("/logout")
def logout():
    resp = RedirectResponse("/login" if auth.MODE == "local" else "/", status_code=303)
    resp.delete_cookie(auth.COOKIE)
    return resp
_by_id = bank.by_id()
# The whole-bank list (GET /api/maps, no query) cost 0.92 s per call on the VM for 11,827 maps and
# queued at 16 people (load test 2026-10-19). Built once as encoded JSON, dropped by _bank_changed().
_maps_json: bytes | None = None


def _bank_changed() -> None:
    """Call under _BANK_LOCK after any change to bank.maps."""
    global _dups, _by_id, _maps_json
    _dups = {i: g for g in bank.duplicate_groups() for i in g}
    _by_id = bank.by_id()
    _maps_json = None


def _summary(i: int) -> dict:
    if i < 0:
        return _draft_summary(i)
    m = bank.maps[i]
    return {
        "key": i,
        "plasmid_id": m.plasmid_id,
        "kind": m.kind,
        "id_source": m.id_source,
        "file": _display(m),
        "name": m.name,
        "description": m.description,
        "length": m.length,
        "topology": m.topology,
        "format": m.format,
        "status": m.status,
        "n_features": len(m.features),
        "xrefs": m.xrefs,
        "same_sequence_as": [k for k in _dups.get(i, []) if k != i],
        "files_with_same_id": [k for k in _by_id.get(m.plasmid_id or "", []) if k != i],
        "uploaded_by": _uploader(m),
    }


def _uploader(m) -> str | None:
    p = Path(m.path)
    return p.parent.name if p.parent.parent == UPLOADS else None


# A draft in a project opens in the map viewer under a NEGATIVE key, -item id:
# item ids are global and never reused, so the bank's own keys (list positions, >= 0) never collide,
# and every /api/maps/{key}/... tab works on a draft unchanged. The user comes from the request
# (a ContextVar set by the auth middleware), so a draft opens only for its project's members.
_USER: ContextVar[str] = ContextVar("user", default="local")


def _get(key: int):
    if key < 0:
        return _draft_record(key)[1]
    if not 0 <= key < len(bank.maps):
        raise HTTPException(404, "no such map")
    return bank.maps[key]


def _draft_record(key: int) -> tuple[dict, MapRecord]:
    it = _pj(projects.find, -key, _USER.get())
    if it["kind"] != "draft":
        raise HTTPException(404, "no such map")
    m = read_map(projects.root / str(it["project_id"]) / it["draft_file"])
    m.display_path = None
    return it, m


def _draft_summary(key: int) -> dict:
    it, m = _draft_record(key)
    p = projects._project(it["project_id"])
    return {"key": key, "plasmid_id": m.plasmid_id, "kind": m.kind, "id_source": m.id_source,
            "file": "/".join(x for x in (p["name"], it["folder"], it["label"] or it["draft_file"]) if x),
            "name": m.name, "description": m.description, "length": m.length, "topology": m.topology,
            "format": m.format, "status": m.status, "n_features": len(m.features), "xrefs": m.xrefs,
            "same_sequence_as": [], "files_with_same_id": [], "uploaded_by": None,
            "draft": {"project": it["project_id"], "item": it["id"], "folder": it["folder"], "label": it["label"],
                      "personal": bool(p["personal"])}}


@app.get("/")
def home():
    """The lab dashboard. Old links of the form /#<key>/<tab> are forwarded to /bank by its script."""
    return FileResponse(STATIC / "home.html", headers={"Cache-Control": "no-cache"})


@app.get("/bank")
def index():
    return FileResponse(STATIC / "index.html", headers={"Cache-Control": "no-cache"})  # always revalidate: new versions show on reload


@app.get("/changelog")
def changelog():
    """Every release from whatsnew.js on one page, for anyone signed in."""
    return FileResponse(STATIC / "changelog.html", headers={"Cache-Control": "no-cache"})


@app.get("/whatsnew.js")
def whatsnew_js():
    """The What's new card both pages load (once per browser per release, 2026-10-28)."""
    return FileResponse(STATIC / "whatsnew.js", media_type="text/javascript", headers={"Cache-Control": "no-cache"})


@app.get("/tour.js")
def tour_js():
    """The quick guide both pages load (pointer bubbles, shown once per browser)."""
    return FileResponse(STATIC / "tour.js", media_type="text/javascript", headers={"Cache-Control": "no-cache"})


@app.get("/api/links")
def links():
    return {"genomes": GENOMES_URL}


FEEDBACK = STATE / "feedback.jsonl"


@app.post("/api/feedback")
def feedback(request: Request, text: str = Body(..., embed=True), page: str = Body("", embed=True)):
    """A wish or a gripe from the dashboard box: one JSON line each, plus a Discord note if configured."""
    import json
    from datetime import datetime
    text = text.strip()
    if not text:
        raise HTTPException(400, "empty message")
    row = {"time": datetime.now().isoformat(timespec="seconds"), "user": request.state.user,
           "page": page[:200], "text": text[:4000]}
    with FEEDBACK.open("a", encoding="utf-8") as f:
        f.write(json.dumps(row, ensure_ascii=False) + "\n")
    up.discord(f"Feedback from {row['user']}: {row['text'][:300]}")
    return {"ok": True}


@app.get("/api/feedback")
def feedback_list(limit: int = 30):
    import json
    if not FEEDBACK.exists():
        return []
    rows = [json.loads(x) for x in FEEDBACK.read_text(encoding="utf-8").splitlines() if x.strip()]
    return rows[::-1][: max(1, min(limit, 200))]


@app.get("/api/stats")
def stats():
    st: dict[str, int] = {}
    for m in bank.maps:
        st[m.status] = st.get(m.status, 0) + 1
    return {
        "data": str(DATA),
        "maps": len(bank.maps),
        "by_status": st,
        "duplicate_groups": len(bank.duplicate_groups()),
        "oligos": len(oligos),
        "oligo_sequences_matched": len(oligo_groups),
        "oligos_skipped": oligo_skipped,
    }


@app.get("/api/maps")
def maps(q: str = ""):
    global _maps_json
    if q.strip():
        return [_summary(i) for i in bank.search(q)]
    if _maps_json is None:
        with _BANK_LOCK:  # so a list built during an upload is never kept
            if _maps_json is None:
                _maps_json = json.dumps([_summary(i) for i in range(len(bank.maps))]).encode()
    return Response(_maps_json, media_type="application/json")


@app.get("/api/maps/{key}")
def map_detail(key: int):
    m = _get(key)
    return {
        **_summary(key),
        "opened_at": time.time(),  # sent back on save: newer versions by others are flagged
        "seq": m.seq,
        "messages": m.messages,
        "features": [{**asdict(f), "wraps": f.wraps} for f in m.features],
    }


@app.get("/api/maps/{key}/primers")
def primers(key: int, max_mismatches: int = 2, min_binding: int = 15):
    m = _get(key)
    if not m.seq:
        raise HTTPException(400, "map has no sequence")
    return [asdict(h) for h in match_map(m, oligo_groups, max_mismatches, min_binding)]


def _collection_hits(key: int):
    if key < 0:  # a draft changes when it is saved: never cached
        return match_map(_get(key), oligo_groups)
    if key not in _hits_cache:
        _hits_cache[key] = match_map(bank.maps[key], oligo_groups)
    return _hits_cache[key]


@app.get("/api/design/presets")
def design_presets():
    """Golden Gate tail presets, the Tm model, and a suggested (not reserved) next PM- number."""
    nums = [int(o.oligo_id[3:]) for o in oligos if o.oligo_id.startswith("PM-") and o.oligo_id[3:].isdigit()]
    top = max(nums) if nums else 0
    import json
    labels = json.loads((Path(__file__).parents[1] / "oligo_mods.json").read_text(encoding="utf-8"))
    return {"enzymes": GG_PREFIX, "fusion_sites": FUSION_SITES, "tm_model": thermo.MODEL,
            "oligo_mods": labels["mods"], "order": labels["order"],
            "next_pm": f"PM-{top + 1:05d}" if nums else None, "highest_pm": f"PM-{top:05d}" if nums else None,
            "fusion_source": "PhytoBrick common syntax, Patron et al. 2015 New Phytol. doi:10.1111/nph.13532"}


@app.get("/api/maps/{key}/design")
def design(key: int, task: str, start: int, end: int, opt_tm: float = 60.0, max_extra: int = 300,
           tail_f: str = "", tail_r: str = "", sense: str = "+", enzyme: str = "", oh5: str = "", oh3: str = "",
           min_len: int | None = None, max_len: int | None = None):
    m = _get(key)
    if not m.seq:
        raise HTTPException(400, "map has no sequence")
    if not (1 <= start <= m.length and 1 <= end <= m.length):
        raise HTTPException(400, f"positions must be within 1..{m.length}")
    if m.topology == "linear" and end < start:
        raise HTTPException(400, "end before start on a linear map")
    tail_f, tail_r = "".join(tail_f.upper().split()), "".join(tail_r.upper().split())
    if set(tail_f + tail_r) - set("ACGTRYSWKMBDHVN"):
        raise HTTPException(400, "tails may only contain IUPAC bases")
    size = None
    if min_len or max_len:  # annealing part only; tails come on top
        size = (min_len or 15, max_len or 35)
        if not 12 <= size[0] <= size[1] <= 40:
            raise HTTPException(400, "primer length must be 12-40 nt, min no larger than max")
    existing = _collection_hits(key)
    try:
        if task == "pcr":
            return design_pcr(m, start, end, existing, opt_tm, max_extra, tail_f, tail_r, size=size)
        if task == "cloning":
            return design_cloning(m, start, end, existing, opt_tm, tail_f, tail_r, sense, enzyme, oh5, oh3, size=size)
        if task == "sequencing":
            return design_sequencing(m, start, end, existing, opt_tm, size=size)
        if task == "walk":
            return design_walk(m, start, end, existing, opt_tm=opt_tm, size=size)
    except ValueError as e:
        raise HTTPException(400, str(e))
    raise HTTPException(400, "task must be pcr, cloning, sequencing or walk")


# ── BLAST against the bank ─────────────────────────────────────
BLAST_DIR = STATE / "blast"


def _warm_blast() -> None:
    """Build the BLAST database in the background, so the first search does not wait for it."""
    from ..tools import blast
    if blast.available():
        try:
            blast.ensure_db(bank.maps, BLAST_DIR)
        except Exception as e:  # never stop the app over this; a search will report it
            print(f"BLAST database not built: {e}")


import threading as _threading  # noqa: E402

_threading.Thread(target=_warm_blast, daemon=True).start()
_threading.Thread(target=_watch_index, daemon=True).start()


@app.get("/api/blast/status")
def blast_status():
    from ..tools import blast
    return {"blast": blast.available()}


@app.post("/api/blast")
def blast_search(key: int | None = Body(None, embed=True), start: int | None = Body(None, embed=True),
                 end: int | None = Body(None, embed=True), seq: str = Body("", embed=True),
                 max_hits: int = Body(25, embed=True)):
    """BLAST a region of a map (key + start/end, may wrap), a whole map (key), or a pasted sequence."""
    from ..tools import blast
    if not blast.available():
        raise HTTPException(503, "BLAST+ is not installed on this server")
    if not seq:
        if key is None:
            raise HTTPException(400, "give a sequence or a map")
        m = _get(key)
        if not m.seq:
            raise HTTPException(400, "map has no sequence")
        if start and end:
            L = m.length
            n = end - start + 1 if end >= start else end + L - start + 1
            seq = (m.seq + m.seq)[start - 1 : start - 1 + n] if m.circular else m.seq[start - 1 : end]
        else:
            seq = m.seq
    if len(seq) > MAX_QUERY_BP:
        raise HTTPException(413, f"BLAST queries are limited to {MAX_QUERY_BP:,} bp")
    try:
        with _heavy_slot():
            r = blast.search(seq, bank.maps, BLAST_DIR, max(1, min(max_hits, 100)))
    except ValueError as e:
        raise HTTPException(400, str(e))
    except RuntimeError as e:
        raise HTTPException(500, str(e))
    keys = {m.path: i for i, m in enumerate(bank.maps)}
    for h in r["hits"]:
        h["key"] = keys.get(h.pop("path"))
        h["self"] = key is not None and h["key"] == key
    return r


@app.post("/api/build")
def build_map(body: dict = Body(...)):
    """Vector + inserts -> a new map: Gibson primers and fragments, or just the map of a finished cloning.

    body: {vector: key, open: {enzymes: [..]} | {start, end}, mode: gibson|document,
           inserts: [{key, start, end, strand, name} | {seq, name}], overlap, opt_tm, min_len, max_len, name}
    """
    from ..tools import build as B
    try:
        vec = _get(int(body["vector"]))
        keys = {}
        inserts = []
        for i, it in enumerate(body.get("inserts") or []):
            if it.get("seq"):
                inserts.append(B.pasted_insert(it["seq"], it.get("name") or f"fragment {i + 1}"))
            else:
                ref = _get(int(it["key"]))
                keys[id(ref)] = int(it["key"])
                inserts.append(B.region_insert(ref, int(it["start"]), int(it["end"]), it.get("strand", "+"), it.get("name", "")))
        keys[id(vec)] = int(body["vector"])
        size = None
        if body.get("min_len") or body.get("max_len"):
            size = (int(body.get("min_len") or 15), int(body.get("max_len") or 35))
            if not 12 <= size[0] <= size[1] <= 40:
                raise ValueError("primer length must be 12-40 nt, min no larger than max")
        op = body.get("open") or {}
        vec_open = {"enzymes": [e for e in op.get("enzymes", []) if e]} if op.get("enzymes") else \
            {"start": int(op.get("start", 0)), "end": int(op.get("end", 0))}
        return B.build(vec, vec_open, inserts, body.get("mode", "gibson"), int(body.get("overlap") or B.OVERLAP),
                       float(body.get("opt_tm") or 60), size, lambda ref: _collection_hits(keys[id(ref)]) if id(ref) in keys else [],
                       body.get("name", ""))
    except (KeyError, TypeError) as e:
        raise HTTPException(400, f"incomplete request: {e}")
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.post("/api/build/from_sequence")
def build_from_sequence(text: str = Body(..., embed=True), name: str = Body("", embed=True),
                        circular: bool = Body(True, embed=True)):
    """A pasted sequence (bare, FASTA or GenBank) -> a GenBank file, ready to upload or save."""
    from ..tools.build import map_from_text
    try:
        return map_from_text(text, name, circular)
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.get("/api/maps/{key}/sdm")
def design_mutagenesis(key: int, start: int = 0, end: int = 0, new: str = "", feature: int = -1, mutation: str = "",
                       opt_tm: float = 60.0, min_len: int | None = None, max_len: int | None = None):
    """Q5-style site-directed mutagenesis: a DNA change (start, end, new) or a protein change in a CDS."""
    from ..tools.mutagenesis import design_sdm, protein_change
    m = _get(key)
    if not m.seq:
        raise HTTPException(400, "map has no sequence")
    size = None
    if min_len or max_len:
        size = (min_len or 15, max_len or 35)
        if not 12 <= size[0] <= size[1] <= 40:
            raise HTTPException(400, "primer length must be 12-40 nt, min no larger than max")
    try:
        protein = None
        if mutation:
            if not 0 <= feature < len(m.features) or not m.features[feature].strand:
                raise ValueError("pick a coding feature (with a strand) for a protein change")
            protein = protein_change(m, m.features[feature], mutation)
            start, end, new = protein["start"], protein["end"], protein["new"]
        return design_sdm(m, start, end, new, _collection_hits(key), opt_tm, size, protein)
    except ValueError as e:
        raise HTTPException(400, str(e))


def _clean(key: int, rules: str):
    m = _get(key)
    if m.format != "genbank":
        raise HTTPException(400, "the cleaner works on GenBank files only")
    chosen = [r for r in rules.split(",") if r] if rules else DEFAULT_ON
    try:
        return clean_record(read_genbank(m.path), chosen, _oligo_names)
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.get("/api/clean/rules")
def clean_rules():
    return {"default_on": DEFAULT_ON, "opt_in": OPT_IN}


@app.get("/api/maps/{key}/clean")
def clean_preview(key: int, rules: str = ""):
    """Dry run: what the cleaner would change. Nothing is written anywhere."""
    r = _clean(key, rules)
    return {"n_before": r.n_before, "n_after": r.n_after,
            "changes": [asdict(c) for c in r.changes], "flags": [asdict(c) for c in r.flags]}


@app.get("/api/maps/{key}/clean.gb")
def clean_download(key: int, rules: str = ""):
    """The cleaned copy as a download; the original file is not modified."""
    r = _clean(key, rules)
    name = Path(_get(key).path).stem.split("__", 1)[-1] + "_cleaned.gb"
    return PlainTextResponse(r.genbank(), headers={"Content-Disposition": f'attachment; filename="{name}"'})


@app.get("/api/maps/{key}/orfs")
def orfs(key: int, min_aa: int = 100):
    m = _get(key)
    return find_orfs(m, min_aa) if m.seq else []


@app.get("/api/pcr/polymerases")
def pcr_polymerases():
    from ..tools.pcr import POLYMERASES, TEMPLATES
    return {"polymerases": [{"id": k, "name": v["name"], "confidence": v["confidence"]} for k, v in POLYMERASES.items()],
            "templates": list(TEMPLATES)}


@app.post("/api/maps/{key}/pcr")
def pcr_product(key: int, primers: list[dict] = Body(..., embed=True), polymerase: str = Body("Q5", embed=True),
                template: str = Body("plasmid", embed=True)):
    """Two primers (lab oligo IDs or typed sequences) -> the product and a cycling program."""
    from ..tools.pcr import pcr
    from ..tools.build import map_from_text
    m = _get(key)
    by_id = {o.oligo_id: o for o in oligos}
    resolved = []
    for p in primers:
        if p.get("oligo_id"):
            o = by_id.get(p["oligo_id"])
            if o is None:
                raise HTTPException(400, f"no oligo {p['oligo_id']} in the lab's collection")
            resolved.append({"name": f"{o.oligo_id} {o.name}".strip(), "seq": o.seq})
        else:
            seq = "".join(c for c in str(p.get("seq", "")).upper() if c.isalpha())
            if not 12 <= len(seq) <= 200:
                raise HTTPException(400, "a typed primer must be 12–200 bases")
            resolved.append({"name": (p.get("name") or "").strip()[:40] or f"primer {len(resolved) + 1}", "seq": seq})
    try:
        r = pcr(m, resolved, polymerase, template)
    except ValueError as e:
        raise HTTPException(400, str(e))
    name = f"PCR_{r['forward']['name'].split()[0]}_{r['reverse']['name'].split()[0]}"
    r["genbank"] = map_from_text(r["product_seq"], name, circular=False)["genbank"]
    r["name"] = name
    return r


@app.post("/api/maps/{key}/join")
def join(key: int, indices: list[int] = Body(..., embed=True)):
    """Exons (features) joined into one coding sequence: frame, translation, stops."""
    from ..tools.analysis import join_features
    try:
        return join_features(_get(key), indices)
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.get("/api/maps/{key}/translate")
def translate(key: int, start: int, end: int, strand: str = "+"):
    m = _get(key)
    if not (1 <= start <= m.length and 1 <= end <= m.length):
        raise HTTPException(400, f"positions must be within 1..{m.length}")
    return translate_region(m, start, end, strand)


@app.get("/api/maps/{key}/digest")
def digest_map(key: int, enzymes: str):
    m = _get(key)
    names = [e for e in enzymes.replace(",", " ").split() if e]
    try:
        frags, cuts = digest(m, names)
    except ValueError as e:
        raise HTTPException(400, str(e))
    return {"cuts": cuts, "fragments": [{"length": f.length, "left": asdict(f.left), "right": asdict(f.right),
                                          "start": f.src_start % m.length + 1,
                                          "features": [x["label"] for x in f.features]} for f in frags]}


@app.post("/api/gel")
def gel(lanes: list[dict] = Body(...), ladder: str = Body("1 kb (NEB N3232)")):
    out = []
    for lane in lanes[:15]:
        try:
            out.append(gel_lane(_get(int(lane["key"])), [e for e in lane.get("enzymes", []) if e]))
        except ValueError as e:
            out.append({"label": str(lane), "bands": [], "error": str(e)})
    if ladder not in LADDERS:
        raise HTTPException(400, f"ladder must be one of {', '.join(LADDERS)}")
    return {"ladder": {"label": ladder, "bands": LADDERS[ladder]}, "lanes": out, "ladders": list(LADDERS),
            "note": "ladder sizes UNVERIFIED against vendor sheets; band positions are a log-scale approximation"}


def _fusion_names() -> dict[str, list[str]]:
    """Overhang -> where it sits in the PhytoBrick common syntax ('end of 5'UTR', 'start of CDS'...)."""
    out: dict[str, list[str]] = {}
    for part, (a, b) in FUSION_SITES.items():
        out.setdefault(a, []).append(f"start of {part}")
        out.setdefault(b, []).append(f"end of {part}")
    return out


def _input_record(request: Request, ref: dict) -> MapRecord:
    """One assembly input: {"key": bank key} or {"project": id, "item": id} (a draft or a link)."""
    if "key" in ref:
        return _get(int(ref["key"]))
    pid, iid = int(ref["project"]), int(ref["item"])
    it = _pj(projects.item, pid, request.state.user, iid)
    if it["kind"] == "draft":
        return read_map(_pj(projects.draft_path, pid, request.state.user, iid))
    k = _resolve(it)
    if k is None:
        raise HTTPException(404, f"{it['plasmid_id'] or it['path']} is not in the bank (yet)")
    return bank.maps[k]


@app.post("/api/assemble/golden_gate")
def assemble_gg(request: Request, acceptor: int | None = Body(None), parts: list[int] = Body([]),
                inputs: list[dict] = Body([]), enzyme: str = Body("BsaI")):
    """Inputs are symmetric: every one gives the fragment that no longer carries a site, and the
    fragments are tried in every order and orientation. 'acceptor' is kept for older pages."""
    refs = ([{"key": acceptor}] if acceptor is not None else []) + [{"key": k} for k in parts] + list(inputs)
    if len(refs) < 2:
        raise HTTPException(400, "choose at least two plasmids")
    if len(refs) > MAX_GG_INPUTS:
        raise HTTPException(400, f"at most {MAX_GG_INPUTS} plasmids in one Golden Gate assembly")
    recs = [_input_record(request, x) for x in refs]
    for m in recs:
        if not m.seq:
            raise HTTPException(400, f"{m.plasmid_id or m.name} has no sequence (CLC file?)")
    try:
        with _heavy_slot():
            r = golden_gate(recs[0], recs[1:], enzyme)
    except ValueError as e:
        raise HTTPException(400, str(e))
    names = _fusion_names()
    for i, p in enumerate(r["products"]):
        # Read the product on the strand where its junctions match the standard fusion sites,
        # so overhangs show as GGAG/AATG..., not their reverse complements.
        fwd = sum(j["overhang"].upper() in names for j in p.junctions)
        rev = sum(revcomp(j["overhang"]).upper() in names for j in p.junctions)
        if rev > fwd or (rev == fwd and p.parts[0]["orientation"] == "-"):
            p = r["products"][i] = flip_product(p)
        for k, j in enumerate(p.junctions):
            j["site"] = " / ".join(names.get(j["overhang"].upper(), []))
            j["pos"] = p.parts[(k + 1) % len(p.parts)]["start"]  # 0-based start of the overhang on the product
        flipped = [x["source"] for x in p.parts if x["orientation"] == "-"]
        if flipped:
            p.warnings.append("used reverse-complemented relative to how its map is drawn: " + ", ".join(flipped))
        p.fusions = fusion_frames(p)  # CDS parts joined end to end must stay in frame
        p.warnings += [f["message"] for f in p.fusions]
    for d in r.get("diagnosis", []):  # a part read the other way round shows the reverse complement
        oh = d["overhang"].upper()
        d["site"] = " / ".join(names.get(oh, []) or [f"{x} (read reversed)" for x in names.get(revcomp(oh), [])])
    name = "GG_" + "_".join((m.plasmid_id or m.name or "x") for m in recs)
    return {**{k: v for k, v in r.items() if k != "products"},
            "products": [{"length": len(p.seq), "parts": p.parts, "junctions": p.junctions, "warnings": p.warnings,
                          "fusions": getattr(p, "fusions", []),
                          "features": p.features, "seq": p.seq,
                          "genbank": product_genbank(p, name[:16], f"Golden Gate ({enzyme}) {name}")}
                         for p in r["products"]]}


LAB_ENZYMES = STATE / "lab_enzymes.json"
_ENZ_LOCK = threading.Lock()


def _lab_enzymes() -> dict:
    import json
    if not LAB_ENZYMES.exists():
        return {"enzymes": [], "updated_by": None, "updated": None, "version": 0}
    return json.loads(LAB_ENZYMES.read_text(encoding="utf-8"))


@app.get("/api/lab/enzymes")
def lab_enzymes_get():
    """The enzymes in the lab's freezer: one shared list, editable by anyone signed in."""
    return _lab_enzymes()


@app.put("/api/lab/enzymes")
def lab_enzymes_put(request: Request, enzymes: str = Body(..., embed=True), version: int | None = Body(None, embed=True)):
    """`version` is the one the editor loaded; a stale one is refused so nobody's change is lost silently."""
    import json
    from datetime import datetime
    try:
        names = enzymes_for("common", enzymes) if enzymes.strip() else []
    except ValueError as e:
        raise HTTPException(400, str(e))
    with _ENZ_LOCK:
        cur = _lab_enzymes()
        if version is not None and version != cur.get("version", 0):
            raise HTTPException(409, f"{cur.get('updated_by') or 'Someone'} changed the list at {cur.get('updated')} "
                                     "after you opened it. Nothing was saved: reload to see their change, then add yours.")
        row = {"enzymes": sorted(set(names)), "updated_by": request.state.user,
               "updated": datetime.now().isoformat(timespec="seconds"), "version": cur.get("version", 0) + 1}
        tmp = LAB_ENZYMES.with_suffix(".tmp")
        tmp.write_text(json.dumps(row, indent=1), encoding="utf-8")
        tmp.replace(LAB_ENZYMES)
    return row


@app.post("/api/enzymes/select")
def enzymes_select(request: Request, conditions: list[dict] = Body(...), set: str = Body("commercial"),
                   enzymes: str = Body(""), one_per_site: bool = Body(True)):
    """ApE-style finder: e.g. once in the vector AND absent from the insert, from the lab's enzymes."""
    if not conditions:
        raise HTTPException(400, "add at least one map")
    lab = _lab_enzymes()["enzymes"]
    if enzymes.strip():
        try:
            names = enzymes_for("common", enzymes)
        except ValueError as e:
            raise HTTPException(400, str(e))
    elif set == "lab":
        if not lab:
            raise HTTPException(400, "the lab enzyme list is empty: add the freezer's enzymes first")
        names = lab
    elif set == "commercial":
        names = commercial()
    elif set in SETS:
        names = SETS[set]
    else:
        raise HTTPException(400, f"unknown set {set!r}")
    conds, labels = [], []
    for c in conditions:
        m = _input_record(request, c)
        if not m.seq:
            raise HTTPException(400, f"{m.plasmid_id or m.name} has no sequence (CLC file?)")
        lo, hi = int(c.get("min", 0)), c.get("max")
        reg = (int(c["region_start"]), int(c["region_end"])) if c.get("region_start") and c.get("region_end") else None
        conds.append({"map": m, "min": lo, "max": None if hi in (None, "") else int(hi), "region": reg})
        labels.append(m.plasmid_id or m.name)
    rows = select_enzymes(conds, names, one_per_site, prefer=set_prefer(lab))
    return {"maps": labels, "searched": len(names), "enzymes": rows}


def set_prefer(lab: list[str]) -> list[str]:
    """Which name to show for a group of isoschizomers: the lab's own first, then the common set."""
    return list(lab) + [n for n in COMMON if n not in lab]


@app.get("/api/align/available")
def align_available():
    return {"mafft": mafft_available()}


MAX_UPLOAD_FILES = int(os.environ.get("PLASMIDBANK_MAX_UPLOAD_FILES", 20))
MAX_UPLOAD_MB = int(os.environ.get("PLASMIDBANK_MAX_UPLOAD_MB", 10))  # per file; Caddy caps the request at 50 MB


async def _upload_bytes(f: UploadFile) -> bytes:
    """An uploaded file's bytes, refused (413) past MAX_UPLOAD_MB without reading the rest."""
    limit = MAX_UPLOAD_MB * 1024 * 1024
    data = await f.read(limit + 1)
    if len(data) > limit:
        raise HTTPException(413, f"{f.filename or 'file'} is larger than {MAX_UPLOAD_MB} MB")
    return data


def _check_file_count(n: int) -> None:
    if n > MAX_UPLOAD_FILES:
        raise HTTPException(413, f"at most {MAX_UPLOAD_FILES} files at a time ({n} sent)")


@app.post("/api/align/multiple")
async def align_multiple(keys: str = Form(""), files: list[UploadFile] = File(default=[])):
    recs = []
    keys_ = [k for k in keys.split(",") if k.strip()]
    _check_file_count(len(keys_) + len(files))
    for k in keys_:
        m = _get(int(k))
        if m.seq:
            recs.append((m.plasmid_id or Path(m.path).stem, m.seq, m.circular))
    for f in files:
        for q in load_queries(f.filename or "upload", await _upload_bytes(f), trim=True):
            recs.append((q.name, q.seq, False))
    if any(len(r[1]) > MAX_QUERY_BP for r in recs) or sum(len(r[1]) for r in recs) > MAX_ALIGN_TOTAL_BP:
        raise HTTPException(413, f"alignments are limited to {MAX_QUERY_BP:,} bp per sequence, "
                                 f"{MAX_ALIGN_TOTAL_BP:,} bp in total")
    try:
        with _heavy_slot():
            return await run_in_threadpool(multiple_alignment, recs)
    except (RuntimeError, ValueError) as e:
        raise HTTPException(400, str(e))


OUTBOX = Path(os.environ.get("PLASMIDBANK_OUTBOX", Path.home() / "PlasmidBank" / "outbox"))


@app.post("/api/maps/{key}/save")
def save_map(request: Request, key: int, edits: dict = Body(...)):
    """Apply feature/name/topology edits and save a NEW file to the outbox (never over the original).

    Guards: `expect_file` must still be the file at `key` (keys move when an upload is
    removed); and a version saved by someone else after `opened_at` gets a 409 listing it, unless
    `force` (the user saw it and saves anyway)."""
    expect, opened, force = edits.pop("expect_file", None), edits.pop("opened_at", None), edits.pop("force", False)
    if key < 0:
        return _save_draft(key, edits, opened, force)
    m = _get(key)
    if expect is not None and _display(m) != expect:
        raise HTTPException(409, "The map bank changed since you opened this map, so this number now points to another "
                                 "map. Nothing was saved: reload the map and make the edit again.")
    if m.format != "genbank":
        raise HTTPException(400, "editing works on GenBank files (convert FASTA/CLC first)")
    try:
        text, summary = apply_edits(m.path, edits)
    except ValueError as e:
        raise HTTPException(400, str(e))
    from datetime import datetime
    box = _outbox(request)
    prefix = m.plasmid_id or Path(m.path).stem
    if opened is not None and not force:
        newer = _newer_versions(prefix, box, float(opened))
        if newer:
            who = ", ".join(f"{n['user']} at {n['saved']}" for n in newer)
            return JSONResponse({"detail": f"Someone saved a newer version of {prefix} after you opened it: {who}.",
                                 "newer": newer}, status_code=409)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    dest, n = box / f"{prefix}__edited_{stamp}.gb", 1
    while True:  # two saves in one second must not overwrite each other
        try:
            with open(dest, "x", encoding="utf-8") as f:
                f.write(text)
            break
        except FileExistsError:
            n += 1
            dest = box / f"{prefix}__edited_{stamp}-{n}.gb"
    return {"saved_to": f"{request.state.user}/{dest.name}", "summary": summary, "genbank": text,
            "saved_at": time.time()}


def _save_draft(key: int, edits: dict, opened: float | None = None, force: bool = False) -> dict:
    """A draft is the user's own working file: saved in place (the bank is never edited in place).
    The share sync can change it while it is open (an edit made on the share, pulled in): a save based on the
    older copy gets a 409 unless `force`, so the pulled edit is never overwritten unseen."""
    it, m = _draft_record(key)
    if opened is not None and not force:
        from datetime import datetime
        t = Path(m.path).stat().st_mtime
        if t > float(opened):
            return JSONResponse({"detail": "This file changed after you opened it (an edit made on the lab share was "
                                           "brought in by the sync). Saving now replaces that edit.",
                                 "newer": [{"user": "the lab share", "saved": datetime.fromtimestamp(t).strftime("%H:%M %d %b"),
                                            "t": t}]}, status_code=409)
    if m.format != "genbank":
        raise HTTPException(400, "editing works on GenBank files (convert FASTA/CLC first)")
    try:
        text, summary = apply_edits(m.path, edits)
        projects.replace_draft(it["project_id"], _USER.get(), it["id"], text.encode())
    except ValueError as e:
        raise HTTPException(400, str(e))
    return {"saved_to": _draft_summary(key)["file"], "summary": summary, "genbank": text, "saved_at": time.time(),
            "in_place": True}


@app.post("/api/maps/{key}/overwrite")
async def overwrite_on_o(request: Request, key: int, confirm_id: str = Form(...), edits: str = Form(""),
                         file: UploadFile | None = File(None)):
    """Overwrite a bank map on the share (any map, after confirmations; the old file is kept in
    `_replaced/` beside it). The new content is either the edits made in the app or an uploaded GenBank
    file (a CLC export for a .clc). `confirm_id` must repeat the map's ID (the second confirmation is
    checked here, not only in the browser). The hourly sync carries it out (uploads.apply_overwrites),
    and refuses if the file on the share changed after it was read."""
    import hashlib
    import uuid
    from datetime import datetime
    data = await _upload_bytes(file) if file is not None and file.filename else None
    if key < 0:
        raise HTTPException(400, "this is your own file: Save writes it in place")
    if auth.role_of(request.state.user) == "student":
        raise HTTPException(403, "students cannot overwrite maps on the share; save a copy to My workspace instead")
    m = _get(key)
    if _uploader(m) is not None:
        raise HTTPException(400, "this upload is not on the share yet; remove it from the upload queue and upload the new version")
    rel = _display(m)
    if Path(rel).is_absolute() or Path(rel).parts[0] in ("Workspaces", "uploads"):
        raise HTTPException(400, "only maps in the bank on the share can be overwritten")
    expected = m.plasmid_id or Path(rel).name
    if confirm_id.strip().upper() != expected.upper():
        raise HTTPException(400, f"type {expected} exactly to confirm the overwrite")
    if any(r.get("status") == "queued" and r["rel"] == rel for r in _overwrite_requests()):
        raise HTTPException(409, "a replacement for this map is already waiting for the hourly sync")
    if data is not None:
        text = data.decode("utf-8", errors="replace")
    elif edits:
        if m.format != "genbank":
            raise HTTPException(400, "this map is not GenBank (e.g. a .clc): upload its GenBank export instead")
        try:
            text, _ = apply_edits(m.path, json.loads(edits))
        except (ValueError, json.JSONDecodeError) as e:
            raise HTTPException(400, str(e))
    else:
        raise HTTPException(400, "nothing to write: make edits or choose a GenBank file")
    rid = datetime.now().strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
    OVERWRITES.mkdir(parents=True, exist_ok=True)
    (OVERWRITES / f"{rid}.gb").write_text(text, encoding="utf-8")
    try:
        new = read_map(OVERWRITES / f"{rid}.gb")
        if new.format != "genbank" or not new.seq:
            raise ValueError("not a GenBank file with a sequence")
    except Exception as e:
        (OVERWRITES / f"{rid}.gb").unlink(missing_ok=True)
        raise HTTPException(400, f"not a usable map: {e}")
    name = Path(rel).name
    new_name = Path(name).stem + ".gbk" if name.lower().endswith(".clc") else name
    req = {"id": rid, "rel": rel, "new_name": new_name, "file": f"{rid}.gb", "user": request.state.user,
           "requested": datetime.now().isoformat(timespec="seconds"), "status": "queued",
           "sha_before": hashlib.sha256(Path(m.path).read_bytes()).hexdigest(), "plasmid_id": m.plasmid_id}
    (OVERWRITES / f"{rid}.json").write_text(json.dumps(req, indent=1))
    with _BANK_LOCK:
        _overlay_overwrites(bank)
        _bank_changed()
    return {"queued": rid, "rel": rel, "new_name": new_name, "key": key}


@app.post("/api/maps/{key}/copy")
def copy_to_workspace(request: Request, key: int, body: dict = Body(...)):
    """Save this map, with the edits not yet saved, as a new file in My workspace.
    A copy: the bank map is never touched."""
    m = _get(key)
    edits = body.get("edits") or {}
    if m.format == "genbank":
        try:
            text, _ = apply_edits(m.path, edits)
        except ValueError as e:
            raise HTTPException(400, str(e))
    elif not edits:
        text = Path(m.path).read_text(encoding="utf-8", errors="replace")
    else:
        raise HTTPException(400, "editing works on GenBank files (convert FASTA/CLC first)")
    user = request.state.user
    pid = projects.personal(user)
    label = (body.get("label") or "").strip() or (m.plasmid_id or Path(m.path).stem.split("__", 1)[-1]) + " copy"
    ext = ".gb" if m.format == "genbank" else Path(m.path).suffix or ".txt"
    iid = _pj(projects.add_draft, pid, user, label + ext, text.encode(), label, body.get("folder") or "")
    return {"project": pid, "item": iid, "key": -iid}


@app.post("/api/projects/{pid}/items/{iid}/move")
def project_move_item(request: Request, pid: int, iid: int, folder: str | None = Body(None, embed=True),
                      label: str | None = Body(None, embed=True)):
    """Move a draft to another folder of My workspace and/or rename it."""
    _pj(projects.move, pid, request.state.user, iid, folder, label)
    return {"ok": True}


@app.post("/api/projects/{pid}/folders/rename")
def project_rename_folder(request: Request, pid: int, old: str = Body(..., embed=True), new: str = Body(..., embed=True)):
    """Rename (or move) a My workspace folder with everything in it."""
    return {"moved": _pj(projects.rename_folder, pid, request.state.user, old, new)}


def _newer_versions(prefix: str, own_box: Path, since: float) -> list[dict]:
    """Versions of this map in OTHER users' outboxes saved after `since` (none in single-user mode)."""
    import glob
    from datetime import datetime
    if own_box == OUTBOX:
        return []
    out = []
    for p in OUTBOX.glob(f"*/{glob.escape(prefix)}__edited_*.gb"):
        t = p.stat().st_mtime
        if p.parent != own_box and t > since:
            out.append({"user": p.parent.name, "saved": datetime.fromtimestamp(t).strftime("%H:%M %d %b"), "t": t})
    return sorted(out, key=lambda x: x["t"])


def _outbox(request: Request) -> Path:
    """Each user's own folder; in single-user mode ('local') the outbox itself."""
    user = request.state.user
    box = OUTBOX if user == "local" else OUTBOX / auth.safe_name(user)
    box.mkdir(parents=True, exist_ok=True)
    return box


def _upload_user(request: Request) -> str:
    user = request.state.user
    return "local" if user == "local" else auth.safe_name(user)


@app.post("/api/upload")
async def upload_map(request: Request, file: UploadFile = File(...), plasmid_id: str = Form(""),
                     construct_name: str = Form(""), folder: str = Form("")):
    """A new map: stored on the server at once, copied to the share by the hourly sync (add-only). Without a
    lab ID (blank, 'Test', 'TMP-3', a project code) it becomes a private draft in the user's
    "My drafts" project instead: never in the bank, never on the share."""
    data = await _upload_bytes(file)
    if _is_draft_tag(plasmid_id):
        return await run_in_threadpool(_store_draft, request, data, file.filename or "map.gb", plasmid_id,
                                       construct_name, folder)
    return await run_in_threadpool(_store_new_map, request, data, plasmid_id, construct_name)


def _is_draft_tag(text: str) -> bool:
    """Not a lab ID. Something that starts like one ('EC123', 'PMC-1234567') is a typo, not a
    draft: it stays an error, so a mistyped number never quietly misses the share."""
    try:
        up.normalise_id(text)
        return False
    except ValueError:
        return not re.match(r"\s*(PMC|EC)[-_ ]?\d", text, re.I)


def _store_draft(request: Request, data: bytes, filename: str, tag: str, construct_name: str,
                 folder: str = "") -> dict:
    user = request.state.user
    pid = projects.personal(user)
    label = " ".join(x for x in (tag.strip(), construct_name.strip()) if x)
    iid = _add_checked_draft(pid, user, filename, data, label, folder)
    return {"draft": True, "project": pid, "item": iid, "key": -iid, "file": label or Path(filename).stem}


# Every change to bank.maps (upload, promote, removing an upload) holds this lock, so two uploads of
# one new ID cannot both pass the "already there?" checks, and keys don't move mid-change.
# One process only: a second uvicorn worker would need a file lock.
_BANK_LOCK = threading.Lock()


def _store_new_map(request: Request, data: bytes, plasmid_id: str, construct_name: str) -> dict:
    with _BANK_LOCK:
        return _store_new_map_locked(request, data, plasmid_id, construct_name)


def _store_new_map_locked(request: Request, data: bytes, plasmid_id: str, construct_name: str) -> dict:
    if len(data) > 5 * 1024 * 1024:
        raise HTTPException(400, "file too large for a plasmid map (5 MB max)")
    from_o = {m.plasmid_id for m in bank.maps if m.plasmid_id and _uploader(m) is None}
    try:
        dest = up.store_upload(UPLOADS, _upload_user(request), data, plasmid_id, construct_name, from_o)
    except ValueError as e:
        raise HTTPException(400, str(e))
    annotated = 0
    if annotations.on_upload(request.state.user):  # the user's profile setting
        annotated = _annotate_file(dest, request.state.user)
    m = _read_upload(dest)
    for i, old in enumerate(bank.maps):  # the user's own earlier upload of this ID is replaced
        if _uploader(old) and old.plasmid_id == m.plasmid_id:
            bank.maps[i] = m
            key = i
            break
    else:
        bank.maps.append(m)
        key = len(bank.maps) - 1
    _bank_changed()
    same = [bank.maps[k].plasmid_id or Path(bank.maps[k].path).name for k in _dups.get(key, []) if k != key]
    return {"key": key, "file": dest.name, "same_sequence_as": same, "plasmid_id": m.plasmid_id, "annotated": annotated}


def _annotate_file(path: Path, user: str) -> int:
    """Add the user's common + private annotations to an uploaded GenBank file, in place. FASTA has no
    features, so it is left alone. Returns how many features were added."""
    from ..annotations import annotate
    m = read_map(path)
    if m.format != "genbank":
        return 0
    found = [p for p in annotate(m, annotations.entries(user)) if not p["same_span_annotated"]]
    if not found:
        return 0
    text, _ = apply_edits(str(path), {"add": [{k: p[k] for k in ("start", "end", "strand", "type", "label")} for p in found]})
    path.write_text(text, encoding="utf-8")
    return len(found)


# ── annotation lists ─────────────────────────────────────────────
@app.get("/api/annotations")
def annotation_list(request: Request):
    from ..annotations import TYPES
    u = request.state.user
    return {"entries": annotations.entries(u), "types": TYPES, "annotate_on_upload": annotations.on_upload(u)}


@app.post("/api/annotations")
def annotation_add(request: Request, label: str = Body(..., embed=True), type: str = Body("misc_feature", embed=True),
                   seq: str = Body("", embed=True), scope: str = Body("mine", embed=True), note: str = Body("", embed=True),
                   key: int | None = Body(None, embed=True), start: int | None = Body(None, embed=True),
                   end: int | None = Body(None, embed=True), strand: int = Body(1, embed=True), maps: int | None = Body(None, embed=True)):
    """A new entry from a sequence, or from a stretch of a map (key + start/end + strand)."""
    if not seq and key is not None and start and end:
        m = _get(key)
        L = m.length
        n = end - start + 1 if end >= start else end + L - start + 1
        seq = ((m.seq + m.seq) if m.circular else m.seq)[start - 1 : start - 1 + n]
        if strand < 0:
            seq = revcomp(seq)
    try:
        return {"id": annotations.add(request.state.user, "common" if scope == "common" else "mine", label, type, seq, note, maps)}
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.put("/api/annotations/settings")
def annotation_settings(request: Request, annotate_on_upload: bool = Body(..., embed=True)):
    annotations.set_on_upload(request.state.user, annotate_on_upload)
    return {"annotate_on_upload": annotate_on_upload}


@app.put("/api/annotations/{entry_id}")
def annotation_update(request: Request, entry_id: int, label: str | None = Body(None, embed=True),
                      type: str | None = Body(None, embed=True), note: str | None = Body(None, embed=True)):
    try:
        annotations.update(request.state.user, entry_id, label, type, note)
    except LookupError as e:
        raise HTTPException(404, str(e))
    except ValueError as e:
        raise HTTPException(400, str(e))
    return {"ok": True}


@app.delete("/api/annotations/{entry_id}")
def annotation_delete(request: Request, entry_id: int):
    try:
        annotations.delete(request.state.user, entry_id)
    except LookupError as e:
        raise HTTPException(404, str(e))
    return {"ok": True}


@app.get("/api/annotations/suggest")
def annotation_suggest(request: Request, min_maps: int = 3, limit: int = 150):
    """Parts the bank uses most often that are not in the lists yet (the seed for the common list)."""
    from ..annotations import suggest_from_bank
    known = {e["seq"] for e in annotations.entries(request.state.user)}
    return suggest_from_bank(bank.maps, known, max(2, min_maps), max(1, min(limit, 500)))


@app.post("/api/maps/{key}/annotate")
def annotate_map(request: Request, key: int):
    """Proposals only: nothing is saved until the user accepts them (Features -> save as new version)."""
    from ..annotations import annotate
    m = _get(key)
    return {"proposals": annotate(m, annotations.entries(request.state.user)), "entries": len(annotations.entries(request.state.user))}


# ── projects (inside the app; the project owner invites members) ──────────
def _pj(fn, *args):
    try:
        return fn(*args)
    except Forbidden as e:
        raise HTTPException(403, str(e))
    except KeyError as e:
        raise HTTPException(404, str(e).strip("'"))
    except ValueError as e:
        raise HTTPException(400, str(e))


def _resolve(it: dict) -> int | None:
    """Bank key of a linked map, preferring the copy on the share over a not-yet-synced upload."""
    if it["plasmid_id"]:
        keys = _by_id.get(it["plasmid_id"], [])
        on_o = [k for k in keys if _uploader(bank.maps[k]) is None]
        return (on_o or keys or [None])[0]
    if it["path"]:
        return next((i for i, m in enumerate(bank.maps) if _display(m) == it["path"]), None)
    return None


def _item_view(pid: int, it: dict, user: str) -> dict:
    out = {k: it[k] for k in ("id", "kind", "plasmid_id", "path", "label", "added_by", "added", "promoted", "folder")}
    if it["kind"] == "ref":
        k = _resolve(it)
        out["key"] = k
        if k is not None:
            s = _summary(k)
            out.update(name=s["name"], length=s["length"], file=s["file"], topology=s["topology"])
    else:
        try:
            m = read_map(projects.draft_path(pid, user, it["id"]))
            out.update(name=m.name, length=m.length, file=it["draft_file"], topology=m.topology, key=-it["id"])
        except Exception as e:  # unreadable draft: still listed, so it can be removed
            out["error"] = str(e)
    return out


@app.get("/api/projects")
def project_list(request: Request, plasmid_id: str = ""):
    """The user's projects; with plasmid_id, each also says which item links that map (or null)."""
    rows = projects.for_user(request.state.user)
    if plasmid_id:
        has = projects.refs_of(request.state.user, plasmid_id)
        for r in rows:
            r["item_for"] = has.get(r["id"])
    return rows


@app.post("/api/projects")
def project_create(request: Request, name: str = Body(..., embed=True), description: str = Body("", embed=True)):
    return {"id": _pj(projects.create, request.state.user, name, description)}


@app.get("/api/projects/{pid}")
def project_detail(request: Request, pid: int):
    user = request.state.user
    p = _pj(projects.get, pid, user)
    p["items"] = [_item_view(pid, it, user) for it in p["items"]]
    p["is_owner"] = p["owner"] == user
    return p


@app.delete("/api/projects/{pid}")
def project_delete(request: Request, pid: int):
    _pj(projects.delete, pid, request.state.user)
    return {"ok": True}


@app.post("/api/projects/{pid}/members")
def project_add_member(request: Request, pid: int, user: str = Body(..., embed=True)):
    name = user.strip().lower()
    if not name:
        raise HTTPException(400, "enter a user name")
    _pj(projects.add_member, pid, request.state.user, name)
    return {"ok": True}


@app.delete("/api/projects/{pid}/members/{member}")
def project_remove_member(request: Request, pid: int, member: str):
    _pj(projects.remove_member, pid, request.state.user, member)
    return {"ok": True}


@app.post("/api/projects/{pid}/refs")
def project_add_ref(request: Request, pid: int, key: int = Body(..., embed=True)):
    m = _get(key)
    return {"id": _pj(projects.add_ref, pid, request.state.user, m.plasmid_id, None if m.plasmid_id else _display(m),
                      m.name or "")}


@app.post("/api/projects/{pid}/drafts")
async def project_add_draft(request: Request, pid: int, file: UploadFile = File(...), label: str = Form(""),
                            folder: str = Form("")):
    data = await _upload_bytes(file)
    return {"id": _add_checked_draft(pid, request.state.user, file.filename or "map.gb", data, label, folder)}


def _add_checked_draft(pid: int, user: str, filename: str, data: bytes, label: str, folder: str = "") -> int:
    iid = _pj(projects.add_draft, pid, user, filename, data, label, folder)
    try:
        read_map(projects.draft_path(pid, user, iid))
    except Exception as e:
        projects.remove_item(pid, user, iid)
        raise HTTPException(400, f"not a readable map: {e}")
    return iid


@app.post("/api/projects/{pid}/drafts/text")
def project_add_draft_text(request: Request, pid: int, filename: str = Body(..., embed=True),
                           text: str = Body(..., embed=True), label: str = Body("", embed=True),
                           folder: str = Body("", embed=True)):
    """Save a generated map (e.g. a Golden Gate product) straight into a project as a draft."""
    return {"id": _pj(projects.add_draft, pid, request.state.user, filename, text.encode(), label, folder)}


@app.get("/api/projects/{pid}/items/{iid}/file")
def project_draft_file(request: Request, pid: int, iid: int):
    p = _pj(projects.draft_path, pid, request.state.user, iid)
    return FileResponse(p, filename=p.name.split("__", 1)[-1], media_type="text/plain")


@app.post("/api/projects/{pid}/items/{iid}/promote")
def project_promote(request: Request, pid: int, iid: int, plasmid_id: str = Body(..., embed=True),
                    construct_name: str = Body("", embed=True)):
    """A draft got its lab ID: upload it to the bank (synced to the share) and keep a link here."""
    p = _pj(projects.draft_path, pid, request.state.user, iid)
    r = _store_new_map(request, p.read_bytes(), plasmid_id, construct_name)
    _pj(projects.promote, pid, request.state.user, iid, r["plasmid_id"])
    return r


@app.delete("/api/projects/{pid}/items/{iid}")
def project_remove_item(request: Request, pid: int, iid: int):
    _pj(projects.remove_item, pid, request.state.user, iid)
    return {"ok": True}


@app.get("/api/uploads")
def list_uploads(request: Request, all: bool = False):
    """Uploads and their sync status. `mine` says which ones this user may delete (before syncing)."""
    me = _upload_user(request)
    return [{**r, "mine": r["user"] == me} for r in up.status(UPLOADS, None if all else me)]


@app.delete("/api/uploads/{plasmid_id}")
def delete_upload(request: Request, plasmid_id: str):
    """Drop one of my own queued uploads before the hourly sync copies it to the share."""
    with _BANK_LOCK:
        try:
            up.delete_queued(UPLOADS, _upload_user(request), plasmid_id.upper())
        except PermissionError as e:
            raise HTTPException(403, str(e))
        except ValueError as e:
            raise HTTPException(400, str(e))
        for i, m in enumerate(bank.maps):  # it was added to the in-memory bank at upload time
            if _uploader(m) and m.plasmid_id == plasmid_id.upper():
                bank.maps.pop(i)
                break
        _bank_changed()
        _hits_cache.clear()
    return {"ok": True, "plasmid_id": plasmid_id.upper()}


@app.get("/api/me")
def me(request: Request):
    box = _outbox(request)
    files = sorted(box.glob("*.gb"), key=lambda p: p.stat().st_mtime, reverse=True)
    u = request.state.user
    return {"user": u, "auth": auth.MODE, "role": auth.role_of(u), "must_change": auth.must_change(u),
            "files": [{"name": p.name, "bytes": p.stat().st_size} for p in files[:200]]}


@app.post("/api/account/password")
def change_password(request: Request, current: str = Body(..., embed=True), new: str = Body(..., embed=True)):
    if auth.MODE != "local":
        raise HTTPException(400, "passwords are only managed here when the app signs people in itself")
    try:
        auth.change_own_password(request.state.user, current, new)
    except ValueError as e:
        raise HTTPException(400, str(e))
    auth.audit(request.state.user, "changed own password", request.state.user)
    return {"ok": True}


# ── admin: accounts and the oligo export ──────────────────────────────────────


def _admin(request: Request) -> str:
    if auth.role_of(request.state.user) != "admin":
        raise HTTPException(403, "admins only")
    return request.state.user


@app.get("/admin")
def admin_page():
    return FileResponse(STATIC / "admin.html", headers={"Cache-Control": "no-cache"})


@app.get("/api/admin/users")
def admin_users(request: Request):
    _admin(request)
    return {"auth": auth.MODE, "roles": auth.ROLES,
            "users": [{**u, "last_active": activity.get(u["name"])} for u in auth.list_users()]}


def _account(request: Request, action: str, fn, name: str, *args, detail: str = ""):
    actor = _admin(request)
    try:
        out = fn(name, *args)
    except ValueError as e:
        raise HTTPException(400, str(e))
    auth.audit(actor, action, name, detail)
    return out


@app.post("/api/admin/users")
def admin_add_user(request: Request, name: str = Body(..., embed=True), role: str = Body("user", embed=True)):
    name = name.strip().lower()
    pw = _account(request, "add", auth.create_user, name, role, detail=role)
    return {"name": name, "role": role, "password": pw}  # shown once, never stored in clear


@app.post("/api/admin/users/{name}/reset")
def admin_reset(request: Request, name: str):
    if auth.MODE != "local":
        raise HTTPException(400, "passwords are only managed here when the app signs people in itself")
    return {"name": name, "password": _account(request, "reset password", auth.reset_password, name)}


@app.put("/api/admin/users/{name}/role")
def admin_role(request: Request, name: str, role: str = Body(..., embed=True)):
    _account(request, "role", auth.set_role, name, role, detail=role)
    return {"ok": True}


@app.delete("/api/admin/users/{name}")
def admin_remove(request: Request, name: str):
    if name == request.state.user:
        raise HTTPException(400, "you cannot remove your own account")
    _account(request, "remove", auth.remove_user, name)
    return {"ok": True}


@app.get("/api/admin/audit")
def admin_audit(request: Request, limit: int = 50):
    import json
    _admin(request)
    if not auth.AUDIT_FILE.exists():
        return []
    rows = [json.loads(x) for x in auth.AUDIT_FILE.read_text(encoding="utf-8").splitlines() if x.strip()]
    return rows[::-1][: max(1, min(limit, 500))]


def _oligo_info() -> dict:
    from datetime import datetime
    src = _oligo_source()
    nums = [int(o.oligo_id[3:]) for o in oligos if o.oligo_id.startswith("PM-") and o.oligo_id[3:].isdigit()]
    return {"source": "uploaded" if src == OLIGOS_STATE else "seed", "oligos": len(oligos),
            "highest_pm": f"PM-{max(nums):05d}" if nums else None,
            "updated": datetime.fromtimestamp(src.stat().st_mtime).isoformat(timespec="minutes") if src.exists() else None}


@app.get("/api/admin/sync")
def admin_sync(request: Request):
    """The lab share sync as the last runs left it: uploads (ledger) and workspaces (status file)."""
    import json
    _admin(request)
    led = up._ledger(UPLOADS)
    try:
        ws = json.loads((STATE / "workspaces_status.json").read_text())
    except (OSError, ValueError):
        ws = None
    states = [e.get("status") for e in led["files"].values()]
    ow = _overwrite_requests()[-20:][::-1]  # newest first
    return {"o_reachable": led.get("o_reachable"), "uploads": {k: states.count(k) for k in set(states)},
            "overwrites": [{k: r.get(k) for k in ("id", "rel", "new_name", "user", "requested", "status", "message")} for r in ow],
            "queued": len([p for p in up.queued_ids(UPLOADS).values()
                           if led["files"].get(f"{p.parent.name}/{p.name}", {}).get("status") != "synced"]),
            "workspaces": ws}


@app.delete("/api/admin/overwrites/{rid}")
def admin_delete_overwrite(request: Request, rid: str):
    """The admin's fix button: delete one overwrite request from <state>/overwrites/.
    A queued one is then never written to the share, and the bank shows the map from the index again; a done or
    refused one only leaves the list (the file on the share and its _replaced/ copy are not touched)."""
    actor = _admin(request)
    if not re.fullmatch(r"[\w-]{1,60}", rid):
        raise HTTPException(400, "bad request id")
    meta = OVERWRITES / f"{rid}.json"
    if not meta.is_file():
        raise HTTPException(404, "no such overwrite request")
    req = json.loads(meta.read_text())
    for f in (meta, OVERWRITES / req.get("file", f"{rid}.gb")):
        f.unlink(missing_ok=True)
    if req.get("status") == "queued":
        fresh = _load_bank()  # the overlay is gone: show the map as the index has it
        with _BANK_LOCK:
            bank.maps[:] = fresh.maps
            _bank_changed()
            _hits_cache.clear()
    auth.audit(actor, "deleted overwrite request", req.get("rel", rid), f"{req.get('status')}, asked by {req.get('user')}")
    return {"ok": True, "was": req.get("status")}


@app.get("/api/admin/oligos")
def admin_oligos(request: Request):
    _admin(request)
    return _oligo_info()


@app.post("/api/admin/oligos")
async def admin_upload_oligos(request: Request, file: UploadFile = File(...), apply: bool = Form(False)):
    """Merge a new oligo export into the collection. Without apply=true it is only a preview."""
    from ..formats import decode
    from ..oligos import merge_export
    actor = _admin(request)
    src = _oligo_source()
    current = src.read_text(encoding="utf-8", errors="replace") if src.exists() else ""
    try:
        text, diff = merge_export(current, decode(await _upload_bytes(file))[0])
    except ValueError as e:
        raise HTTPException(400, str(e))
    apply = apply and bool(diff["new"] or diff["changed"] or src != OLIGOS_STATE)  # nothing new: nothing written
    if apply:
        if OLIGOS_STATE.exists():
            OLIGOS_STATE.replace(OLIGOS_STATE.with_suffix(".prev.tab"))  # one step of undo
        tmp = OLIGOS_STATE.with_suffix(".tmp")
        tmp.write_text(text, encoding="utf-8", newline="")
        tmp.replace(OLIGOS_STATE)
        _load_oligos()
        auth.audit(actor, "oligo export", file.filename or "", f"+{len(diff['new'])} new, {len(diff['changed'])} changed")
    return {**diff, "applied": apply, **_oligo_info()}


@app.post("/api/oligos/lookup")
def oligo_lookup(sequences: list[str] = Body(..., embed=True)):
    """For each sequence as ordered (tail and [label] codes included): the oligo IDs that have it.

    Both the bases and the labels must agree, so a plain primer never picks up the PM number of its
    FAM-labelled twin. Case and whitespace are ignored.
    """
    from ..oligos import label_codes, normalise
    by_key: dict[tuple, list[str]] = {}
    for o in oligos:
        if o.seq:
            by_key.setdefault((o.seq, label_codes(o.raw_seq)), []).append(o.oligo_id)
    return {s: by_key.get((normalise(s)[0], label_codes(s)), []) for s in sequences[:500]}


@app.get("/api/me/files/{name}")
def my_file(request: Request, name: str):
    p = _outbox(request) / Path(name).name  # never a path outside the user's own folder
    if not p.is_file():
        raise HTTPException(404, "no such file")
    return FileResponse(p, filename=p.name, media_type="text/plain")


@app.get("/api/maps/{key}/sites")
def sites(key: int, set: str = "unique", enzymes: str = ""):
    m = _get(key)
    if not m.seq:
        raise HTTPException(400, "map has no sequence")
    try:
        return find_sites(m, set, enzymes)
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.post("/api/maps/{key}/verify")
async def verify_clone(
    key: int, files: list[UploadFile] = File(...), mode: str = Form("auto")
):
    m = _get(key)
    if not m.seq:
        raise HTTPException(400, "map has no sequence")
    _check_file_count(len(files))
    parsed = []  # (filename, queries or error)
    for f in files:
        data = await _upload_bytes(f)
        try:
            parsed.append((f.filename, load_queries(f.filename or "upload", data)))
        except Exception as e:
            parsed.append((f.filename, str(e)))
    with _heavy_slot():
        return await run_in_threadpool(_verify_reads, m, parsed, mode)


def _verify_reads(m, parsed, mode):
    out, done = [], []  # done: (row, query, result) for every read that aligned
    for fname, queries in parsed:
        if isinstance(queries, str):
            out.append({"file": fname, "error": queries})
            continue
        for q in queries:
            if len(q.seq) > MAX_READ_BP:
                out.append({"file": fname, "error": f"{q.name}: reads are limited to {MAX_READ_BP:,} bp"})
                continue
            r = verify(m, q, None if mode == "auto" else mode)
            row = {"file": fname, **asdict(r), "verdict": r.verdict, "has_trace": bool(q.peaks)}
            out.append(row)
            done.append((row, q, r))
    # Several reads at once: how many reads cover each change and how many agree on it,
    # plus the chromatogram at each change (when the .ab1 file has a peak table)
    support = read_support([r for _, _, r in done], m.length)
    for (row, q, r), sup in zip(done, support):
        for d, s in zip(row["discrepancies"], sup):
            d.update(s)
            d["trace"] = trace_window(q, r.strand, d["query_pos"]) if q.peaks else None
    return out


@app.post("/api/trace")
async def read_trace(file: UploadFile = File(...), strand: str = Form("+"), index: int = Form(0)):
    """The whole chromatogram of one uploaded read, in map orientation (strand from /verify)."""
    data = await _upload_bytes(file)  # outside the try: its 413 must not become a 400
    try:
        queries = load_queries(file.filename or "upload", data)
    except Exception as e:
        raise HTTPException(400, str(e))
    if not 0 <= index < len(queries):
        raise HTTPException(400, "no such read in that file")
    t = trace_full(queries[index], "-" if strand == "-" else "+")
    if t is None:
        raise HTTPException(400, "this file has no peak table")
    return t


# ── label images for a lab database or label printer ──────────────────────────────────
# Insert from URL with cURL option  -H "Authorization: Bearer <token>"  (see _LABEL_ROUTE above).
def _label(make, text: str, px: int | None) -> Response:
    if not labels.TEXT.match(text):
        raise HTTPException(400, "text must be 1-40 of A-Z a-z 0-9 . _ -")
    return Response(make(text, px), media_type="image/png", headers={"Cache-Control": "private, max-age=86400"})


@app.get("/api/qr/{text}.png")
def label_qr(text: str, px: int | None = None):
    """QR code of exactly `text` (EC level M, 4-module quiet zone), px x px; px clamped to 200-1200."""
    return _label(labels.qr_png, text, px)


@app.get("/api/code128/{text}.png")
def label_code128(text: str, px: int | None = None):
    """Code 128 of `text`, bars only, quiet zones included; px = width (400-1800), height = 2/9 of it."""
    return _label(labels.code128_png, text, px)


def main():
    import uvicorn

    uvicorn.run(
        app,
        host=os.environ.get("PLASMIDBANK_HOST", "127.0.0.1"),
        port=int(os.environ.get("PLASMIDBANK_PORT", 8765)),
        proxy_headers=True,
    )


if __name__ == "__main__":
    main()
