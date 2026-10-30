"""Extra modules beside the plasmid bank, mounted on the main app with `app.include_router(router)`.

What lives here: /ui.js (Basic/Advanced and the colour theme on every page) and the admin's Basic-view
settings. It is also where your own lab modules plug in.

Adding a module (e.g. a safety quiz or a chemical register):
  1. Write the logic in its own file under plasmidbank/ (no FastAPI import there; raise ValueError or
     KeyError and translate them to HTTP errors in the route).
  2. Add its page (`static/<module>.html`) and its routes to this router; build its store in setup(),
     which the server calls once with the state folder.
  3. Turn its dashboard tile on in `static/home.html` (the tiles marked "add your own module"), and add
     an admin card in `static/admin.html` if admins upload or edit its data.
  4. Sign-in is already enforced by the server's middleware for every route here; check roles with
     `_admin(request)` where needed.
"""

from __future__ import annotations

import json
from pathlib import Path

from fastapi import APIRouter, Body, HTTPException, Request
from fastapi.responses import FileResponse

from . import auth

STATIC = Path(__file__).parent / "static"
router = APIRouter()
_state: Path = Path(".")


def setup(state: Path) -> None:
    """Called once by the server with the state folder. Build your modules' stores here."""
    global _state
    _state = state


def _admin(request: Request) -> str:
    if auth.role_of(request.state.user) != "admin":
        raise HTTPException(403, "admins only")
    return request.state.user


def _page(name: str):
    return FileResponse(STATIC / name, headers={"Cache-Control": "no-cache"})


@router.get("/ui.js")
def ui_js():
    return FileResponse(
        STATIC / "ui.js",
        media_type="text/javascript",
        headers={"Cache-Control": "no-cache"},
    )


# ── Basic view: which bank tools show in Basic, set by admins with checkboxes on /admin ──
# Everything not ticked here is Advanced-only. One shared file with a version, like the enzyme list.
VIEW_ITEMS = [
    {"id": "tab:features", "label": "Features tab", "basic": True},
    {"id": "tab:verify", "label": "Verify clone tab", "basic": True},
    {
        "id": "tab:primers",
        "label": "Primers tab (the lab's oligos on the map)",
        "basic": True,
    },
    {"id": "tab:sequence", "label": "Sequence tab (rows) and FASTA tab", "basic": True},
    {"id": "tab:design", "label": "Design primers tab", "basic": False},
    {
        "id": "tab:assemble",
        "label": "Assemble tab (Gibson, document a cloning)",
        "basic": False,
    },
    {"id": "tab:gel", "label": "Gel tab", "basic": False},
    {"id": "tab:align", "label": "Align tab", "basic": False},
    {"id": "tab:clean", "label": "Clean tab", "basic": False},
    {"id": "newmap", "label": "＋ New map button", "basic": True},
    {"id": "gg", "label": "Golden Gate button", "basic": True},
    {"id": "enzymes", "label": "Enzyme finder button", "basic": False},
    {"id": "annotations", "label": "Annotations button", "basic": False},
    {"id": "blast", "label": "BLAST buttons on a map", "basic": False},
]


def _view_file() -> Path:
    return _state / "ui_settings.json"


def _view() -> dict:
    try:
        saved = json.loads(_view_file().read_text())
    except (OSError, ValueError):
        saved = {"version": 0}
    basic = (
        set(saved["basic"])
        if "basic" in saved
        else {i["id"] for i in VIEW_ITEMS if i["basic"]}
    )
    return {
        "version": saved.get("version", 0),
        "items": [{**i, "basic": i["id"] in basic} for i in VIEW_ITEMS],
    }


@router.get("/api/ui/view")
def ui_view():
    return _view()


@router.put("/api/admin/ui/view")
def admin_ui_view(
    request: Request,
    basic: list[str] = Body(..., embed=True),
    version: int = Body(..., embed=True),
):
    actor = _admin(request)
    cur = _view()
    if version != cur["version"]:
        raise HTTPException(409, "the Basic view was changed by someone else; reload")
    known = {i["id"] for i in VIEW_ITEMS}
    keep = sorted(set(basic) & known)
    p = _view_file()
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps({"version": cur["version"] + 1, "basic": keep}))
    tmp.replace(p)
    auth.audit(actor, "basic view", "bank", ", ".join(keep))
    return _view()
