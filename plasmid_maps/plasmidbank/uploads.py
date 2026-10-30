"""New maps uploaded on the server, and their add-only hourly copy to the lab share.

    python -m plasmidbank.uploads --uploads DIR --dest "/data/bank-rw" [--write]

Upload (server side, `store_upload`): a new GenBank/FASTA file lands at once in
`<uploads>/<user>/<ID>_<construct name>.<ext>`, so it is in the bank immediately. IDs are
**not allocated here**: PMC and EC numbers both come from the lab's own records (a database, a
spreadsheet), so the user types the one they were given and
the server only validates it. An ID already in the bank is refused: edits to an existing map go
through Features → Edit → Save as new version (the outbox), never to the share directly.

Sync (`sync`, hourly from a systemd timer on the VM): copies every not-yet-synced upload to the share.
Layout (the example scheme; adapt `bin_for` and friends to your own share):
- PMC: 1000-wide bins at the top of the bank, named `PMC7001-8000`. A number with no bin gets its
  bin created on demand; when a file lands in the second-to-last bin of a run, the next five bins
  are created (few folders, but always room ahead).
- EC: flat in `ENSA/`, no bins.
The old, emptied folders (`PMC-07001-7500/`, `ENSA genbank/`) are ignored. The share counts as reachable
only when `ENSA/` exists, so an unmounted mount point is never mistaken for the share.
🛑 ADD-ONLY outside `Workspaces/` by design: never overwrites, renames, moves or deletes anything on the share (creating
the bin folders above is the one other write). A name
that already exists with other content, or an ID already present in the bin, is refused and
reported. Dry run is the default; `--write` copies.

Workspaces (`sync_workspaces`, same timer): every user's My workspace is synced
BOTH WAYS with `Workspaces/<user>/<folder>/` under the bank.
Inside `Workspaces/` only, the sync may overwrite and moving the files the sync tracks. The
rules are in `sync_workspaces`'s docstring: three-way on the ledger's base hash, conflicts keep both,
deletes are honoured one-sidedly and never destroy anything. The scan and the bank skip
`Workspaces/`, so these files never become bank maps.

State lives in `<uploads>/ledger.json` (written only by `sync`; the server reads it for status).
An upload with no ledger entry is queued. Discord is pinged (`PLASMIDBANK_DISCORD_WEBHOOK`) when O
becomes unreachable, when it recovers, and when a file is refused — not every hour.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import urllib.request
from datetime import datetime
from pathlib import Path

from .formats import read_map
from .ids import identify

ID_RE = re.compile(r"^(PMC-\d{5}|EC\d{5})$")


def normalise_id(text: str) -> str:
    """'pmc 7383' / 'PMC-07383' / 'EC-58489' / 'ec58489' → canonical; ValueError otherwise."""
    r = identify(text.strip())
    if not r.plasmid_id or not ID_RE.match(r.plasmid_id):
        raise ValueError(
            f"not a PMC or EC number: {text!r} (use PMC-07383 or EC58489, from your lab's records)"
        )
    return r.plasmid_id


def safe_construct_name(text: str) -> str:
    s = re.sub(r"[^\w.+()\-]+", "_", text.strip()).strip("._")
    return s[:80]


def suggest_filename(plasmid_id: str, construct: str, fmt: str) -> str:
    ext = ".gb" if fmt == "genbank" else ".fasta"
    name = safe_construct_name(construct)
    return f"{plasmid_id}_{name}{ext}" if name else f"{plasmid_id}{ext}"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def queued_ids(uploads: Path) -> dict[str, Path]:
    """Every upload on the server by plasmid ID (synced or not)."""
    out: dict[str, Path] = {}
    if uploads.is_dir():
        for p in sorted(uploads.glob("*/*")):
            if p.is_file() and not p.name.startswith("."):
                pid = identify(p.name).plasmid_id
                if pid:
                    out[pid] = p
    return out


def store_upload(
    uploads: Path,
    user: str,
    data: bytes,
    plasmid_id: str,
    construct: str,
    existing_ids: set[str],
) -> Path:
    """Validate and write one upload; returns its path. ValueError with a user-facing reason."""
    pid = normalise_id(plasmid_id)
    if pid in existing_ids:
        raise ValueError(
            f"{pid} is already in the map bank; to change it, open it and use "
            "Features → Edit features → Save as new version"
        )
    mine = uploads / user
    other = queued_ids(uploads).get(pid)
    if other and other.parent != mine:
        raise ValueError(f"{pid} was already uploaded by {other.parent.name}")
    if (
        other
        and _ledger(uploads)["files"].get(_rel(uploads, other), {}).get("status")
        == "synced"
    ):
        raise ValueError(
            f"{pid} is already on the share; edit it instead of uploading it again"
        )
    mine.mkdir(parents=True, exist_ok=True)
    tmp = mine / f".upload-{os.getpid()}.tmp"
    tmp.write_bytes(data)
    try:
        m = read_map(tmp)
        if m.format not in ("genbank", "fasta") or not m.seq:
            raise ValueError(
                f"not a readable GenBank or FASTA file ({m.format}, {m.status})"
            )
        dest = mine / suggest_filename(pid, construct or m.name, m.format)
        if other and other != dest:
            other.unlink()  # the user's own earlier, still-queued upload of this ID
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)
    return dest


def delete_queued(uploads: Path, user: str, plasmid_id: str) -> Path:
    """Remove a user's own upload that has NOT been copied to the share yet.

    Only the person who uploaded it may delete it, and only before the sync has copied it: once a
    file is on the share it stays there (the share is add-only, and the copy is not ours to remove).
    """
    p = queued_ids(uploads).get(plasmid_id)
    if p is None:
        raise ValueError(f"no upload for {plasmid_id}")
    if p.parent.name != user:
        raise PermissionError(f"{plasmid_id} was uploaded by {p.parent.name}")
    led = _ledger(uploads)
    entry = led["files"].get(_rel(uploads, p), {})
    if entry.get("status") == "synced":
        raise ValueError(f"{plasmid_id} is already on the share; it cannot be removed from here")
    p.unlink()
    led["files"].pop(_rel(uploads, p), None)
    _save_ledger(uploads, led)
    return p


# ── ledger ────────────────────────────────────────────────────────────────


def _rel(uploads: Path, p: Path) -> str:
    return str(p.relative_to(uploads))


def _ledger(uploads: Path) -> dict:
    p = uploads / "ledger.json"
    try:
        d = json.loads(p.read_text())
    except (OSError, ValueError):
        d = {}
    d.setdefault("files", {})
    d.setdefault("o_reachable", None)
    return d


def _save_ledger(uploads: Path, d: dict) -> None:
    p = uploads / "ledger.json"
    p.parent.mkdir(parents=True, exist_ok=True)  # no uploads yet: the folder does not exist
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(d, indent=1, sort_keys=True))
    tmp.replace(p)


def status(uploads: Path, user: str | None = None) -> list[dict]:
    led = _ledger(uploads)
    rows = []
    for pid, p in queued_ids(uploads).items():
        if user is not None and p.parent.name != user:
            continue
        e = led["files"].get(_rel(uploads, p), {})
        rows.append(
            {
                "plasmid_id": pid,
                "file": p.name,
                "user": p.parent.name,
                "uploaded": datetime.fromtimestamp(p.stat().st_mtime).isoformat(
                    timespec="minutes"
                ),
                "status": e.get("status", "queued"),
                "dest": e.get("dest"),
                "error": e.get("error"),
                "synced_at": e.get("at"),
            }
        )
    return sorted(rows, key=lambda r: r["uploaded"], reverse=True)


# ── sync to the share ─────────────────────────────────────────────────────────────


# Example layout: PMC maps in 1000-wide bins named like
# 'PMC7001-8000' at the top of the bank, and every EC map flat in 'ENSA/'. The old 500-wide folders
# ('PMC-07001-7500') and 'ENSA genbank/' still exist, emptied, and are deliberately ignored: only
# names in the new style count as bins, so a file can never land in a dead folder.
PMC_STEP = 1000
BUFFER = 5  # PMC bins created ahead once a file lands in the second-to-last bin
ENSA = "ENSA"
_PMC_BIN = re.compile(
    r"^PMC(\d{4,5})-(\d{4,5})$"
)  # new names only: no hyphen after 'PMC'


def _kind(plasmid_id: str) -> str:
    return "PMC" if plasmid_id.startswith("PMC") else "EC"


def _num(plasmid_id: str) -> int:
    return int(re.sub(r"\D", "", plasmid_id))


def grid_bin(n: int) -> tuple[int, int]:
    """The PMC bin for number n: 1-1000, 1001-2000, …"""
    a = (n - 1) // PMC_STEP * PMC_STEP + 1
    return a, a + PMC_STEP - 1


def bin_name(a: int, b: int) -> str:
    """Named like the bins of 2026-10-07: 'PMC0001-1000', 'PMC7001-8000', 'PMC10001-11000'."""
    return f"PMC{a:04d}-{b}"


def _pmc_bins(dest: Path) -> list[tuple[int, int, Path]]:
    out = []
    if dest.is_dir():
        for d in dest.iterdir():
            m = _PMC_BIN.match(d.name)
            if m and d.is_dir() and int(m[1]) <= int(m[2]):
                out.append((int(m[1]), int(m[2]), d))
    return out


def bin_for(dest: Path, plasmid_id: str) -> Path | None:
    """Where a map goes: ENSA/ for EC; for PMC the narrowest existing bin covering it, or None."""
    if _kind(plasmid_id) == "EC":
        return dest / ENSA
    n, best = _num(plasmid_id), None
    for a, b, d in _pmc_bins(dest):
        if a <= n <= b and (best is None or b - a < best[0]):
            best = (b - a, d)
    return best[1] if best else None


def series(dest: Path) -> list[tuple[int, int, Path]]:
    """PMC bins of the standard size, by start. Odd-sized ones don't count."""
    return sorted(x for x in _pmc_bins(dest) if PMC_STEP - 2 <= x[1] - x[0] <= PMC_STEP)


def second_to_last_in_run(dest: Path, folder: Path) -> int | None:
    """If `folder` is the second-to-last bin of a run of consecutive bins, the run's last number.

    Runs, not the whole series: a stray high bin must not make every older run's bins look
    "second-to-last"."""
    s = series(dest)
    i = next((k for k, x in enumerate(s) if x[2] == folder), None)
    if i is None:
        return None
    end = i
    while end + 1 < len(s) and s[end + 1][0] == s[end][1] + 1:
        end += 1
    return s[end][1] if end == i + 1 else None


def bins_ahead(dest: Path, after: int, count: int = BUFFER) -> list[Path]:
    """The next `count` PMC bins after number `after` (skipping any that exist)."""
    out, n = [], after + 1
    while len(out) < count:
        a, b = grid_bin(n)
        if not any(x <= a <= y for x, y, _ in _pmc_bins(dest)):
            out.append(dest / bin_name(a, b))
        n = b + 1
    return out


def discord(msg: str) -> bool:
    url = os.environ.get("PLASMIDBANK_DISCORD_WEBHOOK", "")
    if not url:
        return False
    mention = os.environ.get("PLASMIDBANK_DISCORD_MENTION", "")  # e.g. <@123456789>
    body = json.dumps({"content": f"{mention} {msg}".strip()}).encode()
    req = urllib.request.Request(
        url, body, {"Content-Type": "application/json", "User-Agent": "plasmidbank"}
    )
    try:
        urllib.request.urlopen(req, timeout=15).close()
        return True
    except OSError as e:
        print(f"discord alert failed: {e}", file=sys.stderr)
        return False


def sync(uploads: Path, dest: Path, write: bool = False, alert=discord) -> dict:
    led = _ledger(uploads)
    todo = [
        p
        for p in queued_ids(uploads).values()
        if led["files"].get(_rel(uploads, p), {}).get("status") != "synced"
        or led["files"][_rel(uploads, p)].get("sha256") != sha256(p)
    ]
    # An unmounted mount point is an empty local folder: require the bank's own layout, or files
    # would be "synced" onto the server's disk.
    # Since 2026-10-07 that is ENSA/; the old 'ENSA genbank/' alone does not count.
    reachable = (dest / ENSA).is_dir() and os.access(dest, os.W_OK)
    out = {
        "queued": len(todo),
        "copied": [],
        "refused": [],
        "already": [],
        "folders": [],
        "reachable": reachable,
        "write": write,
    }
    if write:
        if reachable is not led["o_reachable"]:
            if not reachable:
                alert(
                    f"lab-server: the lab share is not reachable at {dest}; {len(todo)} upload(s) wait in the queue."
                )
            elif led["o_reachable"] is False:
                alert(
                    f"lab-server: the lab share is reachable again; syncing {len(todo)} upload(s)."
                )
            led["o_reachable"] = reachable
        led["last_run"] = datetime.now().isoformat(timespec="seconds")
    if not reachable:
        if write:
            _save_ledger(uploads, led)
        return out
    for p in todo:
        pid = identify(p.name).plasmid_id
        kind = _kind(pid)
        folder = bin_for(dest, pid)
        if (
            folder is None
        ):  # PMC with no bin covering it: the convention's bin, created on demand
            folder = dest / bin_name(*grid_bin(_num(pid)))
        target = folder / p.name
        entry = {
            "sha256": sha256(p),
            "dest": str(target),
            "at": datetime.now().isoformat(timespec="seconds"),
        }
        clash = (
            [
                q.name
                for q in folder.iterdir()
                if q.is_file()
                and identify(q.name).plasmid_id == pid
                and q.name != p.name
            ]
            if folder.is_dir()
            else []
        )
        if target.exists():
            if sha256(target) == entry["sha256"]:
                out["already"].append(str(target))
                entry["status"] = "synced"
            else:
                entry |= {
                    "status": "refused",
                    "error": f"{target.name} already exists on the share with other content",
                }
        elif clash:
            entry |= {
                "status": "refused",
                "error": f"{pid} is already on the share as {clash[0]}",
            }
        else:
            if not folder.is_dir() and str(folder) not in out["folders"]:
                out["folders"].append(str(folder))
            if write:
                folder.mkdir(
                    exist_ok=True
                )  # never replaces anything: an existing folder is reused
                with open(target, "xb") as fh:  # x: fail rather than overwrite
                    fh.write(p.read_bytes())
                if sha256(target) != entry["sha256"]:
                    entry |= {
                        "status": "error",
                        "error": "copy on the share does not match the upload",
                    }
                else:
                    entry["status"] = "synced"
            out["copied"].append(str(target))
            # EC maps go flat into ENSA/: no bins, so nothing to create ahead
            run_end = second_to_last_in_run(dest, folder) if kind == "PMC" else None
            if (
                run_end is not None
            ):  # landed in the second-to-last bin of its run: add the next five
                ahead = [
                    d for d in bins_ahead(dest, run_end) if str(d) not in out["folders"]
                ]
                out["folders"] += [str(d) for d in ahead]
                if write and ahead:
                    for d in ahead:
                        d.mkdir(exist_ok=True)
                    alert(
                        f"lab-server: {pid} went into {folder.name}, the second-to-last bin, so created "
                        f"{len(ahead)} new bins on the share: {', '.join(d.name for d in ahead)}"
                    )
        if entry.get("status") in ("refused", "error"):
            out["refused"].append(f"{p.parent.name}/{p.name}: {entry['error']}")
            prev = led["files"].get(_rel(uploads, p), {})
            if write and prev.get("error") != entry["error"]:
                alert(
                    f"lab-server: upload {p.parent.name}/{p.name} was not copied to the share: {entry['error']}"
                )
        if write:
            led["files"][_rel(uploads, p)] = entry
    if write:
        _save_ledger(uploads, led)
    return out


WORKSPACES = "Workspaces"
DELETED = "_deleted"


def _ws_name(text: str) -> str:
    return re.sub(r"[^\w .+()\-]+", "_", text).strip(" ._")[:120] or "map"


def sync_workspaces(state: Path, dest: Path, write: bool = False, alert=discord, stable_seconds: int = 120,
                    now: float | None = None) -> dict:
    """Two-way sync of every My workspace with `<dest>/Workspaces/<user>/`.

    A three-way comparison against the ledger's base (the hash both sides had at the last sync):
    only the app changed → push; only the share changed → pull into the app; both → conflict, both kept (the
    share version becomes a new file "… (changed on the share <date>)" in the app, and the app's version goes to
    the share). A deletion is honoured only when the other side did not change since the base: a file
    deleted in the app goes to `Workspaces/<user>/_deleted/` on the share; one deleted on the share leaves the app,
    but a copy is kept in `<state>/workspaces_trash/`. A new readable file (GenBank/FASTA, ≤ 5 MB,
    unchanged for `stable_seconds`) put into `Workspaces/<user>/…` on the share is imported into that
    user's My workspace, in the same folder. Unreadable files (e.g. `.clc`) are listed, never
    touched. Only one run at a time (`<state>/workspaces.lock`); each run writes
    `<state>/workspaces_status.json`.
    """
    import fcntl
    import shutil
    import time as _time

    from .projects import MAX_DRAFT, Projects

    now = _time.time() if now is None else now
    out = {"written": [], "moved": [], "deleted": [], "pulled": [], "imported": [], "removed_in_app": [],
           "conflicts": [], "unreadable": [], "waiting": [], "errors": [], "reachable": False, "busy": False}
    if not (state / "projects.sqlite").is_file() and not (state / "users.json").is_file():
        return out
    out["reachable"] = reachable = (dest / ENSA).is_dir() and os.access(dest, os.W_OK)
    if not reachable:
        return out
    lock = open(state / "workspaces.lock", "a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        out["busy"] = True
        lock.close()
        return out
    try:
        return _sync_workspaces_locked(state, dest, write, alert, stable_seconds, now, out, Projects, MAX_DRAFT,
                                       shutil)
    finally:
        fcntl.flock(lock, fcntl.LOCK_UN)
        lock.close()


def _sync_workspaces_locked(state, dest, write, alert, stable_seconds, now, out, Projects, MAX_DRAFT, shutil):
    led_path = state / "workspaces_ledger.json"
    try:
        led = json.loads(led_path.read_text())
    except (OSError, ValueError):
        led = {}
    if "items" not in led:  # the one-way ledger of 2026-10-27 afternoon: {item id: {rel, sha256}}
        led = {"items": {k: v for k, v in led.items() if isinstance(v, dict) and "rel" in v}, "unreadable": {}}
    items, unreadable = led["items"], led.setdefault("unreadable", {})
    root, stamp = dest / WORKSPACES, datetime.now().strftime("%Y%m%d-%H%M%S")
    pj = Projects(state)
    drafts = {str(d["id"]): d for d in pj.workspace_drafts() if d["path"].is_file()}

    # Who owns which folder: everyone with an account or a workspace. Two names that sanitise to one
    # folder are refused rather than guessed.
    try:
        users = set(json.loads((state / "users.json").read_text()))
    except (OSError, ValueError, TypeError):
        users = set()
    users |= {d["owner"] for d in drafts.values()}
    owner_of: dict[str, str] = {}
    clash: set[str] = set()
    for u in sorted(users):
        f = _ws_name(u)
        if f in owner_of and owner_of[f] != u:
            clash.add(f)
        owner_of.setdefault(f, u)
    for f in clash:
        out["errors"].append(f"folder Workspaces/{f}/ matches more than one user; skipped")
        owner_of.pop(f)

    def stable(p: Path) -> bool:
        return now - p.stat().st_mtime >= stable_seconds

    def readable(p: Path) -> bool:
        try:
            m = read_map(p)
            return bool(m.seq) and m.format in ("genbank", "fasta") and p.stat().st_size <= MAX_DRAFT
        except Exception:
            return False

    def aside(rel: str, why: str, kind: str = "deleted") -> None:
        src = root / rel
        tgt = root / rel.split("/", 1)[0] / DELETED / f"{stamp}_{src.stem}{why}{src.suffix}"
        out[kind].append(f"{rel} -> {tgt.relative_to(root)}")
        if write:
            tgt.parent.mkdir(parents=True, exist_ok=True)
            os.replace(src, tgt)

    def push(key: str, d: dict, rel: str, digest: str) -> None:
        out["written"].append(rel)
        if write:
            t = root / rel
            t.parent.mkdir(parents=True, exist_ok=True)
            tmp = t.with_name(t.name + ".part")
            tmp.write_bytes(d["path"].read_bytes())
            os.replace(tmp, t)
            items[key] = {"rel": rel, "sha256": digest}

    def import_file(owner: str, rel: str, label: str) -> None:
        f = root / rel
        folder = "/".join(rel.split("/")[1:-1])
        out["imported"].append(rel)
        if write:
            pid = pj.personal(owner)
            iid = pj.add_draft(pid, owner, f.name, f.read_bytes(), label, folder)
            items[str(iid)] = {"rel": rel, "sha256": sha256(f)}

    # Where each app draft belongs on the share; tracked paths are reserved first so a new name never lands on them
    taken = {e["rel"] for e in items.values()}
    want: dict[str, str] = {}
    for key, d in drafts.items():
        folder = "/".join(x for x in (_ws_name(d["owner"]), d["folder"]) if x)
        ext = d["path"].suffix or ".gb"
        stem = _ws_name(d["label"] or d["draft_file"].split("__", 1)[-1].rsplit(".", 1)[0])
        rel = f"{folder}/{stem}{ext}"
        cur = items.get(key, {}).get("rel")
        if rel != cur and (rel in taken or rel in want.values() or (root / rel).exists()):
            rel = f"{folder}/{stem}_{key}{ext}"  # a name clash, or a file the sync did not write
        want[key] = cur if cur == rel else rel

    # 1. files the sync already knows: three-way on the base hash
    for key, e in list(items.items()):
        d, f, base = drafts.get(key), root / e["rel"], e["sha256"]
        o_sha = sha256(f) if f.is_file() else None
        if d is None:  # removed in the app
            if o_sha == base:
                aside(e["rel"], "")
            if write:
                del items[key]  # a copy on the share changed since is untracked now, and step 3 imports it
            continue
        a_sha = sha256(d["path"])
        if o_sha is None:  # removed (or moved away) on the share
            if a_sha == base:
                out["removed_in_app"].append(f"{d['owner']}: {d['label']}")
                if write:
                    trash = state / "workspaces_trash"
                    trash.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(d["path"], trash / f"{stamp}_{_ws_name(d['owner'])}_{d['path'].name}")
                    pj.remove_item(d["project_id"], d["owner"], d["id"])
                    del items[key]
                drafts.pop(key)  # gone from the app now: step 2 must not push it back
                continue
            push(key, d, want[key], a_sha)  # edited in the app after the delete on the share: the edit wins
            continue
        o_changed, a_changed = o_sha != base, a_sha != base
        if o_changed and not stable(f):
            out["waiting"].append(e["rel"])  # someone may still be saving it
            continue
        if o_changed and a_changed and o_sha != a_sha:  # both changed: keep both
            out["conflicts"].append(e["rel"])
            if readable(f):
                import_file(d["owner"], e["rel"], f"{d['label']} (changed on the share {stamp[:8]})")
                if write:  # the import must not claim the share path: the app's version goes there
                    items.pop(next(k for k, v in items.items() if v["rel"] == e["rel"] and k != key), None)
            else:
                aside(e["rel"], "_changed-on-share", "conflicts")
            push(key, d, e["rel"], a_sha)
        elif o_changed and not a_changed:  # only the share changed: pull
            if readable(f):
                out["pulled"].append(e["rel"])
                if write:
                    pj.replace_draft(d["project_id"], d["owner"], d["id"], f.read_bytes())
                    items[key] = {"rel": e["rel"], "sha256": o_sha}
            else:
                out["unreadable"].append(e["rel"])
                aside(e["rel"], "_changed-on-share", "conflicts")
                push(key, d, e["rel"], a_sha)
        elif a_changed and not o_changed:
            push(key, d, e["rel"], a_sha)
        elif o_sha == a_sha != base and write:  # the same edit on both sides
            items[key] = {"rel": e["rel"], "sha256": a_sha}
        # moved or renamed in the app
        e = items.get(key, e)
        if key in items and e["rel"] != want[key] and (root / e["rel"]).is_file() and not (root / want[key]).exists():
            out["moved"].append(f"{e['rel']} -> {want[key]}")
            if write:
                (root / want[key]).parent.mkdir(parents=True, exist_ok=True)
                os.replace(root / e["rel"], root / want[key])
                items[key] = {**e, "rel": want[key]}

    # 2. new in the app
    for key, d in drafts.items():
        if key not in items:
            if (root / want[key]).exists():
                continue  # step 1 left a file there this run (dry run); next run sorts it out
            push(key, d, want[key], sha256(d["path"]))

    # 3. new on the share: every user folder exists, and readable files nobody tracks are imported
    tracked = {e["rel"] for e in items.values()} | set(out["written"])
    for folder, owner in sorted(owner_of.items()):
        udir = root / folder
        if write:
            udir.mkdir(parents=True, exist_ok=True)
        if not udir.is_dir():
            continue
        for dirpath, dirnames, filenames in os.walk(udir):
            dirnames[:] = sorted(x for x in dirnames if not x.startswith(".") and x != DELETED)
            for fn in sorted(filenames):
                f = Path(dirpath) / fn
                rel = str(f.relative_to(root))
                if fn.startswith(".") or fn.endswith(".part") or rel in tracked:
                    continue
                if len(rel.split("/")) - 2 > 4:
                    continue  # deeper than My workspace allows
                if not stable(f):
                    out["waiting"].append(rel)
                    continue
                if not readable(f):
                    if unreadable.get(rel) != f.stat().st_size:
                        out["unreadable"].append(rel)
                        unreadable[rel] = f.stat().st_size
                    continue
                unreadable.pop(rel, None)
                import_file(owner, rel, f.stem)

    if write:
        tmp = led_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(led, indent=1, sort_keys=True))
        tmp.replace(led_path)
        status = {"finished": datetime.now().isoformat(timespec="seconds"),
                  **{k: (v if isinstance(v, bool) else len(v)) for k, v in out.items()},
                  "conflict_files": out["conflicts"][:50], "unreadable_files": out["unreadable"][:50],
                  "error_list": out["errors"][:50]}
        tmp = state / "workspaces_status.tmp"
        tmp.write_text(json.dumps(status, indent=1))
        tmp.replace(state / "workspaces_status.json")
        if out["conflicts"] or out["errors"] or out["unreadable"]:
            alert("lab-server workspaces: " + "; ".join(
                x for x in (f"{len(out['conflicts'])} conflict(s), both copies kept" if out["conflicts"] else "",
                            f"{len(out['unreadable'])} file(s) the app cannot read (e.g. .clc)" if out["unreadable"] else "",
                            "; ".join(out["errors"][:3])) if x))
    return out


REPLACED = "_replaced"


def apply_overwrites(state: Path, dest: Path, write: bool = False) -> dict:
    """Carry out the overwrites people asked for in the app (any bank map, after
    confirmations; the old file is kept). The app (read-only on the share) queues each one in
    `<state>/overwrites/<id>.json` + `<id>.gb`; this run, which can write, does it:
    - the file on the share must still be the one the user saw (its sha256 then): otherwise refused, so an edit
      someone made on the share in the meantime is never lost;
    - the old file moves to `_replaced/<date>_<user>_<name>` beside it, and so does a file already at
      the new name (a .clc becomes `<stem>.gbk`, so its exported .gbk twin is replaced too);
    - the new file is written under the old name (a .clc as `.gbk`).
    Only inside the bank, never inside `Workspaces/` (which has its own sync), never `..`."""
    out = {"done": [], "refused": [], "reachable": False}
    q = state / "overwrites"
    if not q.is_dir():
        return out
    out["reachable"] = reachable = (dest / ENSA).is_dir() and os.access(dest, os.W_OK)
    if not reachable:
        return out
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    for meta_path in sorted(q.glob("*.json")):
        try:
            req = json.loads(meta_path.read_text())
        except (OSError, ValueError):
            continue
        if req.get("status") != "queued":
            continue
        rel = Path(req["rel"])

        def finish(status: str, msg: str = "") -> None:
            out["done" if status == "done" else "refused"].append(f"{req['rel']}: {msg or status}")
            if write:
                req.update(status=status, message=msg, finished=datetime.now().isoformat(timespec="seconds"))
                tmp = meta_path.with_suffix(".tmp")
                tmp.write_text(json.dumps(req, indent=1))
                tmp.replace(meta_path)

        if rel.is_absolute() or ".." in rel.parts or not rel.parts or rel.parts[0] == WORKSPACES or REPLACED in rel.parts:
            finish("refused", "not a path inside the bank")
            continue
        target = dest / rel
        if not target.is_file():
            finish("refused", "the file is no longer on the share")
            continue
        if sha256(target) != req["sha_before"]:
            finish("refused", "someone changed the file on the share after it was opened; nothing was overwritten")
            continue
        new = target.with_name(req["new_name"])
        if new.parent != target.parent or REPLACED in new.parts:
            finish("refused", "bad new name")
            continue
        if not (q / req["file"]).is_file():  # an admin deleted the request while this run was going
            continue
        if write:
            keep = target.parent / REPLACED
            keep.mkdir(exist_ok=True)
            who = re.sub(r"[^\w.-]+", "_", req["user"])[:40]
            for old in {target, new}:
                if old.is_file():
                    os.replace(old, keep / f"{stamp}_{who}_{old.name}")
            tmp = new.with_name(new.name + ".part")
            tmp.write_bytes((q / req["file"]).read_bytes())
            os.replace(tmp, new)
        finish("done", f"-> {new.name}, old file in {REPLACED}/")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--uploads", required=True, type=Path)
    ap.add_argument(
        "--dest",
        required=True,
        type=Path,
        help="the 'Plasmid map bank' folder on the share (writable mount)",
    )
    ap.add_argument("--write", action="store_true", help="copy (default: dry run)")
    ap.add_argument("--state", type=Path, help="the app's state folder, for My workspace (default: the parent of --uploads)")
    a = ap.parse_args(argv)
    r = sync(a.uploads, a.dest, write=a.write)
    try:
        o = apply_overwrites(a.state or a.uploads.parent, a.dest, write=a.write)
        print(f"overwrites ({'done' if a.write else 'dry run'}): done {len(o['done'])} · refused {len(o['refused'])}")
        for k in ("done", "refused"):
            for x in o[k]:
                print(f"  {k}: {x}")
        if o["refused"] and a.write:
            discord(f"lab-server: {len(o['refused'])} overwrite(s) on the share refused: " + "; ".join(o["refused"][:3]))
    except (OSError, ValueError) as e:
        print(f"overwrites: ERROR {e}")
    try:  # a workspace problem must not fail the uploads' hourly run
        w = sync_workspaces(a.state or a.uploads.parent, a.dest, write=a.write)
    except (OSError, ValueError) as e:
        print(f"workspaces: ERROR {e}")
        w = {}
    verb = "copied" if a.write else "would copy"
    print(
        f"share reachable: {r['reachable']} · queued {r['queued']} · {verb} {len(r['copied'])} · "
        f"already there {len(r['already'])} · refused {len(r['refused'])}"
    )
    for c in r["copied"]:
        print(f"  {verb}: {c}")
    for x in r["refused"]:
        print(f"  refused: {x}")
    for d in r["folders"]:
        print(f"  {'created' if a.write else 'would create'} bin folder: {d}")
    if w.get("busy"):
        print("workspaces: another run holds the lock; skipped")
    elif w:
        keys = ("written", "pulled", "imported", "moved", "deleted", "removed_in_app", "conflicts", "unreadable",
                "waiting", "errors")
        print(f"workspaces ({'done' if a.write else 'dry run'}): " + " · ".join(f"{k.replace('_', ' ')} {len(w[k])}" for k in keys))
        for k in keys:
            for x in w[k]:
                print(f"  {k}: {x}")
    return 0 if r["reachable"] and not r["refused"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
