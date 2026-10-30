"""Shared lab projects: a named set of map references plus draft files, private to its members.

Projects live inside the web app (not a network folder), and the project owner invites
members. An item is either
- a **ref**: a stable pointer to a bank map, by PMC/EC ID or, for ID-less files, by display path.
  Bank keys are list positions and change between restarts, so they are never stored;
- a **draft**: a file kept in the project store, visible to members only. Promoting a draft with
  the lab ID it was given sends it through the normal upload pipeline (``uploads.store_upload``,
  hourly add-only sync to the share) and turns the item into a ref: "a link replaces the file".

Drafts never enter ``Bank.maps``: search, stats and duplicate detection would leak them to every
user. Store: ``<state>/projects.sqlite`` and ``<state>/projects/<project id>/``.

The personal project is "My workspace" (stores from older versions called it "My drafts"): drafts there sit
in folders (``items.folder``, e.g. ``RinRK1/constructs``), open in the map viewer, are edited in
place, and are mirrored to the share at ``Plasmid map bank/Workspaces/<user>/<folder>/`` by
``uploads.sync_workspaces``. Shared projects keep no folders and no copy on the share.
"""

from __future__ import annotations

import re
import sqlite3
from datetime import datetime
from pathlib import Path

SCHEMA = """
CREATE TABLE IF NOT EXISTS projects (
  id INTEGER PRIMARY KEY, name TEXT NOT NULL, description TEXT NOT NULL DEFAULT '',
  owner TEXT NOT NULL, created TEXT NOT NULL, personal INTEGER NOT NULL DEFAULT 0);
CREATE TABLE IF NOT EXISTS members (
  project_id INTEGER NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
  user TEXT NOT NULL, added TEXT NOT NULL, PRIMARY KEY (project_id, user));
CREATE TABLE IF NOT EXISTS items (
  id INTEGER PRIMARY KEY,
  project_id INTEGER NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
  kind TEXT NOT NULL CHECK (kind IN ('ref', 'draft')),
  plasmid_id TEXT, path TEXT, draft_file TEXT, label TEXT NOT NULL DEFAULT '',
  added_by TEXT NOT NULL, added TEXT NOT NULL, promoted TEXT, folder TEXT NOT NULL DEFAULT '');
"""
MAX_DRAFT = 5 * 1024 * 1024
WORKSPACE = "My workspace"
WORKSPACE_NOTE = "Your own files and folders, copied to the lab share under Plasmid map bank/Workspaces/."
MAX_FOLDER_DEPTH = 4


class Forbidden(Exception):
    """The user is not allowed to do this (not a member, or not the owner)."""


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _safe_file(name: str) -> str:
    stem = re.sub(r"[^\w.+()\-]+", "_", Path(name).name).strip("._")
    return stem[:120] or "map.gb"


def safe_folder(text: str) -> str:
    """'RinRK1 / constructs' -> 'RinRK1/constructs'. Each part is a plain name (no '..', no hidden
    names; a leading '_' is stripped, so no folder can be the sync's own '_deleted'); at most MAX_FOLDER_DEPTH levels."""
    parts = []
    for raw in (text or "").replace("\\", "/").split("/"):
        part = re.sub(r"[^\w .+()\-]+", "_", raw).strip(" ._")[:60]
        if part:
            parts.append(part)
    if len(parts) > MAX_FOLDER_DEPTH:
        raise ValueError(f"folders go at most {MAX_FOLDER_DEPTH} levels deep")
    return "/".join(parts)


class Projects:
    def __init__(self, state: Path):
        self.root = Path(state) / "projects"
        self.root.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(
            Path(state) / "projects.sqlite", check_same_thread=False
        )
        self.db.row_factory = sqlite3.Row
        self.db.execute("PRAGMA foreign_keys = ON")
        self.db.executescript(SCHEMA)
        cols = {r["name"] for r in self.db.execute("PRAGMA table_info(projects)")}
        if "personal" not in cols:  # a store from before 2026-10-22
            self.db.execute(
                "ALTER TABLE projects ADD COLUMN personal INTEGER NOT NULL DEFAULT 0"
            )
        if "folder" not in {r["name"] for r in self.db.execute("PRAGMA table_info(items)")}:
            self.db.execute("ALTER TABLE items ADD COLUMN folder TEXT NOT NULL DEFAULT ''")
        with self.db:  # "My drafts" became "My workspace" on 2026-10-27
            self.db.execute("UPDATE projects SET name = ?, description = ? WHERE personal = 1 AND name = 'My drafts'",
                            (WORKSPACE, WORKSPACE_NOTE))
        self.db.execute(  # one "My drafts" per user, even when two uploads race to create it
            "CREATE UNIQUE INDEX IF NOT EXISTS one_personal ON projects(owner) WHERE personal = 1"
        )

    # ── access ──────────────────────────────────────────────────────────
    def _project(self, pid: int) -> sqlite3.Row:
        row = self.db.execute("SELECT * FROM projects WHERE id = ?", (pid,)).fetchone()
        if row is None:
            raise KeyError("no such project")
        return row

    def is_member(self, pid: int, user: str) -> bool:
        q = "SELECT 1 FROM members WHERE project_id = ? AND user = ?"
        return self.db.execute(q, (pid, user)).fetchone() is not None

    def require_member(self, pid: int, user: str) -> sqlite3.Row:
        p = self._project(pid)
        if not self.is_member(pid, user):
            raise Forbidden("you are not a member of this project")
        return p

    def require_owner(self, pid: int, user: str) -> sqlite3.Row:
        p = self._project(pid)
        if p["owner"] != user:
            raise Forbidden("only the project owner can do this")
        return p

    # ── projects and members ────────────────────────────────────────────
    def create(self, user: str, name: str, description: str = "") -> int:
        name = name.strip()
        if not name:
            raise ValueError("a project needs a name")
        with self.db:
            cur = self.db.execute(
                "INSERT INTO projects (name, description, owner, created) VALUES (?, ?, ?, ?)",
                (name[:120], description.strip()[:1000], user, _now()),
            )
            self.db.execute(
                "INSERT INTO members VALUES (?, ?, ?)", (cur.lastrowid, user, _now())
            )
        return cur.lastrowid

    def personal(self, user: str) -> int:
        """The user's own "My workspace" project, created on first use: where an
        upload without a lab ID goes. Never in the bank; mirrored to the share's Workspaces/."""
        with self.db:
            self.db.execute(
                "INSERT OR IGNORE INTO projects (name, description, owner, created, personal) "
                "VALUES (?, ?, ?, ?, 1)",
                (
                    WORKSPACE,
                    WORKSPACE_NOTE,
                    user,
                    _now(),
                ),
            )
            pid = self.db.execute(
                "SELECT id FROM projects WHERE owner = ? AND personal = 1", (user,)
            ).fetchone()["id"]
            self.db.execute(
                "INSERT OR IGNORE INTO members VALUES (?, ?, ?)", (pid, user, _now())
            )
        return pid

    def for_user(self, user: str) -> list[dict]:
        rows = self.db.execute(
            """SELECT p.*, (SELECT COUNT(*) FROM items i WHERE i.project_id = p.id) AS n_items,
                      (SELECT COUNT(*) FROM members m2 WHERE m2.project_id = p.id) AS n_members
               FROM projects p JOIN members m ON m.project_id = p.id
               WHERE m.user = ? ORDER BY p.created DESC""",
            (user,),
        ).fetchall()
        return [dict(r) for r in rows]

    def refs_of(self, user: str, plasmid_id: str) -> dict[int, int]:
        """{project id: item id} for the user's projects that link this plasmid ID (the map view's star)."""
        rows = self.db.execute(
            """SELECT i.project_id, i.id FROM items i JOIN members m ON m.project_id = i.project_id
               WHERE m.user = ? AND i.kind = 'ref' AND i.plasmid_id = ?""",
            (user, plasmid_id),
        ).fetchall()
        return {r[0]: r[1] for r in rows}

    def get(self, pid: int, user: str) -> dict:
        p = dict(self.require_member(pid, user))
        p["members"] = [
            r["user"]
            for r in self.db.execute(
                "SELECT user FROM members WHERE project_id = ? ORDER BY added", (pid,)
            )
        ]
        p["items"] = [
            dict(r)
            for r in self.db.execute(
                "SELECT * FROM items WHERE project_id = ? ORDER BY id", (pid,)
            )
        ]
        return p

    def add_member(self, pid: int, owner: str, member: str) -> None:
        self.require_owner(pid, owner)
        with self.db:
            self.db.execute(
                "INSERT OR IGNORE INTO members VALUES (?, ?, ?)", (pid, member, _now())
            )

    def remove_member(self, pid: int, owner: str, member: str) -> None:
        p = self.require_owner(pid, owner)
        if member == p["owner"]:
            raise ValueError("the owner cannot be removed")
        with self.db:
            self.db.execute(
                "DELETE FROM members WHERE project_id = ? AND user = ?", (pid, member)
            )

    def delete(self, pid: int, owner: str) -> None:
        self.require_owner(pid, owner)
        for f in (self.root / str(pid)).glob("*"):
            f.unlink()
        with self.db:
            self.db.execute("DELETE FROM projects WHERE id = ?", (pid,))

    # ── items ───────────────────────────────────────────────────────────
    def add_ref(
        self,
        pid: int,
        user: str,
        plasmid_id: str | None,
        path: str | None,
        label: str = "",
    ) -> int:
        self.require_member(pid, user)
        if not plasmid_id and not path:
            raise ValueError("a reference needs an ID or a file path")
        dup = self.db.execute(
            "SELECT id FROM items WHERE project_id = ? AND kind = 'ref' AND "
            "((? IS NOT NULL AND plasmid_id = ?) OR (? IS NULL AND path = ?))",
            (pid, plasmid_id, plasmid_id, plasmid_id, path),
        ).fetchone()
        if dup:
            return dup["id"]
        with self.db:
            cur = self.db.execute(
                "INSERT INTO items (project_id, kind, plasmid_id, path, label, added_by, added) "
                "VALUES (?, 'ref', ?, ?, ?, ?, ?)",
                (
                    pid,
                    plasmid_id,
                    None if plasmid_id else path,
                    label[:200],
                    user,
                    _now(),
                ),
            )
        return cur.lastrowid

    def add_draft(
        self, pid: int, user: str, filename: str, data: bytes, label: str = "", folder: str = ""
    ) -> int:
        self.require_member(pid, user)
        sub = self._folder_for(pid, folder)
        if len(data) > MAX_DRAFT:
            raise ValueError("file too large for a plasmid map (5 MB max)")
        if not data.strip():
            raise ValueError("empty file")
        folder = self.root / str(pid)
        folder.mkdir(parents=True, exist_ok=True)
        with self.db:
            cur = self.db.execute(
                "INSERT INTO items (project_id, kind, label, added_by, added, folder) "
                "VALUES (?, 'draft', ?, ?, ?, ?)",
                (pid, (label or Path(filename).stem)[:200], user, _now(), sub),
            )
            name = f"{cur.lastrowid}__{_safe_file(filename)}"
            (folder / name).write_bytes(data)
            self.db.execute(
                "UPDATE items SET draft_file = ? WHERE id = ?", (name, cur.lastrowid)
            )
        return cur.lastrowid

    def _folder_for(self, pid: int, folder: str) -> str:
        folder = safe_folder(folder)
        if folder and not self._project(pid)["personal"]:
            raise ValueError("folders are for My workspace; shared projects have none")
        return folder

    def find(self, item_id: int, user: str) -> dict:
        """An item by its id alone (the map viewer's key for a draft), if the user is a member."""
        row = self.db.execute("SELECT project_id FROM items WHERE id = ?", (item_id,)).fetchone()
        if row is None:
            raise KeyError("no such item")
        return self.item(row["project_id"], user, item_id)

    def move(self, pid: int, user: str, item_id: int, folder: str | None = None, label: str | None = None) -> None:
        """Put a draft in another folder ('' = the top) and/or rename it."""
        it = self.item(pid, user, item_id)
        if it["kind"] != "draft":
            raise ValueError("only drafts can be moved or renamed")
        new_folder = it["folder"] if folder is None else self._folder_for(pid, folder)
        new_label = it["label"] if label is None else (label.strip()[:200] or it["label"])
        with self.db:
            self.db.execute("UPDATE items SET folder = ?, label = ? WHERE id = ?", (new_folder, new_label, item_id))

    def rename_folder(self, pid: int, user: str, old: str, new: str) -> int:
        """Rename a My workspace folder, its subfolders with it ('A' -> 'B' makes 'A/x' 'B/x'); moving it
        into another folder is the same call ('A' -> 'Done/A'). Returns the number of files moved."""
        p = self.require_owner(pid, user)
        if not p["personal"]:
            raise ValueError("folders are for My workspace; shared projects have none")
        old, new = safe_folder(old), safe_folder(new)
        if not old or not new:
            raise ValueError("give the folder and its new name")
        if new == old:
            return 0
        if new.startswith(old + "/"):
            raise ValueError("a folder cannot go inside itself")
        rows = self.db.execute(
            "SELECT id, folder FROM items WHERE project_id = ? AND (folder = ? OR folder LIKE ? ESCAPE '\\')",
            (pid, old, old.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_") + "/%"),
        ).fetchall()
        if not rows:
            raise KeyError("no such folder")
        moved = [(safe_folder(new + r["folder"][len(old):]), r["id"]) for r in rows]  # depth checked here
        with self.db:
            self.db.executemany("UPDATE items SET folder = ? WHERE id = ?", moved)
        return len(moved)

    def replace_draft(self, pid: int, user: str, item_id: int, data: bytes) -> None:
        """Save over a draft (My workspace is edited in place; the bank never is)."""
        if len(data) > MAX_DRAFT:
            raise ValueError("file too large for a plasmid map (5 MB max)")
        p = self.draft_path(pid, user, item_id)
        tmp = p.with_name(p.name + ".tmp")
        tmp.write_bytes(data)
        tmp.replace(p)

    def workspace_drafts(self) -> list[dict]:
        """Every draft in every My workspace, for the share sync: owner, folder, label, file."""
        rows = self.db.execute(
            """SELECT i.id, i.folder, i.label, i.draft_file, p.id AS project_id, p.owner
               FROM items i JOIN projects p ON p.id = i.project_id
               WHERE p.personal = 1 AND i.kind = 'draft' AND i.draft_file IS NOT NULL ORDER BY i.id"""
        ).fetchall()
        return [{**dict(r), "path": self.root / str(r["project_id"]) / r["draft_file"]} for r in rows]

    def item(self, pid: int, user: str, item_id: int) -> dict:
        self.require_member(pid, user)
        row = self.db.execute(
            "SELECT * FROM items WHERE id = ? AND project_id = ?", (item_id, pid)
        ).fetchone()
        if row is None:
            raise KeyError("no such item")
        return dict(row)

    def draft_path(self, pid: int, user: str, item_id: int) -> Path:
        it = self.item(pid, user, item_id)
        if it["kind"] != "draft":
            raise ValueError("this item is a link to a bank map, not a draft")
        return self.root / str(pid) / it["draft_file"]

    def promote(self, pid: int, user: str, item_id: int, plasmid_id: str) -> None:
        """The draft was uploaded to the bank under ``plasmid_id``: keep a link, drop the file."""
        p = self.draft_path(pid, user, item_id)
        with self.db:
            self.db.execute(
                "UPDATE items SET kind = 'ref', plasmid_id = ?, draft_file = NULL, promoted = ? "
                "WHERE id = ?",
                (plasmid_id, _now(), item_id),
            )
        p.unlink(missing_ok=True)

    def remove_item(self, pid: int, user: str, item_id: int) -> None:
        it = self.item(pid, user, item_id)
        with self.db:
            self.db.execute("DELETE FROM items WHERE id = ?", (item_id,))
        if it["kind"] == "draft" and it["draft_file"]:
            (self.root / str(pid) / it["draft_file"]).unlink(missing_ok=True)
