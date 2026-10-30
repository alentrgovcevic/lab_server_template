"""SQLite index of the map bank (v1). Replaces the in-memory `bank.Bank` for anything folder-sized.

One row per file ever seen, keyed by its NFC-normalised absolute path. Rows are never deleted: a file
that disappears gets `missing_since`, and every helper query ignores missing rows unless asked.

Incremental: a file whose size and mtime (ns) match its row is not opened again; the cached parse
is reused. `PARSER` is stored per row, so bumping it forces a re-parse of everything.
"""

from __future__ import annotations

import json
import re
import sqlite3
import unicodedata
from collections import defaultdict
from datetime import datetime
from pathlib import Path, PurePath

from .formats import MapRecord

DEFAULT_HOME = Path.home() / "PlasmidBank"
DEFAULT_DB = DEFAULT_HOME / "index.sqlite"
PARSER = 1  # bump when formats.read_map / ids.identify change what they return

SCHEMA = """
CREATE TABLE IF NOT EXISTS file (
    id            INTEGER PRIMARY KEY,
    path          TEXT NOT NULL UNIQUE,      -- absolute, NFC
    mirror_path   TEXT,                      -- absolute path of the current mirror copy
    kind          TEXT,                      -- PMC | EC | NONE
    plasmid_id    TEXT,
    id_source     TEXT,                      -- filename | locus | none
    format        TEXT,                      -- genbank | fasta | raw | clc | unknown
    status        TEXT,                      -- ok | warn | error | clc_binary
    size          INTEGER,
    mtime         INTEGER,                   -- st_mtime_ns; NULL forces a re-parse
    sha256        TEXT,                      -- of the raw bytes
    content_hash  TEXT,                      -- sequtil.canonical_hash (rotation/strand invariant)
    length        INTEGER,
    topology      TEXT,
    name          TEXT,
    definition    TEXT,
    parser        INTEGER,
    first_seen    TEXT NOT NULL,
    last_seen     TEXT NOT NULL,
    seen_run      INTEGER,                   -- run.id that last saw it (drives missing_since)
    missing_since TEXT,
    messages      TEXT NOT NULL DEFAULT '[]' -- JSON list
);
CREATE INDEX IF NOT EXISTS file_pid ON file(plasmid_id);
CREATE INDEX IF NOT EXISTS file_hash ON file(content_hash);
CREATE TABLE IF NOT EXISTS feature (
    file_id INTEGER NOT NULL REFERENCES file(id),
    type    TEXT,
    label   TEXT,
    start   INTEGER,
    "end"   INTEGER,
    strand  INTEGER
);
CREATE INDEX IF NOT EXISTS feature_file ON feature(file_id);
CREATE TABLE IF NOT EXISTS xref (
    file_id       INTEGER NOT NULL REFERENCES file(id),
    referenced_id TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS xref_file ON xref(file_id);
CREATE INDEX IF NOT EXISTS xref_ref ON xref(referenced_id);
CREATE TABLE IF NOT EXISTS run (
    id       INTEGER PRIMARY KEY,
    started  TEXT NOT NULL,
    finished TEXT,
    root     TEXT NOT NULL,
    counts   TEXT,                           -- JSON
    errors   TEXT                            -- JSON list
);
"""

# Bin folders: "PMC-03001- PMC-03500", "PMC-04001-4500", "PMC-06500-7000", "PMC-00001 - PMC-00500".
BIN_RE = re.compile(
    r"\s*PMC[-_ ]?0*(\d{1,5})\s*[-–—_]+\s*(?:PMC[-_ ]?)?0*(\d{1,5})(?!\d)",
    re.IGNORECASE,
)


def nfc(path: str | PurePath) -> str:
    return unicodedata.normalize("NFC", str(path))


def now() -> str:
    return datetime.now().isoformat(timespec="microseconds")


def parse_bin(folder: str) -> tuple[int, int] | None:
    """(lo, hi) for a PMC bin folder name, matched loosely; None if it is not a bin."""
    m = BIN_RE.match(unicodedata.normalize("NFC", folder))
    if not m:
        return None
    lo, hi = int(m.group(1)), int(m.group(2))
    return (lo, hi) if lo < hi else None


def bin_of(path: str, root: str | None = None) -> tuple[str, int, int] | None:
    """The nearest ancestor folder (below `root`, if given) that names a PMC bin."""
    p = PurePath(path)
    parents = p.parents
    if root:
        try:
            parents = [PurePath(root) / q for q in p.relative_to(root).parents]
        except ValueError:
            pass
    for d in parents:
        b = parse_bin(d.name)
        if b:
            return d.name, b[0], b[1]
    return None


def _under(root: str | None) -> tuple[str, list]:
    """SQL fragment limiting file.path to `root` (substr, not LIKE: paths contain `_` and `%`)."""
    if not root:
        return "1", []
    prefix = nfc(root).rstrip("/") + "/"
    return "substr(file.path, 1, ?) = ?", [len(prefix), prefix]


class Index:
    def __init__(self, path: str | Path = DEFAULT_DB):
        self.path = Path(path).expanduser()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.path)
        self.db.row_factory = sqlite3.Row
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.execute("PRAGMA foreign_keys=ON")
        self.db.executescript(SCHEMA)
        self.run_id: int | None = None  # set by start_run; stamped on every file it sees

    def close(self) -> None:
        self.db.commit()
        self.db.close()

    def commit(self) -> None:
        self.db.commit()

    # ── runs ──────────────────────────────────────────────────────────────────

    def start_run(self, root: str) -> tuple[int, str]:
        started = now()
        cur = self.db.execute(
            "INSERT INTO run(started, root) VALUES (?, ?)", (started, nfc(root))
        )
        self.db.commit()
        self.run_id = cur.lastrowid
        return cur.lastrowid, started

    def finish_run(
        self, run_id: int, counts: dict, errors: list[str], mark_missing: bool
    ) -> int:
        """Close a run. Files under its root that this run did not see get `missing_since`.

        Decided by run id, not by comparing timestamps (two runs can start in the same tick).

        `mark_missing` must be False when the walk was incomplete (an SMB drop mid-walk would
        otherwise flag thousands of present files as gone).
        """
        run = self.db.execute("SELECT * FROM run WHERE id = ?", (run_id,)).fetchone()
        missing = 0
        if mark_missing:
            where, args = _under(run["root"])
            cur = self.db.execute(
                f"UPDATE file SET missing_since = ? WHERE {where} "
                "AND (seen_run IS NULL OR seen_run != ?) AND missing_since IS NULL",
                [run["started"], *args, run_id],
            )
            missing = cur.rowcount
        counts = {**counts, "newly_missing": missing}
        self.db.execute(
            "UPDATE run SET finished = ?, counts = ?, errors = ? WHERE id = ?",
            (now(), json.dumps(counts), json.dumps(errors), run_id),
        )
        self.db.commit()
        return missing

    def last_run(self, root: str | None = None) -> sqlite3.Row | None:
        if root:
            return self.db.execute(
                "SELECT * FROM run WHERE root = ? ORDER BY id DESC LIMIT 1",
                (nfc(root),),
            ).fetchone()
        return self.db.execute("SELECT * FROM run ORDER BY id DESC LIMIT 1").fetchone()

    # ── files ─────────────────────────────────────────────────────────────────

    def get(self, path: str) -> sqlite3.Row | None:
        return self.db.execute(
            "SELECT * FROM file WHERE path = ?", (nfc(path),)
        ).fetchone()

    def cached(self, path: str, size: int, mtime_ns: int) -> sqlite3.Row | None:
        """The row if its parse is still valid for this size+mtime, else None."""
        row = self.get(path)
        if (
            row is not None
            and row["size"] == size
            and row["mtime"] == mtime_ns
            and row["parser"] == PARSER
        ):
            return row
        return None

    def touch(self, file_id: int, seen: str) -> None:
        self.db.execute(
            "UPDATE file SET last_seen = ?, seen_run = ?, missing_since = NULL WHERE id = ?",
            (seen, self.run_id, file_id),
        )

    def set_mirror(self, file_id: int, mirror_path: str | None) -> None:
        self.db.execute(
            "UPDATE file SET mirror_path = ? WHERE id = ?", (mirror_path, file_id)
        )

    def upsert(
        self,
        rec: MapRecord,
        *,
        size: int,
        mtime_ns: int | None,
        sha256: str | None,
        seen: str,
    ) -> int:
        """Insert or replace the parse for `rec.path`; returns the file id. Features/xrefs are rewritten."""
        path = nfc(rec.path)
        vals = dict(
            kind=rec.kind,
            plasmid_id=rec.plasmid_id,
            id_source=rec.id_source,
            format=rec.format,
            status=rec.status,
            size=size,
            mtime=mtime_ns,
            sha256=sha256,
            content_hash=rec.content_hash,
            length=rec.length if rec.seq else None,
            topology=rec.topology,
            name=rec.name,
            definition=rec.description,
            parser=PARSER,
            last_seen=seen,
            seen_run=self.run_id,
            messages=json.dumps(rec.messages),
        )
        row = self.get(path)
        if row is None:
            cols = ["path", "first_seen", *vals]
            cur = self.db.execute(
                f"INSERT INTO file({', '.join(cols)}) VALUES ({', '.join('?' * len(cols))})",
                [path, seen, *vals.values()],
            )
            file_id = cur.lastrowid
        else:
            file_id = row["id"]
            sets = ", ".join(f"{k} = ?" for k in vals)
            self.db.execute(
                f"UPDATE file SET {sets}, missing_since = NULL WHERE id = ?",
                [*vals.values(), file_id],
            )
            self.db.execute("DELETE FROM feature WHERE file_id = ?", (file_id,))
            self.db.execute("DELETE FROM xref WHERE file_id = ?", (file_id,))
        self.db.executemany(
            'INSERT INTO feature(file_id, type, label, start, "end", strand) VALUES (?,?,?,?,?,?)',
            [
                (file_id, f.type, f.label, f.start, f.end, f.strand)
                for f in rec.features
            ],
        )
        self.db.executemany(
            "INSERT INTO xref(file_id, referenced_id) VALUES (?, ?)",
            [(file_id, x) for x in dict.fromkeys(rec.xrefs)],
        )
        return file_id

    # ── queries (present files only unless include_missing) ───────────────────

    def _files(
        self,
        cond: str = "1",
        args: list | None = None,
        root: str | None = None,
        include_missing: bool = False,
        order: str = "path",
    ) -> list[sqlite3.Row]:
        where, rargs = _under(root)
        miss = "1" if include_missing else "missing_since IS NULL"
        return self.db.execute(
            f"SELECT * FROM file WHERE {where} AND {miss} AND ({cond}) ORDER BY {order}",
            [*rargs, *(args or [])],
        ).fetchall()

    def files(self, root: str | None = None, include_missing: bool = False):
        return self._files(root=root, include_missing=include_missing)

    def by_id(
        self, plasmid_id: str, include_missing: bool = False
    ) -> list[sqlite3.Row]:
        return self._files(
            "plasmid_id = ?", [plasmid_id], include_missing=include_missing
        )

    def features(self, file_id: int) -> list[sqlite3.Row]:
        return self.db.execute(
            'SELECT * FROM feature WHERE file_id = ? ORDER BY start, "end"', (file_id,)
        ).fetchall()

    def xrefs(self, file_id: int) -> list[str]:
        return [
            r[0]
            for r in self.db.execute(
                "SELECT referenced_id FROM xref WHERE file_id = ?", (file_id,)
            )
        ]

    def search(
        self, query: str, root: str | None = None, limit: int = 500
    ) -> list[sqlite3.Row]:
        """Every whitespace token must appear in the ID, path, name, definition, an xref or a feature label."""
        toks = query.lower().split()
        conds, args = [], []
        for t in toks:
            esc = t.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
            like = "%" + esc + "%"
            conds.append(
                "(lower(coalesce(plasmid_id,'') || ' ' || path || ' ' || coalesce(name,'') || ' ' "
                "|| coalesce(definition,'')) LIKE ? ESCAPE '\\' "
                "OR EXISTS (SELECT 1 FROM xref x WHERE x.file_id = file.id "
                "AND lower(x.referenced_id) LIKE ? ESCAPE '\\') "
                "OR EXISTS (SELECT 1 FROM feature f WHERE f.file_id = file.id "
                "AND lower(f.label) LIKE ? ESCAPE '\\'))"
            )
            args += [like, like, like]
        rows = self._files(" AND ".join(conds) or "1", args, root=root)
        return rows[:limit]

    def duplicates(self, root: str | None = None) -> list[list[sqlite3.Row]]:
        """Groups of 2+ present files with the same canonical sequence hash."""
        groups: dict[str, list] = defaultdict(list)
        for r in self._files(
            "content_hash IS NOT NULL", root=root, order="content_hash, path"
        ):
            groups[r["content_hash"]].append(r)
        return [g for g in groups.values() if len(g) > 1]

    def ambiguous_ids(self, root: str | None = None) -> dict[str, list[sqlite3.Row]]:
        """IDs held by more than one present file (never synced; includes .clc/.gb pairs)."""
        groups: dict[str, list] = defaultdict(list)
        for r in self._files(
            "plasmid_id IS NOT NULL", root=root, order="plasmid_id, path"
        ):
            groups[r["plasmid_id"]].append(r)
        return {k: v for k, v in groups.items() if len(v) > 1}

    def unidentified(self, root: str | None = None) -> list[sqlite3.Row]:
        return self._files("plasmid_id IS NULL", root=root)

    def from_locus(self, root: str | None = None) -> list[sqlite3.Row]:
        """IDs that came from the LOCUS line, not the filename: flagged for review."""
        return self._files("id_source = 'locus'", root=root)

    def parse_errors(self, root: str | None = None) -> list[sqlite3.Row]:
        return self._files("status = 'error'", root=root)

    def clc_only(self, root: str | None = None) -> list[str]:
        """IDs whose every present file is a CLC binary (nothing readable to sync)."""
        where, args = _under(root)
        return [
            r[0]
            for r in self.db.execute(
                f"SELECT plasmid_id FROM file WHERE {where} AND missing_since IS NULL "
                "AND plasmid_id IS NOT NULL GROUP BY plasmid_id "
                "HAVING SUM(format != 'clc') = 0 ORDER BY plasmid_id",
                args,
            )
        ]

    def missing(self, root: str | None = None) -> list[sqlite3.Row]:
        return self._files("missing_since IS NOT NULL", root=root, include_missing=True)

    def out_of_bin(
        self, root: str | None = None
    ) -> list[tuple[sqlite3.Row, str, int, int]]:
        """PMC files whose number lies outside the range of the bin folder they sit in."""
        out = []
        for r in self._files("kind = 'PMC'", root=root):
            b = bin_of(r["path"], root)
            if b is None:
                continue
            folder, lo, hi = b
            n = int(r["plasmid_id"].split("-")[1])
            if not lo <= n <= hi:
                out.append((r, folder, lo, hi))
        return out

    def counts(self, column: str, root: str | None = None) -> dict[str, int]:
        if column not in {"status", "format", "kind", "id_source"}:
            raise ValueError(column)
        where, args = _under(root)
        return {
            (k if k is not None else "-"): n
            for k, n in self.db.execute(
                f"SELECT {column}, COUNT(*) FROM file WHERE {where} AND missing_since IS NULL "
                f"GROUP BY {column} ORDER BY COUNT(*) DESC",
                args,
            )
        }
