"""In-memory map bank over a folder tree (v0; the SQLite index replaces this in v1)."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .formats import MapRecord, read_map
from .ids import find_ids

# CLC workspace/index files, OS junk and non-sequence documents that live in the map folders.
SKIP_SUFFIXES = {
    ".bucket",
    ".clcinfo",
    ".orderlist2",
    ".fileinfocache",
    ".cfs",
    ".cfe",
    ".si",
    ".fdt",
    ".fdx",
    ".fnm",
    ".frq",
    ".prx",
    ".tii",
    ".tis",
    ".nrm",
    ".del",
    ".gen",
    ".ds_store",
    ".db",
    ".lnk",
    ".zip",
    ".pdf",
    ".doc",
    ".docx",
    ".xls",
    ".xlsx",
    ".txt",
    ".png",
    ".jpg",
    ".pptx",
    # BLAST databases (formatdb/makeblastdb output; e.g. 'Clones_blast_db.db/' in the map bank)
    ".nhr",
    ".nin",
    ".nsq",
    ".nsd",
    ".nsi",
    ".nog",
    ".nal",
    ".ndb",
    ".not",
    ".ntf",
    ".nto",
    ".njs",
    ".phr",
    ".pin",
    ".psq",
    ".psd",
    ".psi",
    ".pog",
    ".pal",
    ".pdb",
    ".pot",
    ".ptf",
    ".pto",
    ".pjs",
    ".log",
}

MAX_BYTES = (
    50 * 1024 * 1024
)  # plasmid maps are KB; anything this big is a genome or archive


@dataclass
class Bank:
    root: Path
    maps: list[MapRecord] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)  # files over MAX_BYTES

    @classmethod
    def load(cls, root: str | Path) -> "Bank":
        root = Path(root)
        bank = cls(root)
        for p in sorted(root.rglob("*")):
            if not p.is_file() or p.name.startswith("."):
                continue
            if {"Workspaces", "_replaced"} & set(p.relative_to(root).parts[:-1]):  # as in scan.EXCLUDE_DIRS
                continue
            if p.suffix.lower() in SKIP_SUFFIXES or p.name.lower().startswith(
                "segments"
            ):
                continue
            if (
                p.stat().st_size > MAX_BYTES
            ):  # e.g. "Morex genome v1.0 data files/" sits in the map bank
                bank.skipped.append(str(p))
                continue
            bank.maps.append(read_map(p))
        return bank

    @classmethod
    def from_index(cls, db_path: str | Path, root: str | None = None) -> "Bank":
        """Current (not missing) files from the v1 index, read from their local mirror copies.

        No file on the share is opened: maps are parsed from the mirror, and CLC binaries
        (never mirrored) become metadata-only records. `display_path` keeps the original path.
        """
        from .index import Index

        idx = Index(db_path)
        bank = cls(Path(root) if root else Path("/"))
        try:
            for row in idx.files(root=root):
                if row["mirror_path"] and Path(row["mirror_path"]).exists():
                    m = read_map(row["mirror_path"])
                    m.plasmid_id, m.kind, m.id_source = (
                        row["plasmid_id"],
                        row["kind"],
                        row["id_source"],
                    )
                    m.xrefs = idx.xrefs(row["id"])
                else:
                    m = MapRecord(
                        row["mirror_path"] or row["path"],
                        row["plasmid_id"],
                        row["kind"] or "NONE",
                        row["id_source"] or "none",
                        idx.xrefs(row["id"]),
                        row["format"] or "unknown",
                        row["status"] or "error",
                        name=row["name"] or "",
                        description=row["definition"] or "",
                        topology=row["topology"],
                        content_hash=row["content_hash"],
                    )
                    if row["format"] == "clc":
                        m.messages.append("CLC binary: needs export to GenBank")
                    else:
                        m.messages.append(
                            "no local mirror copy (file not readable as text)"
                        )
                m.display_path = row["path"]
                bank.maps.append(m)
        finally:
            idx.close()
        return bank

    def by_id(self) -> dict[str, list[int]]:
        out: dict[str, list[int]] = defaultdict(list)
        for i, m in enumerate(self.maps):
            if m.plasmid_id:
                out[m.plasmid_id].append(i)
        return out

    def duplicate_groups(self) -> list[list[int]]:
        """Maps with identical sequence regardless of rotation/strand, under 2+ files."""
        groups: dict[str, list[int]] = defaultdict(list)
        for i, m in enumerate(self.maps):
            if m.content_hash:
                groups[m.content_hash].append(i)
        return [g for g in groups.values() if len(g) > 1]

    def search(self, query: str) -> list[int]:
        q = query.strip().lower()
        if not q:
            return list(range(len(self.maps)))
        hits = []
        for i, m in enumerate(self.maps):
            hay = " ".join(
                [
                    m.plasmid_id or "",
                    Path(m.display_path or m.path).name,
                    m.name,
                    m.description,
                    " ".join(m.xrefs),
                ]
                + [f.label for f in m.features]
            ).lower()
            if all(tok in hay for tok in q.split()):
                hits.append(i)
        # The map itself before the maps that only name it (56 EC files mention
        # PMC-01118, and pickers that show the first 30 never reached PMC-01118's own file)
        wanted = set(find_ids(query))
        return sorted(hits, key=lambda i: self.maps[i].plasmid_id not in wanted)
