"""Where genomes live: one folder per genome under GENOMES_DIR (see build.py for its contents).

The registry is the directory listing itself: a genome exists when `<slug>/genome.sqlite` does.
Folders starting with '.' are working space (.building, .incoming, .jobs).
"""

from __future__ import annotations

import json
import os
import shutil
import sqlite3
from pathlib import Path

from .build import SLUG_RE


def genomes_dir() -> Path:
    d = Path(
        os.environ.get("GENOMES_DIR", Path.home() / "GenomeBrowser" / "genomes")
    ).expanduser()
    d.mkdir(parents=True, exist_ok=True)
    return d


def genome_path(slug: str) -> Path | None:
    if not SLUG_RE.match(slug or ""):
        return None
    p = genomes_dir() / slug
    return p if (p / "genome.sqlite").is_file() else None


def connect(slug: str) -> sqlite3.Connection | None:
    p = genome_path(slug)
    if not p:
        return None
    c = sqlite3.connect(
        f"file:{p / 'genome.sqlite'}?mode=ro", uri=True, check_same_thread=False
    )
    c.row_factory = sqlite3.Row
    return c


def meta(slug: str) -> dict:
    p = genome_path(slug)
    if not p:
        return {}
    try:
        m = json.loads((p / "meta.json").read_text())
    except (OSError, ValueError):
        m = {}
    m["slug"] = slug
    m.setdefault("name", slug)
    return m


def list_genomes() -> list[dict]:
    out = []
    for p in sorted(genomes_dir().iterdir()):
        if p.name.startswith(".") or not (p / "genome.sqlite").is_file():
            continue
        m = meta(p.name)
        try:
            c = connect(p.name)
            m["n_seqs"], m["length"] = c.execute(
                "SELECT count(*), coalesce(sum(length), 0) FROM seq"
            ).fetchone()
            m["tracks"] = [
                dict(r)
                for r in c.execute("SELECT name, kind, n FROM tracks ORDER BY rowid")
            ]
            c.close()
        except sqlite3.Error:
            continue
        m["size"] = sum(f.stat().st_size for f in p.iterdir() if f.is_file())
        out.append(m)
    return out


def delete_genome(slug: str) -> None:
    p = genome_path(slug)
    if not p:
        raise KeyError(slug)
    trash = genomes_dir() / ".building" / f"deleted-{slug}"
    shutil.rmtree(trash, ignore_errors=True)
    trash.parent.mkdir(exist_ok=True)
    p.rename(trash)  # gone from the listing at once; the files are removed after
    shutil.rmtree(trash, ignore_errors=True)


def disk() -> dict:
    u = shutil.disk_usage(genomes_dir())
    return {"total": u.total, "used": u.used, "free": u.free}
