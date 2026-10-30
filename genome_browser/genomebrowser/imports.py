"""Getting genomes in: chunked uploads, files already on the server, and the build queue.

Uploads: the page cuts a file into chunks (CHUNK bytes, under the proxy's request limit) and sends
them in order; the server appends each one only if it starts where the file ends, so a broken
upload resumes from `received`. State: `.incoming/<id>.part` + `<id>.json`.

Server files: an admin may import from IMPORT_DIR (a read-only mount, e.g. the lab share) without
uploading; paths are resolved and must stay inside it.

Jobs: `.jobs/<id>.json`, run one at a time by a background thread, each in its own process
(`python -m genomebrowser.build --job`), so a failed or heavy build never takes the server down.
"""

from __future__ import annotations

import json
import os
import secrets
import shutil
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

from .store import genomes_dir

CHUNK = 16 * 1024 * 1024
MAX_CHUNK = 32 * 1024 * 1024
ID_LEN = 16
SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".gff", ".gff3", ".bed", ".tsv", ".txt")


class ImportError_(Exception):
    """Refused; the message is shown to the user."""


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _incoming() -> Path:
    d = genomes_dir() / ".incoming"
    d.mkdir(exist_ok=True)
    return d


def _jobs() -> Path:
    d = genomes_dir() / ".jobs"
    d.mkdir(exist_ok=True)
    return d


def max_upload() -> int:
    return int(float(os.environ.get("GENOMES_MAX_UPLOAD_GB", "50")) * 1e9)


def import_dir() -> Path | None:
    v = os.environ.get("GENOMES_IMPORT_DIR", "")
    return Path(v).resolve() if v and Path(v).is_dir() else None


def _valid_id(i: str) -> bool:
    return len(i) == ID_LEN * 2 and all(c in "0123456789abcdef" for c in i)


# ── uploads ───────────────────────────────────────────────────────────────────


def create_upload(user: str, filename: str, size: int) -> dict:
    name = Path(filename or "").name
    base = name[:-3] if name.lower().endswith(".gz") else name
    if not base.lower().endswith(SUFFIXES):
        raise ImportError_(
            f"{name!r}: expected FASTA, GFF3, BED or a descriptions TSV (optionally .gz)"
        )
    if size <= 0 or size > max_upload():
        raise ImportError_(
            f"file size must be between 1 byte and {max_upload() / 1e9:.0f} GB"
        )
    need = size * (
        5 if name.lower().endswith(".gz") else 2.5
    )  # the upload + the built copy
    free = shutil.disk_usage(_incoming()).free
    if free < need:
        raise ImportError_(
            f"not enough disk space on the server: {free / 1e9:.1f} GB free, "
            f"about {need / 1e9:.1f} GB needed. Ask an admin to clean up or extend the disk."
        )
    uid = secrets.token_hex(ID_LEN)
    st = {
        "id": uid,
        "user": user,
        "filename": name,
        "size": size,
        "received": 0,
        "created": _now(),
    }
    (_incoming() / f"{uid}.part").write_bytes(b"")
    (_incoming() / f"{uid}.json").write_text(json.dumps(st))
    return st | {"chunk": CHUNK}


def upload_status(uid: str, user: str | None = None) -> dict:
    if not _valid_id(uid):
        raise KeyError(uid)
    p = _incoming() / f"{uid}.json"
    if not p.exists():
        raise KeyError(uid)
    st = json.loads(p.read_text())
    if user is not None and st["user"] != user:
        raise KeyError(uid)
    return st


_upload_lock = threading.Lock()


def append_chunk(uid: str, user: str, offset: int, data: bytes) -> dict:
    if len(data) > MAX_CHUNK:
        raise ImportError_(f"chunk too large (max {MAX_CHUNK // 2**20} MB)")
    with _upload_lock:
        st = upload_status(uid, user)
        if offset != st["received"]:
            raise ImportError_(
                f"chunk starts at {offset}, the server has {st['received']} bytes; resume from there"
            )
        if st["received"] + len(data) > st["size"]:
            raise ImportError_("more data than the announced file size")
        with open(_incoming() / f"{uid}.part", "ab") as fh:
            fh.write(data)
        st["received"] += len(data)
        (_incoming() / f"{uid}.json").write_text(json.dumps(st))
    return st


def upload_file(uid: str, user: str) -> Path:
    """The finished upload's path, renamed so the build sees its real extension."""
    st = upload_status(uid, user)
    if st["received"] != st["size"]:
        raise ImportError_(
            f"{st['filename']}: upload incomplete ({st['received']} of {st['size']} bytes)"
        )
    part = _incoming() / f"{uid}.part"
    named = _incoming() / f"{uid}__{st['filename']}"
    if part.exists():
        part.rename(named)
    return named


def discard_upload(uid: str) -> None:
    if not _valid_id(uid):
        return
    for p in _incoming().glob(f"{uid}*"):
        p.unlink(missing_ok=True)


def clean_stale_uploads(days: float = 2) -> int:
    n, cutoff = 0, time.time() - days * 86400
    for p in _incoming().glob("*.json"):
        if p.stat().st_mtime < cutoff:
            discard_upload(p.stem)
            n += 1
    return n


# ── files already on the server ──────────────────────────────────────────────


def server_files(limit: int = 2000) -> list[dict]:
    root = import_dir()
    if not root:
        return []
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        for f in sorted(filenames):
            base = f[:-3] if f.lower().endswith(".gz") else f
            if base.lower().endswith(SUFFIXES):
                p = Path(dirpath) / f
                out.append({"path": str(p.relative_to(root)), "size": p.stat().st_size})
                if len(out) >= limit:
                    return out
    return out


def server_file(rel: str) -> Path:
    root = import_dir()
    if not root:
        raise ImportError_("no server import folder is configured (GENOMES_IMPORT_DIR)")
    p = (root / rel).resolve()
    if root not in p.parents or not p.is_file():
        raise ImportError_(f"{rel!r} is not a file in the import folder")
    return p


# ── jobs ──────────────────────────────────────────────────────────────────────


def submit(job: dict) -> dict:
    jid = secrets.token_hex(ID_LEN)
    job = job | {"id": jid, "status": "queued", "created": _now()}
    (_jobs() / f"{jid}.json").write_text(json.dumps(job, indent=1))
    _wake.set()
    return job


def get_job(jid: str) -> dict:
    if not _valid_id(jid) or not (_jobs() / f"{jid}.json").exists():
        raise KeyError(jid)
    job = json.loads((_jobs() / f"{jid}.json").read_text())
    logp = _jobs() / f"{jid}.log"
    job["log"] = logp.read_text(errors="replace")[-4000:] if logp.exists() else ""
    return job


def list_jobs(n: int = 30) -> list[dict]:
    files = sorted(
        _jobs().glob("*.json"), key=lambda p: p.stat().st_mtime, reverse=True
    )[:n]
    out = []
    for p in files:
        try:
            j = json.loads(p.read_text())
        except ValueError:
            continue
        out.append(
            {
                k: j.get(k)
                for k in (
                    "id",
                    "kind",
                    "name",
                    "slug",
                    "track",
                    "user",
                    "status",
                    "created",
                    "finished",
                    "result",
                    "error",
                )
            }
        )
    return out


def _save(job: dict) -> None:
    (_jobs() / f"{job['id']}.json").write_text(json.dumps(job, indent=1))


def _run(job: dict) -> None:
    job |= {"status": "running", "started": _now()}
    _save(job)
    logp = _jobs() / f"{job['id']}.log"
    with open(logp, "w") as log:
        rc = subprocess.call(
            [
                sys.executable,
                "-m",
                "genomebrowser.build",
                "--genomes-dir",
                str(genomes_dir()),
                "--job",
                str(_jobs() / f"{job['id']}.json"),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=Path(__file__).resolve().parents[1],
        )
    text = logp.read_text(errors="replace").strip().splitlines()
    if rc == 0:
        job |= {"status": "done", "result": text[-1] if text else "done"}
    else:
        err = next((ln[7:] for ln in reversed(text) if ln.startswith("ERROR: ")), None)
        job |= {
            "status": "failed",
            "error": err or (text[-1] if text else f"exit code {rc}"),
        }
    job["finished"] = _now()
    for uid in job.get(
        "uploads", []
    ):  # inputs are consumed either way; a retry re-uploads
        discard_upload(uid)
    _save(job)


_wake = threading.Event()
_started = False


def _worker() -> None:
    while True:
        queued = []
        for p in _jobs().glob("*.json"):
            try:
                j = json.loads(p.read_text())
            except ValueError:
                continue
            if j.get("status") == "queued":
                queued.append(j)
        if queued:
            _run(min(queued, key=lambda j: j["created"]))
            continue
        if not _wake.wait(30):  # idle: drop uploads abandoned for days, or they fill the disk
            clean_stale_uploads()
        _wake.clear()


def start_worker() -> None:
    """Once per process: fail jobs a restart interrupted, then run the queue in the background."""
    global _started
    if _started:
        return
    _started = True
    for p in _jobs().glob("*.json"):
        try:
            j = json.loads(p.read_text())
        except ValueError:
            continue
        if j.get("status") == "running":
            j |= {
                "status": "failed",
                "error": "interrupted by a server restart; submit it again",
                "finished": _now(),
            }
            _save(j)
    shutil.rmtree(genomes_dir() / ".building", ignore_errors=True)
    clean_stale_uploads()
    threading.Thread(target=_worker, daemon=True, name="genome-builds").start()


def seed(seed_dir: Path) -> dict | None:
    """First start with no genomes: queue the demo described by `<seed_dir>/seed.json`."""
    spec_p = seed_dir / "seed.json"
    if (
        not spec_p.is_file()
        or any(genomes_dir().glob("*/genome.sqlite"))
        or list_jobs(1)
    ):
        return None
    spec = json.loads(spec_p.read_text())
    return submit(
        {
            "kind": "new_genome",
            "name": spec["name"],
            "description": spec.get("description", ""),
            "user": "seed",
            "fasta": str(seed_dir / spec["fasta"]),
            "tracks": [
                t | {"path": str(seed_dir / t["path"])} for t in spec.get("tracks", [])
            ],
        }
    )
