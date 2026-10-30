"""BLAST a sequence against the map bank: "which plasmids contain this?"

Needs NCBI BLAST+ on PATH (Mac: `brew install blast`; the Docker image installs `ncbi-blast+`).

The database is one FASTA record per map with a sequence, built with makeblastdb into
`<state>/blast/`. A circular map is stored with its first WRAP bases appended, so a hit that spans
the origin is still found; positions are folded back onto the map afterwards. Records are named
`m<n>` and a JSON sidecar maps them to the map's file and ID, because bank keys change on restart.
The database is rebuilt whenever the set of maps changes (a fingerprint of paths and lengths).
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import threading
from pathlib import Path

WRAP = 2000
FIELDS = "sseqid pident length qstart qend sstart send evalue bitscore qlen"
# one lock per database folder: a build for one folder must not hold up a search in another (each app
# start warms its own folder, and tests start many apps; a shared lock made searches wait past 30 s)
_locks: dict[Path, threading.Lock] = {}
_locks_guard = threading.Lock()


def _lock_for(folder: Path) -> threading.Lock:
    with _locks_guard:
        return _locks.setdefault(folder.resolve(), threading.Lock())


def available() -> str | None:
    """The blastn version line, or None when BLAST+ is not installed."""
    exe = shutil.which("blastn")
    if not exe or not shutil.which("makeblastdb"):
        return None
    try:
        return subprocess.run(
            [exe, "-version"], capture_output=True, text=True, timeout=20
        ).stdout.splitlines()[0]
    except (OSError, subprocess.SubprocessError, IndexError):
        return None


def _usable(maps) -> list:
    return [m for m in maps if m.seq and len(m.seq) >= 20]


def fingerprint(maps) -> str:
    h = hashlib.sha256()
    for m in _usable(maps):
        h.update(f"{m.path}\t{len(m.seq)}\t{m.circular}\n".encode())
    return h.hexdigest()[:16]


def ensure_db(maps, folder: Path) -> dict:
    """Build the database if it is missing or the maps changed. Returns the sidecar."""
    folder.mkdir(parents=True, exist_ok=True)
    side = folder / "bank.json"
    fp = fingerprint(maps)
    with _lock_for(folder):
        try:
            meta = json.loads(side.read_text())
            # a big database is split into bank.00.nsq, bank.01.nsq, …
            if meta.get("fingerprint") == fp and any(folder.glob("bank*.nsq")):
                return meta
        except (OSError, ValueError):
            pass
        records = {}
        fasta = folder / "bank.fa"
        with fasta.open("w") as fh:
            for i, m in enumerate(_usable(maps)):
                rid = f"m{i}"
                seq = m.seq.upper()
                stored = seq + seq[: min(WRAP, len(seq))] if m.circular else seq
                fh.write(f">{rid}\n")
                for k in range(0, len(stored), 80):
                    fh.write(stored[k : k + 80] + "\n")
                records[rid] = {
                    "path": m.path,
                    "plasmid_id": m.plasmid_id,
                    "name": m.name,
                    "length": len(seq),
                    "circular": m.circular,
                }
        r = subprocess.run(
            [
                "makeblastdb",
                "-in",
                str(fasta),
                "-dbtype",
                "nucl",
                "-out",
                str(folder / "bank"),
                "-parse_seqids",
            ],
            capture_output=True,
            text=True,
            timeout=600,
        )
        if r.returncode:
            raise RuntimeError(f"makeblastdb failed: {r.stderr.strip()[:300]}")
        fasta.unlink(missing_ok=True)
        meta = {"fingerprint": fp, "records": records, "maps": len(records)}
        side.write_text(json.dumps(meta))
        return meta


def search(
    query: str, maps, folder: Path, max_hits: int = 25, evalue: float | None = None
) -> dict:
    q = "".join(query.upper().split())
    if len(q) < 12:
        raise ValueError("give at least 12 bases to BLAST")
    if set(q) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("the query may only contain DNA letters")
    meta = ensure_db(maps, folder)
    short = len(q) < 50
    cmd = [
        "blastn",
        "-db",
        str(folder / "bank"),
        "-outfmt",
        f"6 {FIELDS}",
        "-task",
        "blastn-short" if short else ("megablast" if len(q) > 2000 else "blastn"),
        "-evalue",
        str(evalue if evalue is not None else (1 if short else 1e-10)),
        "-max_target_seqs",
        str(max(max_hits * 3, 50)),
        "-max_hsps",
        "5",
        "-num_threads",
        "2",
    ]
    r = subprocess.run(
        cmd, input=f">query\n{q}\n", capture_output=True, text=True, timeout=180
    )
    if r.returncode:
        raise RuntimeError(f"blastn failed: {r.stderr.strip()[:300]}")
    best: dict[str, dict] = {}
    for line in r.stdout.splitlines():
        sid, pid, ln, qs, qe, ss, se, ev, bits, qlen = line.split("\t")
        rec = meta["records"].get(sid.split("|")[-1])
        if not rec:
            continue
        L = rec["length"]
        ss, se = int(ss), int(se)
        strand = "+" if se >= ss else "-"
        lo, hi = min(ss, se), max(ss, se)
        if (
            rec["circular"] and lo > L
        ):  # wholly inside the appended copy: the same hit seen twice
            continue
        hit = {
            "plasmid_id": rec["plasmid_id"],
            "name": rec["name"],
            "path": rec["path"],
            "length": L,
            "identity": round(float(pid), 1),
            "aligned": int(ln),
            "coverage": round(100 * int(ln) / int(qlen), 1),
            "query_from": int(qs),
            "query_to": int(qe),
            "strand": strand,
            "map_from": (lo - 1) % L + 1,
            "map_to": (hi - 1) % L + 1,
            "evalue": float(ev),
            "bits": float(bits),
        }
        key = rec["path"]
        if key not in best or hit["bits"] > best[key]["bits"]:
            best[key] = {**hit, "hsps": best.get(key, {}).get("hsps", 0) + 1}
        else:
            best[key]["hsps"] += 1
    hits = sorted(best.values(), key=lambda h: (-h["bits"], h["plasmid_id"] or ""))[
        :max_hits
    ]
    return {
        "query_length": len(q),
        "hits": hits,
        "database_maps": meta["maps"],
        "program": cmd[cmd.index("-task") + 1],
    }
