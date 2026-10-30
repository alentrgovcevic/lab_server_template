"""Usage log. The server only writes; an admin analyses it on their own computer with
``deploy/usage_report.py`` (who uses what and for how long, the bottleneck, cores/RAM/disk).

Two files per month in ``<state>/usage/``, one JSON object per line:

- ``usage-YYYY-MM.jsonl``, one row per request: ``{t, user, m, route, status, ms, inflight}``.
  ``route`` is the route template (``/api/maps/{key}/primers``), never the raw path. Genome-browser
  requests reach this app only as Caddy's ``/auth/check``; they are logged as ``genomes:<path>``
  (Caddy has already stripped ``/genomes``), numbers folded to ``{n}``, and their ``ms`` is the
  sign-in check alone. Their scripts, styles and icons are not logged. ``inflight`` counts the
  requests this process was handling when the row's request started, itself included.
- ``metrics-YYYY-MM.jsonl``, one row a minute: ``{t, cpu, load1, cpus, mem_total, mem_avail,
  app_mem, disk_total, disk_used, inflight}``. Inside the container ``/proc/stat``, ``/proc/meminfo``
  and the load average are the HOST's (checked on an Ubuntu VM: MemTotal matched ``free``);
  ``app_mem`` is this container's cgroup. Memory in MiB, disk in GiB, cpu in % of all cores.
  Nothing is written where ``/proc`` does not exist (a Mac running the app directly).
"""

from __future__ import annotations

import json
import os
import re
import threading
import time
from pathlib import Path

_LOCK = threading.Lock()
_NUM = re.compile(r"/\d+(?=/|$)")
SKIP = ("/healthz", "/favicon.ico")
STATIC_SUFFIXES = (".js", ".css", ".png", ".svg", ".ico", ".woff", ".woff2", ".map")
EVERY_S = 60


def folded(path: str) -> str:
    return _NUM.sub("/{n}", path.split("?", 1)[0])[:200]


def _append(folder: Path, kind: str, row: dict) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    line = json.dumps(row, separators=(",", ":")) + "\n"
    with (
        _LOCK,
        open(
            folder / f"{kind}-{time.strftime('%Y-%m')}.jsonl", "a", encoding="utf-8"
        ) as f,
    ):
        f.write(line)


def record(folder: Path, row: dict) -> None:
    _append(folder, "usage", row)


def sample(
    state: Path,
    prev: tuple | None,
    proc: Path = Path("/proc"),
    cgroup: Path = Path("/sys/fs/cgroup"),
) -> tuple[dict, tuple]:
    """One resource reading. CPU % comes from the change in ``/proc/stat`` since ``prev``."""
    v = [int(x) for x in (proc / "stat").read_text().splitlines()[0].split()[1:]]
    idle, total = v[3] + v[4], sum(v)  # idle + iowait
    cpu = None
    if prev and total > prev[1]:
        cpu = round(100 * (1 - (idle - prev[0]) / (total - prev[1])), 1)
    mem = {}
    for line in (proc / "meminfo").read_text().splitlines():
        k, _, val = line.partition(":")
        mem[k] = int(val.split()[0]) if val.split() else 0
    try:
        app_mem = round(int((cgroup / "memory.current").read_text()) / 1048576)
    except (OSError, ValueError):
        app_mem = None
    st = os.statvfs(state)
    gib = st.f_frsize / 1024**3
    row = {
        "t": round(time.time()),
        "cpu": cpu,
        "load1": float((proc / "loadavg").read_text().split()[0]),
        "cpus": os.cpu_count(),
        "mem_total": mem.get("MemTotal", 0) // 1024,
        "mem_avail": mem.get("MemAvailable", 0) // 1024,
        "app_mem": app_mem,
        "disk_total": round(st.f_blocks * gib, 1),
        "disk_used": round((st.f_blocks - st.f_bfree) * gib, 1),
    }
    return row, (idle, total)


_target: dict = {}
_started = False


def start_sampler(folder: Path, state: Path, inflight) -> None:
    """Write one metrics row a minute. The thread starts once per process; a later call (the tests
    reload the server) only changes where it writes and whom it asks for ``inflight``."""
    global _started
    _target.update(folder=folder, state=state, inflight=inflight)
    if _started or not Path("/proc/stat").exists():
        return
    _started = True

    def loop():
        prev = None
        while True:
            try:
                row, prev = sample(_target["state"], prev)
                if row["cpu"] is not None:
                    _append(
                        _target["folder"],
                        "metrics",
                        {**row, "inflight": _target["inflight"]()},
                    )
            except (OSError, ValueError, IndexError):
                pass  # logging must never take the site down
            time.sleep(EVERY_S)

    threading.Thread(target=loop, name="usage-metrics", daemon=True).start()
