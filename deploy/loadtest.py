#!/usr/bin/env python3
"""Load test for the lab server: N people at once, browsing or running analyses.

Standard library only, so it runs on the Mac and on the VM's own python3. Each virtual person signs
in once (one shared test account) and repeats its scenario with no pause between requests, so
"8 people" means 8 requests always in flight: a worst case, not a typical day. Divide the measured
throughput by what one real person asks for (the usage log on /admin) to get how many fit.

  deploy/loadtest.py --base https://localhost:8443 --user loadtest --password '…' --insecure \\
      --levels 1,2,4,8 --seconds 30 --container lab_server-app-1 --out docs/load_test.html

Scenarios (the one written in the report is what ran):
  browse  map list, one map, its sites and ORFs, the genome browser page
  heavy   BLAST of a 1.5 kb piece of a random map, primer matching against the oligo collection
"""

from __future__ import annotations

import argparse
import html
import json
import random
import ssl
import statistics
import subprocess
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime


class Session:
    def __init__(self, base: str, insecure: bool):
        self.base = base.rstrip("/")
        self.ctx = ssl._create_unverified_context() if insecure else None
        self.cookie = ""

    def req(
        self, method: str, path: str, body: dict | None = None, form: dict | None = None
    ):
        data, headers = None, {"Cookie": self.cookie, "Origin": self.base}
        if body is not None:
            data, headers["Content-Type"] = (
                json.dumps(body).encode(),
                "application/json",
            )
        elif form is not None:
            data = urllib.parse.urlencode(form).encode()
            headers["Content-Type"] = "application/x-www-form-urlencoded"
        r = urllib.request.Request(
            self.base + path, data=data, headers=headers, method=method
        )
        t0 = time.perf_counter()
        try:
            with urllib.request.urlopen(r, timeout=120, context=self.ctx) as resp:
                payload, status = resp.read(), resp.status
        except urllib.error.HTTPError as e:
            payload, status = e.read(), e.code
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            payload, status = str(e).encode(), 0
        return status, (time.perf_counter() - t0) * 1000, payload

    def login(self, user: str, password: str) -> None:
        class NoRedirect(urllib.request.HTTPRedirectHandler):
            def redirect_request(self, *a, **k):
                return None

        opener = urllib.request.build_opener(
            NoRedirect, urllib.request.HTTPSHandler(context=self.ctx)
        )
        r = urllib.request.Request(
            self.base + "/login",
            method="POST",
            headers={"Origin": self.base},
            data=urllib.parse.urlencode({"user": user, "password": password}).encode(),
        )
        try:
            opener.open(r, timeout=30)
        except urllib.error.HTTPError as e:
            self.cookie = (e.headers.get("set-cookie") or "").split(";")[0]
        if not self.cookie:
            raise SystemExit("sign-in failed: check --user/--password")


def browse(s: Session, keys: list[int], log):
    k = random.choice(keys)
    for label, path in (
        ("map list", "/api/maps"),
        ("one map", f"/api/maps/{k}"),
        ("sites", f"/api/maps/{k}/sites"),
        ("ORFs", f"/api/maps/{k}/orfs"),
        ("genome page", "/genomes/"),
    ):
        log(label, *s.req("GET", path)[:2])


def heavy(s: Session, keys: list[int], log):
    k = random.choice(keys)
    st, ms, body = s.req("GET", f"/api/maps/{k}")
    seq = json.loads(body).get("seq", "") if st == 200 else ""
    if len(seq) > 1600:
        i = random.randrange(0, len(seq) - 1500)
        log(
            "BLAST 1.5 kb",
            *s.req("POST", "/api/blast", body={"seq": seq[i : i + 1500]})[:2],
        )
    log("primer match", *s.req("GET", f"/api/maps/{k}/primers")[:2])


def docker_sampler(container: str, stop: threading.Event, out: list):
    while not stop.is_set():
        try:
            line = subprocess.run(
                [
                    "docker",
                    "stats",
                    "--no-stream",
                    "--format",
                    "{{.CPUPerc}}|{{.MemUsage}}",
                    container,
                ],
                capture_output=True,
                text=True,
                timeout=15,
            ).stdout.strip()
            cpu, mem = line.split("|")
            used = mem.split("/")[0].strip()
            mib = (
                float(used[:-3]) * (1024 if used.endswith("GiB") else 1)
                if used[-3:] in ("MiB", "GiB")
                else 0.0
            )
            out.append((float(cpu.rstrip("%")), mib))
        except (ValueError, subprocess.SubprocessError, OSError):
            pass
        stop.wait(1)


def run_level(args, scenario, people: int, keys: list[int], cookie: str) -> dict:
    rows: list[tuple[str, int, float]] = []
    lock = threading.Lock()

    def log(label, status, ms):
        with lock:
            rows.append((label, status, ms))

    stop, deadline = threading.Event(), time.time() + args.seconds
    samples: list = []
    sampler = threading.Thread(
        target=docker_sampler, args=(args.container, stop, samples), daemon=True
    )
    if args.container:
        sampler.start()

    def person():
        s = Session(args.base, args.insecure)
        s.cookie = cookie
        while time.time() < deadline:
            scenario(s, keys, log)

    threads = [threading.Thread(target=person) for _ in range(people)]
    t0 = time.time()
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    wall = time.time() - t0
    stop.set()
    by_label: dict[str, list] = {}
    for label, status, ms in rows:
        by_label.setdefault(label, []).append((status, ms))

    def stats(items):
        ok = sorted(ms for st, ms in items if st == 200)
        q = lambda p: round(ok[min(len(ok) - 1, int(p * len(ok)))]) if ok else None  # noqa: E731
        return {
            "n": len(items),
            "ok": len(ok),
            "busy": sum(st == 429 for st, _ in items),
            "failed": sum(st not in (200, 429) for st, _ in items),
            "p50": q(0.5),
            "p95": q(0.95),
            "max": round(ok[-1]) if ok else None,
            "per_s": round(len(ok) / wall, 2),
        }

    return {
        "people": people,
        "wall_s": round(wall, 1),
        "all": stats([(st, ms) for _, st, ms in rows]),
        "by_label": {k: stats(v) for k, v in sorted(by_label.items())},
        "cpu_peak": max((c for c, _ in samples), default=None),
        "cpu_mean": round(statistics.mean(c for c, _ in samples), 1)
        if samples
        else None,
        "mem_peak_mib": round(max((m for _, m in samples), default=0)) or None,
    }


def render(meta: dict, results: dict) -> str:
    e = html.escape
    sec = []
    for scen, levels in results.items():
        rows = []
        for lv in levels:
            a = lv["all"]
            rows.append(
                f"<tr class='lv'><td>{lv['people']}</td><td>all</td><td class=n>{a['ok']}</td><td class=n>{a['per_s']}</td>"
                f"<td class=n>{a['p50']}</td><td class=n>{a['p95']}</td><td class=n>{a['max']}</td><td class=n>{a['busy']}</td>"
                f"<td class=n>{a['failed']}</td><td class=n>{lv['cpu_peak'] if lv['cpu_peak'] is not None else '–'}</td>"
                f"<td class=n>{lv['mem_peak_mib'] or '–'}</td></tr>"
            )
            for label, b in lv["by_label"].items():
                rows.append(
                    f"<tr><td></td><td>{e(label)}</td><td class=n>{b['ok']}</td><td class=n>{b['per_s']}</td><td class=n>{b['p50']}</td>"
                    f"<td class=n>{b['p95']}</td><td class=n>{b['max']}</td><td class=n>{b['busy']}</td><td class=n>{b['failed']}</td><td></td><td></td></tr>"
                )
        sec.append(
            f"<h2>{e(scen)}</h2><table><thead><tr><th>People</th><th>Request</th><th>OK</th><th>OK/s</th><th>p50 ms</th>"
            f"<th>p95 ms</th><th>max ms</th><th>Busy 429</th><th>Failed</th><th>CPU peak %</th><th>RAM peak MiB</th></tr></thead>"
            f"<tbody>{''.join(rows)}</tbody></table>"
        )
    favicon = (
        "data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Crect x='4' y='14' width='6' "
        "height='14' fill='%232b7a78'/%3E%3Crect x='13' y='8' width='6' height='20' fill='%232b7a78'/%3E%3Crect x='22' y='3' "
        "width='6' height='25' fill='%23e07a5f'/%3E%3C/svg%3E"
    )
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light dark"><title>Load test {e(meta["date"])}</title>
<link rel="icon" href="{favicon}">
<style>
:root {{ --bg:#f7f7f5; --panel:#fff; --ink:#1d2327; --muted:#5f6b73; --line:#dde1e4; --accent:#2b7a78; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#15191c; --panel:#1d2327; --ink:#e6e9eb; --muted:#9aa5ad; --line:#333b41; --accent:#5fb3b0; }} }}
body {{ margin:0; background:var(--bg); color:var(--ink); font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif; }}
main {{ max-width:1000px; margin:0 auto; padding:20px 16px 40px; }}
h1 {{ font-size:20px; }} h2 {{ font-size:15px; color:var(--muted); text-transform:uppercase; letter-spacing:.04em; margin-top:28px; }}
.card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:12px 16px; }}
.wrap {{ overflow-x:auto; }} table {{ border-collapse:collapse; width:100%; font-size:13px; background:var(--panel); }}
th,td {{ padding:4px 8px; border-bottom:1px solid var(--line); text-align:left; }} th {{ color:var(--muted); font-weight:500; }}
td.n {{ text-align:right; font-variant-numeric:tabular-nums; }} tr.lv td {{ font-weight:600; border-top:2px solid var(--line); }}
code {{ font-family:ui-monospace,Menlo,monospace; font-size:12px; }} .muted {{ color:var(--muted); }}
@media print {{ body {{ background:#fff; color:#000; }} .card,table {{ border-color:#999; }} }}
</style></head><body><main>
<h1>Load test, {e(meta["date"])}</h1>
<div class="card"><p><b>Target</b> {e(meta["base"])} · <b>each level</b> {meta["seconds"]} s · <b>maps in the bank</b> {meta["maps"]}
 · <b>container</b> {e(meta["container"] or "not sampled")}</p>
<p class="muted">Every virtual person repeats its scenario with no pause, so "N people" means N requests always in flight, a worst
case. Times are measured at the client, through Caddy. "Busy 429" means the server turned an analysis away because both heavy-job
slots were taken. CPU % is per core (200 % = both cores of a 2-CPU host). The command that produced this page:</p>
<p><code>{e(meta["command"])}</code></p></div>
<div class="wrap">{"".join(sec)}</div>
</main></body></html>"""


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--base", required=True)
    ap.add_argument("--user", required=True)
    ap.add_argument("--password", required=True)
    ap.add_argument(
        "--insecure",
        action="store_true",
        help="accept a self-signed certificate (tls internal)",
    )
    ap.add_argument("--levels", default="1,2,4,8")
    ap.add_argument("--seconds", type=int, default=30)
    ap.add_argument("--scenarios", default="browse,heavy")
    ap.add_argument(
        "--container", default="", help="docker container to sample CPU and RAM from"
    )
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    s = Session(args.base, args.insecure)
    s.login(args.user, args.password)
    st, _, body = s.req("GET", "/api/maps")
    if st != 200:
        raise SystemExit(f"/api/maps answered {st}")
    keys = [
        m["key"]
        for m in json.loads(body)
        if m.get("length") and m["format"] in ("genbank", "fasta")
    ]
    print(f"{len(keys)} maps with a sequence")
    results: dict[str, list] = {}
    for name in args.scenarios.split(","):
        for n in (int(x) for x in args.levels.split(",")):
            lv = run_level(
                args, {"browse": browse, "heavy": heavy}[name], n, keys, s.cookie
            )
            a = lv["all"]
            print(
                f"{name:7} {n:2} people  ok={a['ok']:5} {a['per_s']:6}/s  p50={a['p50']} p95={a['p95']}  "
                f"busy={a['busy']} failed={a['failed']}  cpu_peak={lv['cpu_peak']} mem_peak={lv['mem_peak_mib']}"
            )
            results.setdefault(name, []).append(lv)
    if args.out:
        cmd = " ".join(
            a if a != args.password else "'…'" for a in __import__("sys").argv
        )
        meta = {
            "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
            "base": args.base,
            "seconds": args.seconds,
            "maps": len(keys),
            "container": args.container,
            "command": cmd,
        }
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(render(meta, results))
        with open(args.out.rsplit(".", 1)[0] + ".json", "w", encoding="utf-8") as f:
            json.dump({"meta": meta, "results": results}, f, indent=1)
        print("report:", args.out)


if __name__ == "__main__":
    main()
