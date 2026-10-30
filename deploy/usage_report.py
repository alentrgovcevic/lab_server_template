#!/usr/bin/env python3
"""Usage report for the admin, built on their own computer from the server's usage log.

The server only writes ``/state/usage/{usage,metrics}-YYYY-MM.jsonl`` (``plasmidbank/app/usage.py``).
This fetches them and answers: who uses what and for how long, where people wait (the bottleneck),
and whether the VM needs more cores, RAM or disk. Standard library only.

  deploy/usage_report.py --host lab-server    # fetch over ssh (an ~/.ssh/config host), last 30 days, open the page
  deploy/usage_report.py --days 90
  deploy/usage_report.py --container lab_server-app-1   # a local stack (quickstart) instead
  deploy/usage_report.py --dir ~/PlasmidBank/usage_mirror --no-fetch

The page is personal data (who did what), so it goes to ~/PlasmidBank/reports/, never into the repo.
The fetched copy lives in ~/PlasmidBank/usage_mirror/ and is replaced on every fetch.
"""

from __future__ import annotations

import argparse
import html
import json
import statistics
import subprocess
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

HOME = Path.home() / "PlasmidBank"
SESSION_GAP_S = 30 * 60  # a pause longer than this starts a new session
WINDOW_S = 300  # "at the same time" = active within the same 5 minutes
HEAVY = ("BLAST", "Alignment", "Clone verification", "Golden Gate")

# route template -> tool, first match wins (prefixes)
TOOLS = [
    ("genomes:", "Genome browser"),
    ("/api/blast", "BLAST"),
    ("/api/align", "Alignment"),
    ("/api/maps/{key}/verify", "Clone verification"),
    ("/api/trace", "Clone verification"),
    ("/api/assemble", "Golden Gate"),
    ("/api/build", "Map builder"),
    ("/api/maps/{key}/primers", "Primers"),
    ("/api/oligos", "Primers"),
    ("/api/maps/{key}/design", "Primer design"),
    ("/api/maps/{key}/sdm", "Primer design"),
    ("/api/design", "Primer design"),
    ("/api/maps/{key}/sites", "Enzymes and digests"),
    ("/api/maps/{key}/digest", "Enzymes and digests"),
    ("/api/enzymes", "Enzymes and digests"),
    ("/api/lab/enzymes", "Enzymes and digests"),
    ("/api/gel", "Enzymes and digests"),
    ("/api/maps/{key}/orfs", "ORFs and translation"),
    ("/api/maps/{key}/translate", "ORFs and translation"),
    ("/api/maps/{key}/clean", "Cleaner"),
    ("/api/clean", "Cleaner"),
    ("/api/maps/{key}/save", "Editing and annotation"),
    ("/api/maps/{key}/annotate", "Editing and annotation"),
    ("/api/annotations", "Editing and annotation"),
    ("/api/upload", "Uploads"),
    ("/api/projects", "Projects"),
    ("/api/feedback", "Feedback"),
    ("/api/admin", "Admin"),
    ("/admin", "Admin"),
    ("/api/maps", "Browsing the bank"),
    ("/api/stats", "Browsing the bank"),
    ("/bank", "Browsing the bank"),
    ("/api/me", "Dashboard"),
    ("/api/links", "Dashboard"),
    ("/login", "Sign-in"),
    ("/logout", "Sign-in"),
    ("/api/account", "Sign-in"),
]


def tool_of(route: str) -> str:
    if route == "/":
        return "Dashboard"
    return next((name for prefix, name in TOOLS if route.startswith(prefix)), "Other")


def fetch(dest: Path, container: str, host: str | None) -> None:
    """Replace ``dest`` with the container's /state/usage (ssh to ``host`` when given)."""
    cmd = ["docker", "cp", f"{container}:/state/usage", "-"]
    if host:
        cmd = ["ssh", host, " ".join(cmd)]
    tar = subprocess.run(cmd, capture_output=True, check=True, timeout=300).stdout
    tmp = dest.with_name(dest.name + ".new")
    subprocess.run(["rm", "-rf", str(tmp)], check=True)
    tmp.mkdir(parents=True)
    subprocess.run(
        ["tar", "-x", "-C", str(tmp), "--strip-components", "1"], input=tar, check=True
    )
    subprocess.run(["rm", "-rf", str(dest)], check=True)
    tmp.rename(dest)


def load(folder: Path, kind: str, since: float) -> list[dict]:
    first = datetime.fromtimestamp(since).strftime("%Y-%m")
    out = []
    for p in sorted(folder.glob(f"{kind}-*.jsonl")):
        if p.stem.split("-", 1)[1] < first:
            continue
        for line in p.read_text(encoding="utf-8").splitlines():
            try:
                r = json.loads(line)
            except ValueError:
                continue  # a line cut short by a crash
            if r.get("t", 0) >= since:
                out.append(r)
    return sorted(out, key=lambda r: r["t"])


def pct(xs, q):
    xs = sorted(x for x in xs if x is not None)
    return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else None


def sessions(times: list[float]) -> list[tuple[float, float]]:
    """Split one person's request times at pauses longer than SESSION_GAP_S; (start, end) each."""
    out: list[list[float]] = []
    for t in sorted(times):
        if out and t - out[-1][1] <= SESSION_GAP_S:
            out[-1][1] = t
        else:
            out.append([t, t])
    return [(a, b) for a, b in out]


def session_minutes(s: tuple[float, float]) -> float:
    return (s[1] - s[0]) / 60 + 1  # a single click counts as one minute


def analyse(rows: list[dict], metrics: list[dict]) -> dict:
    people = [r for r in rows if r.get("user")]
    by_user: dict[str, list] = defaultdict(list)
    for r in people:
        by_user[r["user"]].append(r)

    users = []
    for u, rs in by_user.items():
        ss = sessions([r["t"] for r in rs])
        tool_minutes: dict[str, set] = defaultdict(set)
        for r in rs:
            tool_minutes[tool_of(r["route"])].add(int(r["t"] // 60))
        users.append(
            {
                "user": u,
                "requests": len(rs),
                "days": len({datetime.fromtimestamp(r["t"]).date() for r in rs}),
                "sessions": len(ss),
                "hours": round(sum(session_minutes(s) for s in ss) / 60, 1),
                "median_session_min": round(
                    statistics.median(session_minutes(s) for s in ss)
                ),
                "top": [
                    k
                    for k, v in sorted(tool_minutes.items(), key=lambda kv: -len(kv[1]))
                    if k not in ("Sign-in", "Dashboard")
                ][:4],
                "last": datetime.fromtimestamp(rs[-1]["t"]).strftime("%Y-%m-%d %H:%M"),
            }
        )
    users.sort(key=lambda x: -x["hours"])

    per_tool: dict[str, dict] = defaultdict(
        lambda: {"users": set(), "minutes": set(), "ms": [], "busy": 0, "errors": 0}
    )
    for r in people:
        d = per_tool[tool_of(r["route"])]
        d["users"].add(r["user"])
        d["minutes"].add((r["user"], int(r["t"] // 60)))
        d["ms"].append(r.get("ms", 0))
        d["busy"] += r.get("status") == 429
        d["errors"] += r.get("status", 0) >= 500
    tools = sorted(
        (
            {
                "tool": k,
                "users": len(v["users"]),
                "person_hours": round(len(v["minutes"]) / 60, 1),
                "requests": len(v["ms"]),
                "p50": round(pct(v["ms"], 0.5)),
                "p95": round(pct(v["ms"], 0.95)),
                "wait_s": round(sum(v["ms"]) / 1000, 1),
                "busy": v["busy"],
                "errors": v["errors"],
            }
            for k, v in per_tool.items()
        ),
        key=lambda x: -x["person_hours"],
    )

    windows: dict[int, set] = defaultdict(set)
    hours: dict[int, set] = defaultdict(set)  # hour of day -> (user, day)
    for r in people:
        windows[int(r["t"] // WINDOW_S)].add(r["user"])
        dt = datetime.fromtimestamp(r["t"])
        hours[dt.hour].add((r["user"], dt.date()))
    peak_w = max(windows.items(), key=lambda kv: len(kv[1]), default=(0, set()))
    peak_in = max(rows, key=lambda r: r.get("inflight", 0), default={})
    heavy = [r for r in people if tool_of(r["route"]) in HEAVY]
    return {
        "requests": len(people),
        "users": users,
        "tools": tools,
        "peak_users": len(peak_w[1]),
        "peak_users_at": datetime.fromtimestamp(peak_w[0] * WINDOW_S).strftime(
            "%Y-%m-%d %H:%M"
        )
        if peak_w[1]
        else None,
        "peak_inflight": peak_in.get("inflight", 0),
        "heavy": len(heavy),
        "busy": sum(r.get("status") == 429 for r in heavy),
        "by_hour": {h: len(v) for h, v in sorted(hours.items())},
        "resources": resources(metrics),
    }


def resources(m: list[dict]) -> dict | None:
    if not m:
        return None
    used = [
        (r["mem_total"] - r["mem_avail"]) / r["mem_total"] * 100
        for r in m
        if r.get("mem_total")
    ]
    first, last = m[0], m[-1]
    days = max((last["t"] - first["t"]) / 86400, 1e-9)
    growth = (last["disk_used"] - first["disk_used"]) / days if days >= 1 else None

    def days_to(frac):
        left = frac * last["disk_total"] - last["disk_used"]
        return (
            None
            if not growth or growth <= 0
            else round(left / growth)
            if left > 0
            else 0
        )

    return {
        "minutes": len(m),
        "cpus": last.get("cpus"),
        "cpu_p50": pct([r["cpu"] for r in m], 0.5),
        "cpu_p95": pct([r["cpu"] for r in m], 0.95),
        "cpu_max": max(r["cpu"] or 0 for r in m),
        "cpu_busy_share": round(100 * sum((r["cpu"] or 0) > 80 for r in m) / len(m), 1),
        "load_max": max(r["load1"] for r in m),
        "mem_total": last["mem_total"],
        "mem_used_p95": round(pct(used, 0.95), 1),
        "mem_used_max": round(max(used), 1),
        "mem_avail_min": min(r["mem_avail"] for r in m),
        "app_mem_max": max((r.get("app_mem") or 0) for r in m),
        "disk_total": last["disk_total"],
        "disk_used": last["disk_used"],
        "disk_growth_gib_day": round(growth, 2) if growth is not None else None,
        "days_to_75": days_to(0.75),
        "days_to_90": days_to(0.90),
    }


def verdicts(a: dict) -> list[tuple[str, str, str]]:
    """(topic, verdict, why). Rules of thumb, stated as such on the page."""
    out = []
    r = a["resources"]
    busy_share = 100 * a["busy"] / a["heavy"] if a["heavy"] else 0
    if r is None:
        out.append(
            (
                "Cores",
                "no data yet",
                "no metrics rows (the sampler writes one a minute once deployed)",
            )
        )
    else:
        many = r["cpu_p95"] is not None and r["cpu_p95"] > 80
        out.append(
            (
                "Cores",
                "more would help" if many or busy_share > 2 else "enough",
                f"CPU above 80 % in {r['cpu_busy_share']} % of minutes (p95 {r['cpu_p95']} %); "
                f"{busy_share:.1f} % of analyses turned away as busy. Rule of thumb: act when p95 CPU > 80 % or > 2 % busy. "
                "Browsing runs in one process and uses one core, so extra cores speed up analyses, not page loads",
            )
        )
        low = r["mem_avail_min"] < 0.15 * r["mem_total"]
        out.append(
            (
                "RAM",
                "add some" if low else "enough",
                f"at most {r['mem_used_max']} % in use (p95 {r['mem_used_p95']} %), lowest free {r['mem_avail_min']} MiB "
                f"of {r['mem_total']} MiB; the bank app peaked at {r['app_mem_max']} MiB. No swap, so running out kills a "
                "process. Rule of thumb: act when free memory drops below 15 %",
            )
        )
        if r["disk_growth_gib_day"] is None:
            disk = f"{r['disk_used']} of {r['disk_total']} GiB used; growth needs at least a day of data"
        else:
            disk = (
                f"{r['disk_used']} of {r['disk_total']} GiB used, growing {r['disk_growth_gib_day']} GiB/day: "
                f"75 % in {r['days_to_75'] if r['days_to_75'] is not None else '—'} days, 90 % in "
                f"{r['days_to_90'] if r['days_to_90'] is not None else '—'} days"
            )
        near = r["days_to_75"] is not None and r["days_to_75"] < 180
        out.append(("Disk", "plan more" if near else "enough", disk))
    slow = max(a["tools"], key=lambda t: t["wait_s"], default=None)
    if slow:
        out.append(
            (
                "Bottleneck",
                slow["tool"],
                f"people waited {slow['wait_s']} s in total on it (p95 {slow['p95']} ms over {slow['requests']} requests)",
            )
        )
    return out


def render(a: dict, days: int, source: str) -> str:
    e = html.escape
    r = a["resources"]
    users = (
        "".join(
            f"<tr><td>{e(u['user'])}</td><td class=n>{u['hours']}</td><td class=n>{u['sessions']}</td><td class=n>{u['median_session_min']}</td>"
            f"<td class=n>{u['days']}</td><td>{e(', '.join(u['top']))}</td><td class=muted>{u['last']}</td></tr>"
            for u in a["users"]
        )
        or "<tr><td colspan=7 class=muted>Nobody yet.</td></tr>"
    )
    tools = "".join(
        f"<tr><td>{e(t['tool'])}</td><td class=n>{t['users']}</td><td class=n>{t['person_hours']}</td><td class=n>{t['requests']}</td>"
        f"<td class=n>{t['p50']}</td><td class=n>{t['p95']}</td><td class=n>{t['wait_s']}</td><td class=n>{t['busy']}</td>"
        f"<td class=n>{t['errors']}</td></tr>"
        for t in a["tools"]
    )
    top = max(a["by_hour"].values(), default=1)
    hours = "".join(
        f"<tr><td class=n>{h:02d}:00</td><td><span class=bar style='width:{100 * n / top:.0f}%'></span></td>"
        f"<td class=n>{n}</td></tr>"
        for h, n in a["by_hour"].items()
    )
    verd = "".join(
        f"<tr><th>{e(t)}</th><td><b>{e(v)}</b></td><td class=muted>{e(w)}</td></tr>"
        for t, v, w in verdicts(a)
    )
    res = (
        ""
        if r is None
        else (
            f"<p>{r['minutes']:,} one-minute readings · {r['cpus']} CPUs · CPU p50 {r['cpu_p50']} %, p95 {r['cpu_p95']} %, max {r['cpu_max']} % · "
            f"load max {r['load_max']} · memory max {r['mem_used_max']} % of {r['mem_total']} MiB · disk {r['disk_used']} of {r['disk_total']} GiB</p>"
        )
    )
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light dark"><title>Lab server usage {datetime.now():%Y-%m-%d}</title>
<link rel="icon" href="data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Ccircle cx='16' cy='16' r='11' fill='none' stroke='%232b7a78' stroke-width='4'/%3E%3C/svg%3E">
<style>
:root {{ --bg:#f7f7f5; --panel:#fff; --ink:#1d2327; --muted:#5f6b73; --line:#dde1e4; --accent:#2b7a78; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#15191c; --panel:#1d2327; --ink:#e6e9eb; --muted:#9aa5ad; --line:#333b41; --accent:#5fb3b0; }} }}
body {{ margin:0; background:var(--bg); color:var(--ink); font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif; }}
main {{ max-width:1000px; margin:0 auto; padding:20px 16px 40px; }}
h1 {{ font-size:20px; }} h2 {{ font-size:15px; color:var(--muted); text-transform:uppercase; letter-spacing:.04em; margin-top:28px; }}
.card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:12px 16px; overflow-x:auto; }}
table {{ border-collapse:collapse; width:100%; font-size:13px; }} th,td {{ padding:4px 8px; border-bottom:1px solid var(--line); text-align:left; vertical-align:top; }}
th {{ color:var(--muted); font-weight:500; }} td.n {{ text-align:right; font-variant-numeric:tabular-nums; white-space:nowrap; }}
.muted {{ color:var(--muted); }} .bar {{ display:inline-block; height:10px; background:var(--accent); border-radius:2px; }}
@media print {{ body {{ background:#fff; color:#000; }} }}
</style></head><body><main>
<h1>Lab server usage, last {days} days</h1>
<p class="muted">{a["requests"]:,} requests by {len(a["users"])} people · most people active in the same 5 minutes: {a["peak_users"]}
{f"({a['peak_users_at']})" if a["peak_users_at"] else ""} · most requests in progress at once: {a["peak_inflight"]} · source: {e(source)} ·
built {datetime.now():%Y-%m-%d %H:%M}</p>
<h2>Do we need more?</h2><div class="card"><table>{verd}</table>
<p class="muted">The thresholds are rules of thumb, not measurements. {res}</p></div>
<h2>People</h2><div class="card"><table><thead><tr><th>Who</th><th class=n>Hours</th><th class=n>Sessions</th><th class=n>Median session, min</th>
<th class=n>Days</th><th>Uses most</th><th>Last seen</th></tr></thead><tbody>{users}</tbody></table>
<p class="muted">A session ends after {SESSION_GAP_S // 60} minutes without a request; its length runs from the first to the last request, plus one minute.</p></div>
<h2>Tools</h2><div class="card"><table><thead><tr><th>Tool</th><th class=n>People</th><th class=n>Person-hours</th><th class=n>Requests</th>
<th class=n>p50 ms</th><th class=n>p95 ms</th><th class=n>Waiting, s</th><th class=n>Busy</th><th class=n>Errors</th></tr></thead><tbody>{tools}</tbody></table>
<p class="muted">Person-hours count the minutes in which someone used the tool. Waiting adds up every response time: the largest number in that column is where people wait most.
Busy = an analysis turned away because both heavy-job slots were taken. Genome-browser times are the sign-in check only.</p></div>
<h2>Time of day</h2><div class="card"><table><thead><tr><th>Hour</th><th></th><th class=n>Person-days</th></tr></thead><tbody>{hours}</tbody></table></div>
</main></body></html>"""


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--days", type=int, default=30)
    ap.add_argument("--host", default="lab-server", help="ssh host running the stack (an ~/.ssh/config name)")
    ap.add_argument(
        "--container",
        default="",
        help="a LOCAL container instead of the ssh host (e.g. lab_server-app-1)",
    )
    ap.add_argument("--dir", type=Path, default=HOME / "usage_mirror")
    ap.add_argument("--no-fetch", action="store_true")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--no-open", action="store_true")
    args = ap.parse_args()
    if not args.no_fetch:
        try:
            if args.container:
                fetch(args.dir, args.container, None)
            else:
                fetch(args.dir, "lab_server-app-1", args.host)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            err = (getattr(e, "stderr", b"") or b"").decode(errors="replace").strip()
            raise SystemExit(f"could not fetch the usage log ({err or e}). Does `ssh {args.host}` work? "
                             f"--no-fetch builds the page from the last copy in {args.dir}")
    since = time.time() - args.days * 86400
    a = analyse(load(args.dir, "usage", since), load(args.dir, "metrics", since))
    out = args.out or HOME / "reports" / f"usage_{datetime.now():%Y-%m-%d}.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render(a, args.days, args.container or args.host), encoding="utf-8")
    print(f"{a['requests']:,} requests, {len(a['users'])} people -> {out}")
    for t, v, w in verdicts(a):
        print(f"  {t}: {v} — {w}")
    if not args.no_open:
        subprocess.run(["open", str(out)], check=False)


if __name__ == "__main__":
    main()
