"""Batch GenBank cleaner: write cleaned copies plus an HTML report. Originals are never touched.

    ../.venv/bin/python -m plasmidbank.clean ../examples/bank --out cleaned
    ../.venv/bin/python -m plasmidbank.clean "/Volumes/LabShare/Plasmid map bank/PMC7001-8000" --out ~/cleaned \\
        --rules +fragments_to_comment,+drop_colours --oligos ../examples/oligos.tab
    ../.venv/bin/python -m plasmidbank.clean ../examples/bank --dry-run      # report only, no GenBank written

--rules starts from the default set; "+rule" adds an opt-in rule, "-rule" removes a default one.
"""

from __future__ import annotations

import argparse
import html
import sys
from datetime import datetime
from pathlib import Path

from .bank import SKIP_SUFFIXES
from .formats import sniff
from .oligos import load_oligos
from .tools.clean import (
    ALL_RULES,
    DEFAULT_ON,
    OPT_IN,
    CleanResult,
    clean_record,
    read_genbank,
)


def parse_rules(spec: str) -> list[str]:
    rules = list(DEFAULT_ON)
    for tok in filter(None, (t.strip() for t in spec.split(","))):
        name = tok.lstrip("+-")
        if name not in ALL_RULES:
            raise SystemExit(f"unknown rule {name!r}; rules: {', '.join(ALL_RULES)}")
        if tok.startswith("-"):
            rules = [r for r in rules if r != name]
        elif name not in rules:
            rules.append(name)
    return rules


def genbank_files(paths: list[Path], out: Path | None):
    for root in paths:
        files = (
            [root]
            if root.is_file()
            else sorted(p for p in root.rglob("*") if p.is_file())
        )
        for p in files:
            if out and out.resolve() in p.resolve().parents:
                continue
            if (
                p.name.startswith(".")
                or p.suffix.lower() in SKIP_SUFFIXES
                or p.stat().st_size > 50 * 1024 * 1024
            ):
                continue
            with open(p, "rb") as fh:
                if sniff(fh.read(4096)) == "genbank":
                    yield root, p


def report_html(
    rows: list[tuple[Path, CleanResult | None, str]], rules: list[str]
) -> str:
    def table(items):
        return "".join(
            f"<tr><td>{html.escape(c.rule)}</td><td>{html.escape(c.action)}</td>"
            f"<td>{html.escape(c.feature)}</td><td>{html.escape(c.detail)}</td></tr>"
            for c in items
        )

    body, n_changed, n_flagged = [], 0, 0
    for path, res, err in rows:
        if err:
            body.append(
                f"<details open><summary class='bad'>{html.escape(str(path))} — {html.escape(err)}</summary></details>"
            )
            continue
        n_changed += bool(res.changes)
        n_flagged += bool(res.flags)
        if not res.changes and not res.flags:
            continue
        body.append(
            f"<details><summary>{html.escape(str(path))} — {len(res.changes)} change(s), {len(res.flags)} flag(s), "
            f"{res.n_before}→{res.n_after} features</summary>"
            + (
                f"<h4>Changes</h4><table><tr><th>rule</th><th>action</th><th>feature</th><th>detail</th></tr>{table(res.changes)}</table>"
                if res.changes
                else ""
            )
            + (
                f"<h4>Flags (not changed)</h4><table><tr><th>rule</th><th>action</th><th>feature</th><th>detail</th></tr>{table(res.flags)}</table>"
                if res.flags
                else ""
            )
            + "</details>"
        )
    when = datetime.now().strftime("%Y-%m-%d %H:%M")
    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light dark"><title>GenBank clean report</title>
<link rel="icon" href="data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Ccircle cx='16' cy='16' r='11' fill='none' stroke='%232b7a78' stroke-width='4'/%3E%3C/svg%3E">
<style>:root{{--bg:#fff;--ink:#1d2327;--muted:#5f6b73;--line:#dde1e4;--bad:#c0392b}}
@media (prefers-color-scheme: dark){{:root{{--bg:#15191c;--ink:#e6e9eb;--muted:#9aa5ad;--line:#333b41;--bad:#f07b6e}}}}
body{{font:14px/1.45 system-ui,sans-serif;background:var(--bg);color:var(--ink);margin:24px;max-width:1200px}}
summary{{cursor:pointer;padding:4px 0}} .bad{{color:var(--bad)}} table{{border-collapse:collapse;width:100%;font-size:12px;margin-bottom:8px}}
td,th{{border-bottom:1px solid var(--line);padding:3px 6px;text-align:left;vertical-align:top}} th{{color:var(--muted)}}
@media print{{details{{display:block}} details>*{{display:block}}}}</style></head><body>
<h1>GenBank clean report</h1><p>{when} · {len(rows)} file(s) · {n_changed} changed · {n_flagged} with flags · rules: {html.escape(", ".join(rules))}</p>
<p style="color:var(--muted)">Originals are untouched; the sequence and topology of every cleaned copy were verified identical to the input.</p>
{"".join(body)}</body></html>"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("paths", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, default=Path("cleaned"))
    ap.add_argument(
        "--rules",
        default="",
        help=f"defaults: {','.join(DEFAULT_ON)}; opt-in: {','.join(OPT_IN)}",
    )
    ap.add_argument(
        "--oligos",
        type=Path,
        help="oligo export, so clc_primer_hits keeps primers named in the collection",
    )
    ap.add_argument("--dry-run", action="store_true", help="write only the report")
    a = ap.parse_args(argv)
    rules = parse_rules(a.rules)
    names = None
    if a.oligos:
        ol = load_oligos(a.oligos)
        names = {o.name for o in ol} | {o.oligo_id for o in ol}
    out = a.out.expanduser()
    for p in a.paths:
        if p.resolve() == out.resolve():
            raise SystemExit(
                "--out must not be an input folder: originals are never overwritten"
            )
    rows = []
    for root, path in genbank_files(a.paths, out):
        rel = path.relative_to(root) if root.is_dir() else Path(path.name)
        try:
            res = clean_record(read_genbank(path), rules, names)
        except (
            Exception
        ) as e:  # report and continue; one bad file must not stop a 6,000-file run
            rows.append((rel, None, f"{type(e).__name__}: {e}"))
            continue
        rows.append((rel, res, ""))
        if not a.dry_run:
            dest = out / rel.with_suffix(".gb")
            if dest.resolve() == path.resolve():
                raise SystemExit(f"refusing to overwrite the original {path}")
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(res.genbank(), encoding="utf-8")
    out.mkdir(parents=True, exist_ok=True)
    (out / "clean_report.html").write_text(report_html(rows, rules), encoding="utf-8")
    changed = sum(1 for _, r, e in rows if r and r.changes)
    print(
        f"{len(rows)} GenBank files, {changed} changed; report: {out / 'clean_report.html'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
