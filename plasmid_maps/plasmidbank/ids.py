"""PMC / EC identifier normalisation.

Identity comes from the filename; the LOCUS name is only a fallback (it is often
`Exported`, a free-text name, or a different PMC number than the file).
"""

from __future__ import annotations

import re
import unicodedata
from dataclasses import dataclass, field

PMC_RE = re.compile(r"PMC[-_ ]?0*(\d{1,5})(?!\d)", re.IGNORECASE)
EC_RE = re.compile(r"(?<![A-Za-z])EC[-_ ]?(\d{5})(?!\d)", re.IGNORECASE)


@dataclass
class IdResult:
    plasmid_id: str | None
    kind: str  # "PMC" | "EC" | "NONE"
    source: str  # "filename" | "locus" | "none"
    xrefs: list[str] = field(default_factory=list)


def _matches(text: str) -> list[tuple[int, str, str]]:
    text = unicodedata.normalize("NFC", text)
    found = [
        (m.start(), "PMC", f"PMC-{int(m.group(1)):05d}") for m in PMC_RE.finditer(text)
    ]
    found += [(m.start(), "EC", f"EC{m.group(1)}") for m in EC_RE.finditer(text)]
    return sorted(found)


def find_ids(text: str) -> list[str]:
    seen: list[str] = []
    for _, _, pid in _matches(text):
        if pid not in seen:
            seen.append(pid)
    return seen


def identify(filename: str, locus_name: str | None = None) -> IdResult:
    stem = filename.rsplit("/", 1)[-1]
    file_hits = _matches(stem)
    locus_ids = find_ids(locus_name or "")
    if file_hits:
        _, kind, pid = file_hits[0]
        xrefs = [i for i in find_ids(stem) + locus_ids if i != pid]
        return IdResult(pid, kind, "filename", list(dict.fromkeys(xrefs)))
    if locus_ids:
        pid = locus_ids[0]
        kind = "PMC" if pid.startswith("PMC") else "EC"
        return IdResult(pid, kind, "locus", locus_ids[1:])
    return IdResult(None, "NONE", "none")
