"""Annotation library: a common list everyone edits, and a private list per user.

An entry is a named sequence: a label, a feature type and the DNA read 5'->3' on the feature's own
strand. Annotating a map finds every exact occurrence of every entry on both strands (across the
origin on a circular map) and proposes it as a feature, skipping spots already annotated with the
same label and span. Nothing is written by annotating: the proposals are saved through the usual
"save as new version" edit, or, when a user turns on "annotate my uploads", added to their upload.

The common list can be seeded from the bank: `suggest_from_bank` ranks feature sequences by how
many maps carry them, so the parts the lab actually uses come first.
"""

from __future__ import annotations

import re
import sqlite3
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from .sequtil import revcomp

TYPES = (
    "CDS",
    "gene",
    "promoter",
    "terminator",
    "rep_origin",
    "misc_feature",
    "polyA_signal",
    "RBS",
    "protein_bind",
    "regulatory",
    "enhancer",
    "LTR",
    "intron",
    "exon",
    "5'UTR",
    "3'UTR",
    "sig_peptide",
    "primer_bind",
    "mobile_element",
    "misc_recomb",
    "repeat_region",
    "tag",
)
SUGGEST_TYPES = set(TYPES) - {"primer_bind", "exon", "intron", "gene"}
MIN_LEN, MAX_LEN = 15, 10_000
# labels that name a position or a cloning step rather than a part
_GENERIC = re.compile(
    r"^(\d+|insert|vector|fragment.*|feature|new feature|misc_feature|region|selection|"
    r"unnamed|untitled|source|.*mismatch.*|editing history.*|exported)$",
    re.IGNORECASE,
)

SCHEMA = """
CREATE TABLE IF NOT EXISTS entries (
  id INTEGER PRIMARY KEY,
  owner TEXT NOT NULL,                 -- '' = the common list; otherwise the user's private list
  label TEXT NOT NULL, type TEXT NOT NULL, seq TEXT NOT NULL, note TEXT NOT NULL DEFAULT '',
  added_by TEXT NOT NULL, added TEXT NOT NULL, maps INTEGER,
  UNIQUE (owner, seq, label));
CREATE TABLE IF NOT EXISTS settings (
  user TEXT PRIMARY KEY, annotate_on_upload INTEGER NOT NULL DEFAULT 0);
"""


def clean_seq(s: str) -> str:
    s = re.sub(r"\s", "", s).upper()
    if set(s) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("an annotation sequence may only contain DNA letters")
    if not MIN_LEN <= len(s) <= MAX_LEN:
        raise ValueError(f"an annotation must be {MIN_LEN}–{MAX_LEN:,} bases long")
    return s


class Library:
    def __init__(self, state: Path):
        state.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(state / "annotations.sqlite", check_same_thread=False)
        self.db.row_factory = sqlite3.Row
        self.db.executescript(SCHEMA)

    # ── the lists ─────────────────────────────────────────────────────────────
    def entries(self, user: str) -> list[dict]:
        """The common list and this user's private list, each row marked with its scope."""
        rows = self.db.execute(
            "SELECT * FROM entries WHERE owner IN ('', ?) ORDER BY owner != '', label COLLATE NOCASE",
            (user,),
        ).fetchall()
        return [
            {
                **dict(r),
                "scope": "common" if r["owner"] == "" else "mine",
                "length": len(r["seq"]),
            }
            for r in rows
        ]

    def add(
        self,
        user: str,
        scope: str,
        label: str,
        type_: str,
        seq: str,
        note: str = "",
        maps: int | None = None,
    ) -> int:
        label = label.strip()
        if not label:
            raise ValueError("give the annotation a label")
        if type_ not in TYPES:
            raise ValueError(f"type must be one of: {', '.join(TYPES)}")
        owner = "" if scope == "common" else user
        try:
            with self.db:
                cur = self.db.execute(
                    "INSERT INTO entries (owner, label, type, seq, note, added_by, added, maps) VALUES (?,?,?,?,?,?,?,?)",
                    (
                        owner,
                        label,
                        type_,
                        clean_seq(seq),
                        note.strip(),
                        user,
                        datetime.now().isoformat(timespec="minutes"),
                        maps,
                    ),
                )
        except sqlite3.IntegrityError:
            raise ValueError(
                f"'{label}' with that sequence is already in the {'common' if scope == 'common' else 'private'} list"
            )
        return cur.lastrowid

    def _own(self, user: str, entry_id: int) -> sqlite3.Row:
        r = self.db.execute(
            "SELECT * FROM entries WHERE id = ?", (entry_id,)
        ).fetchone()
        if not r or r["owner"] not in (
            "",
            user,
        ):  # common: anyone; private: only its owner
            raise LookupError("no such annotation")
        return r

    def update(
        self,
        user: str,
        entry_id: int,
        label: str | None = None,
        type_: str | None = None,
        note: str | None = None,
    ) -> None:
        r = self._own(user, entry_id)
        if type_ is not None and type_ not in TYPES:
            raise ValueError(f"type must be one of: {', '.join(TYPES)}")
        with self.db:
            self.db.execute(
                "UPDATE entries SET label = ?, type = ?, note = ? WHERE id = ?",
                (
                    (label or r["label"]).strip(),
                    type_ or r["type"],
                    r["note"] if note is None else note,
                    entry_id,
                ),
            )

    def delete(self, user: str, entry_id: int) -> None:
        self._own(user, entry_id)
        with self.db:
            self.db.execute("DELETE FROM entries WHERE id = ?", (entry_id,))

    # ── per-user setting ──────────────────────────────────────────────────────
    def on_upload(self, user: str) -> bool:
        r = self.db.execute(
            "SELECT annotate_on_upload FROM settings WHERE user = ?", (user,)
        ).fetchone()
        return bool(r and r[0])

    def set_on_upload(self, user: str, on: bool) -> None:
        with self.db:
            self.db.execute(
                "INSERT INTO settings (user, annotate_on_upload) VALUES (?, ?) "
                "ON CONFLICT(user) DO UPDATE SET annotate_on_upload = excluded.annotate_on_upload",
                (user, int(on)),
            )


# ── finding entries in a map ──────────────────────────────────────────────────


def _finds(hay: str, needle: str):
    i = hay.find(needle)
    while i != -1:
        yield i
        i = hay.find(needle, i + 1)


def annotate(m, entries: list[dict]) -> list[dict]:
    """Every exact occurrence of every entry, as features (1-based, end may wrap), new ones only."""
    if not m.seq:
        return []
    seq, L = m.seq.upper(), len(m.seq)
    hay = seq + seq if m.circular else seq
    have = {(f.label.lower(), f.start, f.end) for f in m.features}
    have_span = {(f.start, f.end) for f in m.features}
    out, seen = [], set()
    for e in entries:
        s = e["seq"]
        if len(s) > L:
            continue
        for strand, needle in ((1, s), (-1, revcomp(s))):
            if strand == -1 and needle == s:  # palindromic: one strand is enough
                continue
            for i in _finds(hay, needle):
                if i >= L:  # the second copy of the doubled sequence: already seen
                    break
                start0, end0 = i, (i + len(s)) % L or L
                key = (start0, end0, strand, e["label"].lower())
                if key in seen or (e["label"].lower(), start0, end0) in have:
                    continue
                seen.add(key)
                out.append(
                    {
                        "entry": e["id"],
                        "label": e["label"],
                        "type": e["type"],
                        "scope": e.get("scope", ""),
                        "start": start0 + 1,
                        "end": end0,
                        "strand": strand,
                        "length": len(s),
                        "same_span_annotated": (start0, end0) in have_span,
                    }
                )

    # two variants of one part (a 589-bp and a 281-bp 'pBR322 ori'): keep the longer, drop the one inside it
    def inside(a, b):  # a lies within b, on the unwrapped circle
        a0, b0 = a["start"] - 1, b["start"] - 1
        for k in (0, L, -L):
            if b0 <= a0 + k and a0 + k + a["length"] <= b0 + b["length"]:
                return True
        return False

    out = [
        p
        for p in out
        if not any(
            q is not p
            and q["label"].lower() == p["label"].lower()
            and q["length"] > p["length"]
            and inside(p, q)
            for q in out
        )
    ]
    return sorted(out, key=lambda p: (p["start"], -p["length"]))


# ── seeding the common list from the bank ─────────────────────────────────────


def _feature_seq(m, f) -> str:
    L = len(m.seq)
    s = (
        (m.seq + m.seq)[f.start : f.end + L]
        if f.end < f.start
        else m.seq[f.start : f.end]
    )
    s = s.upper()
    return revcomp(s) if f.strand and f.strand < 0 else s


def suggest_from_bank(
    maps, known: set[str], min_maps: int = 3, limit: int = 200
) -> list[dict]:
    """Feature sequences found in at least `min_maps` maps, most widespread first, not yet listed."""
    where: dict[str, set] = defaultdict(set)
    labels: dict[str, Counter] = defaultdict(Counter)
    types: dict[str, Counter] = defaultdict(Counter)
    for m in maps:
        if not m.seq:
            continue
        for f in m.features:
            if (
                f.type not in SUGGEST_TYPES
                or not f.label
                or _GENERIC.match(f.label.strip())
            ):
                continue
            s = _feature_seq(m, f)
            if not MIN_LEN <= len(s) <= MAX_LEN or set(s) - set("ACGT"):
                continue
            where[s].add(m.plasmid_id or m.path)
            labels[s][re.sub(r"\s+", " ", f.label.replace("\\", " ")).strip()] += (
                1  # as the cleaner's unescape_labels
            )
            types[s][f.type] += 1
    out = []
    for s, ms in where.items():
        if len(ms) < min_maps or s in known:
            continue
        top = labels[s].most_common(3)
        out.append(
            {
                "label": top[0][0],
                "also_called": [x for x, _ in top[1:]],
                "type": types[s].most_common(1)[0][0],
                "seq": s,
                "length": len(s),
                "maps": len(ms),
            }
        )
    out.sort(key=lambda c: (-c["maps"], -c["length"]))
    return out[:limit]
