"""Sequence analysis helpers for the map view: ORFs, translation, multiple alignment (MAFFT)."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio.Seq import Seq

from ..formats import MapRecord
from ..sequtil import revcomp
from .align import _anchor

STOPS = {"TAA", "TAG", "TGA"}


def find_orfs(
    ref: MapRecord, min_aa: int = 100, starts: tuple[str, ...] = ("ATG",)
) -> list[dict]:
    """ATG..stop ORFs on both strands, across the origin on circular maps; nested starts collapsed."""
    seq, L = ref.seq.upper(), ref.length
    out = []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        s2 = s + s if ref.circular else s
        for frame in range(3):
            i, open_at = frame, None
            limit = len(s2) - 2
            while i < limit:
                cod = s2[i : i + 3]
                if open_at is None and cod in starts and i < L:
                    open_at = i
                elif open_at is not None and cod in STOPS:
                    n_aa = (i - open_at) // 3
                    if n_aa >= min_aa and i + 3 - open_at <= L:
                        a, b = open_at, i + 3  # strand coordinates, end exclusive
                        if strand == "+":
                            start, end = a % L + 1, (b - 1) % L + 1
                        else:
                            start, end = (L - b) % L + 1, (L - a - 1) % L + 1
                        out.append(
                            {
                                "strand": strand,
                                "start": start,
                                "end": end,
                                "length_aa": n_aa,
                                "frame": frame + 1,
                                "protein": str(Seq(s2[a:b]).translate()).rstrip("*"),
                            }
                        )
                    open_at = None
                i += 3
    seen, uniq = set(), []
    for o in sorted(out, key=lambda o: -o["length_aa"]):
        key = (
            o["strand"],
            o["end"] if o["strand"] == "+" else o["start"],
        )  # same stop = same ORF
        if key not in seen:
            seen.add(key)
            uniq.append(o)
    return sorted(uniq, key=lambda o: o["start"])


def translate_region(ref: MapRecord, start: int, end: int, strand: str = "+") -> dict:
    L = ref.length
    s2 = ref.seq + ref.seq
    length = end - start + 1 if end >= start else end + L - start + 1
    dna = s2[start - 1 : start - 1 + length]
    if strand == "-":
        dna = revcomp(dna)
    usable = len(dna) - len(dna) % 3
    prot = str(Seq(dna[:usable]).translate())
    return {
        "dna_length": len(dna),
        "protein": prot,
        "internal_stops": prot[:-1].count("*"),
        "warning": None
        if len(dna) % 3 == 0
        else f"length {len(dna)} is not a multiple of 3",
    }


def join_features(ref: MapRecord, indices: list[int]) -> dict:
    """Join features (exons) into one coding sequence and check its frame ("CLC
    connects annotated exons"). All must be on one strand; they are joined 5'->3' on that strand,
    so on the minus strand the last exon on the map comes first."""
    if not ref.seq:
        raise ValueError("map has no sequence")
    if len(indices) < 1 or len(set(indices)) != len(indices):
        raise ValueError("pick the features to join, each once")
    feats = []
    for i in indices:
        if not 0 <= i < len(ref.features):
            raise ValueError(f"no feature {i}")
        feats.append(ref.features[i])
    strands = {f.strand for f in feats}
    if len(strands) != 1 or 0 in strands:
        raise ValueError("the features must all be on the same strand (+ or −)")
    strand = strands.pop()
    s = ref.seq
    piece = lambda f: s[f.start:f.end] if f.end >= f.start else s[f.start:] + s[: f.end]
    order = sorted(feats, key=lambda f: f.start, reverse=strand < 0)
    dna = "".join(piece(f) if strand > 0 else revcomp(piece(f)) for f in order)
    usable = len(dna) - len(dna) % 3
    prot = str(Seq(dna[:usable]).translate())
    ends = []  # the frame each part starts in, so a frame shift shows where it happens
    pos = 0
    for f in order:
        ends.append({"label": f.label, "start": f.start + 1, "end": f.end, "length": len(piece(f)),
                     "phase": pos % 3})
        pos += len(piece(f))
    return {"strand": "+" if strand > 0 else "-", "length": len(dna), "multiple_of_3": len(dna) % 3 == 0,
            "starts_with_atg": dna[:3].upper() == "ATG", "ends_with_stop": prot.endswith("*") and len(dna) % 3 == 0,
            "internal_stops": prot.rstrip("*").count("*"), "protein": prot, "parts": ends,
            "locations": [[f.start + 1, f.end] for f in order]}


def mafft_available() -> str | None:
    return shutil.which("mafft")


def multiple_alignment(
    records: list[tuple[str, str, bool]], timeout: int = 300
) -> dict:
    """Align (name, seq, circular) with MAFFT. Circular sequences are first rotated to a shared
    25-mer anchor with the first sequence, and MAFFT adjusts direction (--adjustdirection)."""
    exe = mafft_available()
    if not exe:
        raise RuntimeError(
            "MAFFT is not installed (macOS: brew install mafft; Windows: conda install -c bioconda mafft)"
        )
    if len(records) < 2:
        raise ValueError("need at least two sequences")
    ref_name, ref_seq, _ = records[0]
    prepared, notes = [(ref_name, ref_seq.upper())], []
    for name, seq, circular in records[1:]:
        seq = seq.upper()
        if circular:
            a = _anchor(ref_seq.upper(), seq)
            if a:
                strand, r0, q0, _ = a
                oq = seq if strand == "+" else revcomp(seq)
                oq = oq[q0:] + oq[:q0]
                rot = r0 % len(ref_seq)
                # place the anchor at the same position as in the reference
                oq = oq[-rot:] + oq[:-rot] if rot else oq
                seq = oq
                notes.append(
                    f"{name}: rotated{' and reverse-complemented' if strand == '-' else ''} to match {ref_name}"
                )
        prepared.append((name, seq))
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "in.fa"
        fa.write_text("".join(f">s{i}\n{s}\n" for i, (_, s) in enumerate(prepared)))
        res = subprocess.run(
            [exe, "--auto", "--adjustdirection", "--quiet", str(fa)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    if res.returncode != 0:
        raise RuntimeError(f"MAFFT failed: {res.stderr[:400]}")
    aligned: dict[int, list[str]] = {}
    cur = None
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            tag = line[1:].strip()
            if tag.startswith("_R_"):
                tag = tag[3:]
                notes.append(
                    f"{prepared[int(tag[1:])][0]}: reverse-complemented by MAFFT"
                )
            cur = int(tag[1:])
            aligned[cur] = []
        elif cur is not None:
            aligned[cur].append(line.strip().upper())
    rows = [
        {"name": prepared[i][0], "aligned": "".join(aligned[i])}
        for i in sorted(aligned)
    ]
    ncol = len(rows[0]["aligned"])
    consensus, identical = [], 0
    for c in range(ncol):
        col = [r["aligned"][c] for r in rows]
        bases = [b for b in col if b != "-"]
        top = max(set(bases), key=bases.count) if bases else "-"
        same = len(set(col)) == 1
        identical += same
        consensus.append(top if same else top.lower())
    return {
        "rows": rows,
        "consensus": "".join(consensus),
        "columns": ncol,
        "identity": round(identical / ncol, 4) if ncol else 0.0,
        "notes": notes,
        "engine": "MAFFT --auto --adjustdirection",
    }
