"""Clone verification: compare sequencing results with a plasmid map.

Two modes, chosen by query length:
- read:  a Sanger read (or any partial sequence) -> local alignment against the map,
         both strands, across the origin (map doubled when circular).
- whole: a whole-plasmid consensus (nanopore services etc.) -> rotated to a shared
         k-mer anchor on the right strand, then globally aligned. Indels inside
         homopolymer runs are flagged, since nanopore consensus often gets those wrong.

IUPAC codes count as matches when they are compatible (N in an in-silico map or a
read matches anything). Reads with qualities (.ab1, FASTQ) are Mott-trimmed first.
"""

from __future__ import annotations

import io
import itertools
from dataclasses import dataclass, field
from functools import lru_cache

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

from ..formats import MapRecord, decode, sniff
from ..sequtil import IUPAC, compatible, revcomp

MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2.0, -3.0, -5.0, -2.0
WHOLE_FRACTION = 0.8  # query >= 80% of map length -> whole-plasmid mode
ANCHOR_K = 25
CONTEXT = 12
MAX_QUERIES = 20  # raw nanopore FASTQ can hold thousands of reads; compare a consensus instead


@dataclass
class Query:
    name: str
    seq: str
    qual: list[int] | None = None
    trimmed: tuple[int, int] = (0, 0)  # bases removed from (left, right)
    source_format: str = ""
    # .ab1 only: the four processed trace channels and, per called base (after trimming), the
    # sample index of its peak. None when the file has no peak table (PLOC2/PLOC1).
    trace: dict[str, list[int]] | None = None
    peaks: list[int] | None = None


@dataclass
class Discrepancy:
    kind: str  # substitution | insertion | deletion
    ref_pos: (
        int  # 1-based position on the map (for insertions: base before the insertion)
    )
    ref_bases: str
    query_bases: str
    query_pos: int  # 1-based on the oriented, trimmed query
    homopolymer: bool = False
    quality: int | None = None
    features: list[str] = field(default_factory=list)
    context_ref: str = ""
    context_query: str = ""


@dataclass
class VerifyResult:
    query_name: str
    mode: str  # read | whole
    strand: str  # + | -
    ref_start: int  # 1-based
    ref_end: (
        int  # 1-based, inclusive; < ref_start when the aligned region wraps the origin
    )
    query_length: int
    aligned_columns: int
    identity: float
    score: float
    trimmed: tuple[int, int]
    discrepancies: list[Discrepancy]
    warnings: list[str] = field(default_factory=list)
    # Gapless aligned segments as [map_pos, query_pos, length], both 1-based; query_pos is on the
    # oriented, trimmed query (unrotated). Either may wrap: base k is at (pos - 1 + k) mod length.
    blocks: list[list[int]] = field(default_factory=list)
    # The trimmed query in map orientation, unrotated: what `blocks` count query_pos on. Lets the app
    # draw the reads under the reference, base by base.
    query: str = ""

    @property
    def verdict(self) -> str:
        real = [
            d
            for d in self.discrepancies
            if not (self.mode == "whole" and d.homopolymer)
        ]
        if self.identity < 0.9:
            return "no match"
        return "perfect match" if not real else f"{len(real)} difference(s)"


# ── query loading ──────────────────────────────────────────────────────────────


def mott_trim(qual: list[int], cutoff: float = 0.05) -> tuple[int, int]:
    """Mott's modified algorithm (as in phred/Biopython abi-trim): best-scoring interval."""
    best, best_start, best_end = 0.0, 0, 0
    run, run_start = 0.0, 0
    for i, q in enumerate(qual):
        run += cutoff - 10 ** (-q / 10)
        if run <= 0:
            run, run_start = 0.0, i + 1
        elif run > best:
            best, best_start, best_end = run, run_start, i + 1
    return best_start, best_end


def load_queries(filename: str, data: bytes, trim: bool = True) -> list[Query]:
    lower = filename.lower()
    if lower.endswith((".ab1", ".abi", ".ab!")) or data[:4] == b"ABIF":
        rec = SeqIO.read(io.BytesIO(data), "abi")
        recs, fmt = [rec], "abi"
        trace, peaks = _abi_trace(rec)
    else:
        text, _ = decode(data)
        t = text.lstrip()
        if t.startswith("@"):
            recs = list(itertools.islice(SeqIO.parse(io.StringIO(text), "fastq"), MAX_QUERIES + 1))
            fmt = "fastq"
        else:
            fmt = sniff(data, text)
            if fmt == "raw":
                return [
                    Query(filename, "".join(text.split()).upper(), source_format="raw")
                ]
            if fmt not in ("genbank", "fasta"):
                raise ValueError(f"{filename}: unrecognised sequence format")
            recs = list(SeqIO.parse(io.StringIO(text), fmt))
    if len(recs) > MAX_QUERIES:
        raise ValueError(f"{filename}: more than {MAX_QUERIES} sequences; upload the consensus, not raw reads")
    out = []
    for r in recs:
        seq = str(r.seq).upper()
        qual = r.letter_annotations.get("phred_quality")
        trimmed = (0, 0)
        pk = peaks if fmt == "abi" else None
        if qual and trim:
            s, e = mott_trim(list(qual))
            if e - s >= 20:
                trimmed = (s, len(seq) - e)
                seq, qual = seq[s:e], list(qual[s:e])
                pk = pk[s:e] if pk else None
        out.append(
            # an .ab1 is named after its file: the trace's own sample name is the provider's long number, or
            # Biopython's "<unknown id>"
            Query(filename.rsplit("/", 1)[-1].rsplit(".", 1)[0] if fmt == "abi" else r.id or filename, seq, list(qual) if qual else None, trimmed, fmt,
                  trace if pk else None, pk)
        )
    return out


def _abi_trace(rec) -> tuple[dict[str, list[int]] | None, list[int] | None]:
    """Processed traces (DATA9-12, in the order given by FWO_1) and the peak index of each base."""
    raw = rec.annotations.get("abif_raw", {})
    order = raw.get("FWO_1", b"GATC")
    order = order.decode() if isinstance(order, bytes) else str(order)
    chans = [raw.get(f"DATA{9 + i}") for i in range(4)]
    peaks = raw.get("PLOC2") or raw.get("PLOC1")
    if not all(chans) or not peaks or len(peaks) != len(rec.seq):
        return None, None
    return {b: list(c) for b, c in zip(order, chans)}, list(peaks)


_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def trace_window(q: Query, strand: str, query_pos: int, flank: int = 6) -> dict | None:
    """The chromatogram around one base of an aligned read, in MAP orientation.

    ``query_pos`` is 1-based on the oriented, trimmed read (as in Discrepancy.query_pos). For a read
    that aligned on the minus strand the window is reversed and its channels complemented, so it
    reads like the map. Returns None when the file has no peak table.
    """
    if not q.peaks or not q.trace:
        return None
    n = len(q.peaks)
    i = query_pos - 1 if strand == "+" else n - query_pos
    i = max(0, min(n - 1, i))
    lo, hi = max(0, i - flank), min(n - 1, i + flank)
    s0 = max(0, q.peaks[lo] - (q.peaks[min(lo + 1, n - 1)] - q.peaks[lo]) // 2)
    s1 = q.peaks[hi] + (q.peaks[hi] - q.peaks[max(hi - 1, 0)]) // 2 + 1
    chans = {b: v[s0:s1] for b, v in q.trace.items()}
    peaks = [p - s0 for p in q.peaks[lo : hi + 1]]
    bases = q.seq[lo : hi + 1]
    focus = i - lo
    if strand == "-":
        width = s1 - s0
        chans = {_COMP.get(b, b): v[::-1] for b, v in chans.items()}
        peaks = [width - 1 - p for p in peaks[::-1]]
        bases = "".join(_COMP.get(b, "N") for b in bases[::-1])
        focus = len(bases) - 1 - focus
    return {"channels": chans, "peaks": peaks, "bases": bases, "focus": focus}


def trace_full(q: Query, strand: str) -> dict | None:
    """The whole chromatogram of a read in MAP orientation, plus its called bases and qualities.

    Base i here is base i+1 of the oriented, trimmed query, as used by VerifyResult.blocks. Only
    the stretch from half a base before the first kept peak to half a base after the last is sent.
    """
    if not q.peaks or not q.trace:
        return None
    pk = q.peaks
    n = len(pk)
    pad = max(1, (pk[-1] - pk[0]) // max(1, n - 1) // 2) if n > 1 else 5
    s0 = max(0, pk[0] - pad)
    s1 = pk[-1] + pad + 1
    chans = {b: v[s0:s1] for b, v in q.trace.items()}
    peaks = [p - s0 for p in pk]
    bases, qual = q.seq, list(q.qual or [])
    if strand == "-":
        width = max(len(v) for v in chans.values())
        chans = {_COMP.get(b, b): v[::-1] for b, v in chans.items()}
        peaks = [width - 1 - p for p in peaks[::-1]]
        bases = revcomp(bases)
        qual = qual[::-1]
    return {"channels": chans, "peaks": peaks, "bases": bases, "qual": qual}


def read_support(results: list[VerifyResult], L: int) -> list[list[dict]]:
    """For each read mode result and each discrepancy: how many reads cover that map position, and
    how many of them report the same change. A change seen in one read of several that cover it is
    more likely a misread; one seen in every covering read is more likely real."""
    def covers(r: VerifyResult, pos: int) -> bool:
        if r.mode == "whole":
            return True
        s, e = r.ref_start, r.ref_end
        return s <= pos <= e if s <= e else pos >= s or pos <= e

    out = []
    for r in results:
        row = []
        for d in r.discrepancies:
            cov = [x for x in results if covers(x, d.ref_pos)]
            same = [
                x
                for x in cov
                if any(
                    y.ref_pos == d.ref_pos and y.kind == d.kind and y.query_bases == d.query_bases
                    for y in x.discrepancies
                )
            ]
            row.append({"covered_by": len(cov), "seen_in": len(same)})
        out.append(row)
    return out


# ── alignment core ─────────────────────────────────────────────────────────────


@lru_cache(maxsize=None)
def _matrix():
    alphabet = "".join(IUPAC)
    m = substitution_matrices.Array(alphabet=alphabet, dims=2)
    for a in alphabet:
        for b in alphabet:
            m[a, b] = MATCH if set(IUPAC[a]) & set(IUPAC[b]) else MISMATCH
    return m


def _aligner(mode: str) -> Align.PairwiseAligner:
    return Align.PairwiseAligner(
        mode=mode,
        substitution_matrix=_matrix(),
        open_gap_score=GAP_OPEN,
        extend_gap_score=GAP_EXTEND,
    )


def _clean(seq: str) -> str:
    return "".join(c if c in IUPAC else "N" for c in seq.upper()).replace("U", "T")


def _best(aligner, target: str, query: str):
    return next(iter(aligner.align(target, query)))


def _homorun(seq: str, left: int, base: str) -> int:
    """Length of the run of `base` spanning the boundary between seq[left] and seq[left + 1]."""
    n, count, i = len(seq), 0, left
    while i >= 0 and seq[i] == base:
        count, i = count + 1, i - 1
    j = left + 1
    while j < n and seq[j] == base:
        count, j = count + 1, j + 1
    return count


def _features_at(ref: MapRecord, pos0: int) -> list[str]:
    hits = []
    for f in ref.features:
        inside = (
            (f.start <= pos0 < f.end)
            if not f.wraps
            else (pos0 >= f.start or pos0 < f.end)
        )
        if inside:
            hits.append(f"{f.label} ({f.type})")
    return hits


def _walk(aln, ref: MapRecord, ref_offset: int, query: str, qual, mode: str):
    """Turn an alignment into discrepancies. Target coordinates are shifted by ref_offset mod L."""
    L = ref.length
    target = aln.target
    cols_t: list[str] = []
    cols_q: list[str] = []
    events: list[tuple[int, Discrepancy]] = []
    matches = 0
    coords = aln.coordinates
    for k in range(coords.shape[1] - 1):
        t0, t1 = int(coords[0, k]), int(coords[0, k + 1])
        q0, q1 = int(coords[1, k]), int(coords[1, k + 1])
        if t1 > t0 and q1 > q0:
            for i in range(t1 - t0):
                a, b = target[t0 + i], query[q0 + i]
                if compatible(a, b):
                    matches += 1
                else:
                    pos0 = (t0 + i + ref_offset) % L
                    events.append(
                        (
                            len(cols_t),
                            Discrepancy(
                                "substitution",
                                pos0 + 1,
                                a,
                                b,
                                q0 + i + 1,
                                quality=qual[q0 + i] if qual else None,
                                features=_features_at(ref, pos0),
                            ),
                        )
                    )
                cols_t.append(a)
                cols_q.append(b)
        elif t1 > t0:
            bases = target[t0:t1]
            pos0 = (t0 + ref_offset) % L
            homo = len(set(bases)) == 1 and _homorun(ref.seq, pos0 - 1, bases[0]) >= 3
            events.append(
                (
                    len(cols_t),
                    Discrepancy(
                        "deletion",
                        pos0 + 1,
                        bases,
                        "-" * len(bases),
                        q0 + 1,
                        homopolymer=homo,
                        features=_features_at(ref, pos0),
                    ),
                )
            )
            cols_t.extend(bases)
            cols_q.extend("-" * len(bases))
        elif q1 > q0:
            bases = query[q0:q1]
            pos0 = (t0 - 1 + ref_offset) % L
            homo = len(set(bases)) == 1 and _homorun(ref.seq, pos0, bases[0]) >= 2
            q = min(qual[q0:q1]) if qual else None
            events.append(
                (
                    len(cols_t),
                    Discrepancy(
                        "insertion",
                        pos0 + 1,
                        "-" * len(bases),
                        bases,
                        q0 + 1,
                        homopolymer=homo,
                        quality=q,
                        features=_features_at(ref, pos0),
                    ),
                )
            )
            cols_t.extend("-" * len(bases))
            cols_q.extend(bases)
    for col, d in events:
        lo, hi = max(0, col - CONTEXT), col + CONTEXT + max(len(d.ref_bases), 1)
        d.context_ref = "".join(cols_t[lo:hi])
        d.context_query = "".join(cols_q[lo:hi])
    return [d for _, d in events], matches, len(cols_t)


def _anchor(ref_seq: str, query: str) -> tuple[str, int, int, int] | None:
    """Find strand and a shared exact k-mer: (strand, ref_pos, query_pos, votes)."""
    L, step = len(ref_seq), max(1, len(ref_seq) // 60)
    ref2 = ref_seq + ref_seq[:ANCHOR_K]
    best = None
    for strand, q in (("+", query), ("-", revcomp(query))):
        q2 = q + q[:ANCHOR_K]
        index: dict[str, int] = {}
        dup = set()
        for i in range(len(q)):
            km = q2[i : i + ANCHOR_K]
            if km in index:
                dup.add(km)
            index[km] = i
        votes: dict[int, list[int]] = {}
        for r in range(0, L, step):
            km = ref2[r : r + ANCHOR_K]
            if km in index and km not in dup:
                off = (index[km] - r) % len(q)
                votes.setdefault(off, []).append(r)
        if votes:
            off, rs = max(votes.items(), key=lambda kv: len(kv[1]))
            r0 = rs[0]
            cand = (strand, r0, (r0 + off) % len(q), len(rs))
            if best is None or cand[3] > best[3]:
                best = cand
    return best


def verify(ref: MapRecord, q: Query, mode: str | None = None) -> VerifyResult:
    if not ref.seq:
        raise ValueError("map has no sequence")
    seq = _clean(q.seq)
    L = ref.length
    warnings: list[str] = []
    if mode is None:
        mode = "whole" if len(seq) >= WHOLE_FRACTION * L else "read"
    elif mode == "whole" and len(seq) < WHOLE_FRACTION * L:
        # a Sanger read aligned end to end against the whole map reads as one long deletion
        warnings.append(
            f"{len(seq):,} bp is too short for a whole-plasmid consensus ({L:,} bp map); checked as a read"
        )
        mode = "read"
    if mode == "whole":
        anchor = _anchor(ref.seq, seq) if ref.circular else None
        if anchor is None and ref.circular:
            warnings.append(
                "no shared 25-mer between consensus and map; fell back to read mode"
            )
            mode = "read"
    if mode == "whole":
        if ref.circular:
            strand, r0, q0, _ = anchor
            oq = seq if strand == "+" else revcomp(seq)
            qual = q.qual if strand == "+" or not q.qual else q.qual[::-1]
            target = ref.seq[r0:] + ref.seq[:r0]
            rot = q0
            oq = oq[q0:] + oq[:q0]
            if qual:
                qual = qual[q0:] + qual[:q0]
            offset = r0
        else:
            alns = [
                (s, _best(_aligner("global"), ref.seq, s2))
                for s, s2 in (("+", seq), ("-", revcomp(seq)))
            ]
            strand = max(alns, key=lambda x: x[1].score)[0]
            oq = seq if strand == "+" else revcomp(seq)
            qual = q.qual if strand == "+" or not q.qual else q.qual[::-1]
            target, offset = ref.seq, 0
        aln = _best(_aligner("global"), target, oq)
        if abs(len(oq) - L) > 0.05 * L:
            warnings.append(
                f"consensus length {len(oq)} differs from map length {L} by >5%"
            )
    else:
        target = ref.seq + ref.seq if ref.circular else ref.seq
        a = _aligner("local")
        cands = []
        for strand, s in (("+", seq), ("-", revcomp(seq))):
            cands.append((strand, s, _best(a, target, s)))
        strand, oq, aln = max(cands, key=lambda x: x[2].score)
        qual = q.qual if strand == "+" or not q.qual else q.qual[::-1]
        offset = 0
        aligned_q = int(aln.coordinates[1, -1] - aln.coordinates[1, 0])
        if aligned_q < 0.7 * len(oq):
            warnings.append(
                f"only {aligned_q} of {len(oq)} query bases align to this map"
            )
    discrepancies, matches, ncols = _walk(aln, ref, offset, oq, qual, mode)
    rot = rot if mode == "whole" and ref.circular else 0
    c = aln.coordinates
    blocks = [
        [(int(c[0, k]) + offset) % L + 1, (int(c[1, k]) + rot) % len(oq) + 1, int(c[0, k + 1] - c[0, k])]
        for k in range(c.shape[1] - 1)
        if c[0, k + 1] > c[0, k] and c[1, k + 1] > c[1, k]
    ]
    t_start = int(aln.coordinates[0, 0]) + offset
    t_end = int(aln.coordinates[0, -1]) + offset
    identity = matches / ncols if ncols else 0.0
    if q.qual is None and mode == "read":
        warnings.append("no base qualities: read ends were not trimmed")
    return VerifyResult(
        q.name,
        mode,
        strand,
        t_start % L + 1,
        (t_end - 1) % L + 1,
        len(seq),
        ncols,
        identity,
        float(aln.score),
        q.trimmed,
        discrepancies,
        warnings,
        blocks,
        seq if strand == "+" else revcomp(seq),
    )
