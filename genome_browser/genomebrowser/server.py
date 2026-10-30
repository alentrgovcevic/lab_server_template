"""Genome browser: any genome the lab uploads, with its gene models and extra tracks.

Run from genome_browser/:  ../.venv/bin/python -m genomebrowser.server     (http://127.0.0.1:8766)
Env:  GENOMES_DIR          where genomes live (default ~/GenomeBrowser/genomes)
      GENOMES_AUTH         none (a laptop: everyone is 'local', an admin) | header (behind Caddy,
                           which asks the plasmid bank's /auth/check and passes X-Auth-User/-Role)
      GENOMES_IMPORT_DIR   optional read-only folder admins can import from without uploading
      GENOMES_SEED_DIR     optional folder with a seed.json: built on first start when empty
      GENOMES_UPLOADS      admins (default: only admins add genomes and tracks, so the disk holds one
                           curated copy of each) | everyone
      GENOMES_MAX_UPLOAD_GB, GENOMES_HOST / GENOMES_PORT

Coordinates in and out are 1-based, closed (GFF style).
"""

from __future__ import annotations

import os
import re
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse
from pydantic import BaseModel

from . import imports, store
from .build import TRACK_RE, BuildError, slugify

STATIC = Path(__file__).parent / "static"
DETAIL_SPAN = 2_000_000  # wider than this, tracks come back as density bins
SEQ_SPAN = 2_000  # bases sent with a region when it is at most this wide
MAX_FASTA = 2_000_000
MAX_BINS = 2_000
MAX_FEATURES = 20_000  # per track per region request
AUTH = os.environ.get("GENOMES_AUTH", "none").lower()
UPLOADS = os.environ.get("GENOMES_UPLOADS", "admins").lower()


@asynccontextmanager
async def _lifespan(_app):
    imports.start_worker()
    seed = os.environ.get("GENOMES_SEED_DIR")
    if seed and Path(seed).is_dir():
        imports.seed(Path(seed))
    yield


app = FastAPI(title="Genome browser", lifespan=_lifespan)


# ── who is asking ─────────────────────────────────────────────────────────────


@app.middleware("http")
async def _gate(request: Request, call_next):
    if request.url.path == "/healthz":
        return await call_next(request)
    if AUTH == "header":
        user = request.headers.get("x-auth-user", "")
        role = request.headers.get("x-auth-role", "user")
        if not user:
            return JSONResponse({"detail": "not signed in"}, status_code=401)
    else:
        user, role = "local", "admin"
    origin = request.headers.get("origin")
    if request.method not in ("GET", "HEAD", "OPTIONS") and origin:
        if origin.split("://", 1)[-1].rstrip("/") != request.headers.get("host", ""):
            return JSONResponse(
                {"detail": "cross-site request refused"}, status_code=403
            )
    request.state.user, request.state.admin = user, role == "admin"
    return await call_next(request)


def _may_add(request: Request) -> None:
    if not (request.state.admin or UPLOADS == "everyone"):
        raise HTTPException(403, "only admins add genomes and tracks here; ask one to add what you need")


def _may_change(request: Request, *owners: str) -> None:
    if not (request.state.admin or request.state.user in owners):
        raise HTTPException(
            403, "only the person who added it, or an admin, can do this"
        )


# ── pages and listing ─────────────────────────────────────────────────────────


@app.get("/")
def index():
    return FileResponse(STATIC / "index.html", headers={"Cache-Control": "no-cache"})


@app.get("/healthz")
def healthz():
    return {
        "ok": True,
        "genomes": len(list(store.genomes_dir().glob("*/genome.sqlite"))),
    }


@app.get("/api/me")
def me(request: Request):
    return {
        "user": request.state.user,
        "admin": request.state.admin,
        "can_add": request.state.admin or UPLOADS == "everyone",
        "chunk": imports.CHUNK,
        "max_upload": imports.max_upload(),
        "server_import": bool(imports.import_dir()),
        "disk": store.disk(),
    }


@app.get("/api/genomes")
def genomes():
    return store.list_genomes()


@app.delete("/api/genomes/{slug}")
def delete_genome(slug: str, request: Request):
    m = store.meta(slug)
    if not m.get("added"):
        raise HTTPException(404, f"no genome {slug}")
    _may_change(request, m.get("added_by", ""))
    store.delete_genome(slug)
    return {"deleted": slug}


# ── reading one genome ────────────────────────────────────────────────────────


def _db(slug: str):
    c = store.connect(slug)
    if not c:
        raise HTTPException(404, f"no genome {slug}")
    return c


def _seq_row(c, name: str):
    r = c.execute("SELECT * FROM seq WHERE name = ?", (name,)).fetchone()
    if not r:
        raise HTTPException(404, f"no sequence {name}")
    return r


def fetch_bases(slug: str, c, name: str, start: int, end: int) -> str:
    s = _seq_row(c, name)
    start, end = max(1, start), min(end, s["length"])
    if end < start:
        return ""
    lb, lw, s0 = s["linebases"], s["linewidth"], start - 1
    with open(store.genome_path(slug) / "genome.fa", "rb") as fh:
        fh.seek(s["offset"] + s0 // lb * lw + s0 % lb)
        raw = fh.read((end - s0) // lb * lw + lw + 1)
    return raw.replace(b"\n", b"").replace(b"\r", b"")[: end - s0].decode().upper()


def revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))


def _overlapping(c, num: int, start: int, end: int, where: str = "", args: tuple = ()):
    """Features overlapping [start, end] on sequence `num`. The R*Tree is a float32 prefilter;
    the exact test on feat.start/end follows."""
    return c.execute(
        f"""SELECT f.* FROM feat_rt r JOIN feat f ON f.id = r.id
            WHERE r.c0 = ? AND r.c1 = ? AND r.s <= ? AND r.e >= ? AND f.start <= ? AND f."end" >= ? {where}""",
        (num, num, end, start, end, start, *args),
    ).fetchall()


@app.get("/api/g/{slug}/seqs")
def seqs(slug: str):
    c = _db(slug)
    return {
        **store.meta(slug),
        "seqs": [
            dict(r) for r in c.execute("SELECT num, name, length FROM seq ORDER BY num")
        ],
        "tracks": [dict(r) for r in c.execute("SELECT * FROM tracks ORDER BY rowid")],
    }


@app.get("/api/g/{slug}/search")
def search(slug: str, q: str, limit: int = 20):
    limit = max(1, min(limit, 100))
    q = q.strip()
    if len(q) < 2:
        return []
    c = _db(slug)
    cols = 'f.fid, f.name, f.type, f.track, f.seqid, f.start, f."end", f.strand'
    out = [
        dict(r)
        for r in c.execute(
            f"""SELECT {cols}, g.defline FROM names n JOIN feat f ON f.id = n.feat_id
            LEFT JOIN gene_info g ON g.gene = f.fid
            WHERE n.name LIKE ? ORDER BY length(n.name), n.name LIMIT ?""",
            (q + "%", limit),
        )
    ]
    if (
        not out and len(q) >= 5
    ):  # a number without its prefix, e.g. a mutant line "30055169"
        out = [
            dict(r)
            for r in c.execute(
                f"""SELECT {cols}, NULL AS defline FROM names n JOIN feat f ON f.id = n.feat_id
                WHERE n.name LIKE ? LIMIT ?""",
                (f"%{q}", limit),
            )
        ]
    if len(out) < limit and len(q) >= 3:  # words in the descriptions
        seen = {r["fid"] for r in out}
        for r in c.execute(
            f"""SELECT {cols}, g.defline FROM gene_info g JOIN feat f ON f.fid = g.gene
                WHERE g.defline LIKE ? LIMIT ?""",
            (f"%{q}%", limit - len(out)),
        ):
            if r["fid"] not in seen:
                out.append(dict(r))
                seen.add(r["fid"])
    one, seen_ids = [], set()  # one row per ID: the feature card lists every occurrence
    for r in out:
        if (r["track"], r["fid"]) not in seen_ids:
            seen_ids.add((r["track"], r["fid"]))
            one.append(r)
    return one


@app.get("/api/g/{slug}/region")
def region(slug: str, seq: str, start: int, end: int):
    c = _db(slug)
    s = _seq_row(c, seq)
    start, end = max(1, start), min(end, s["length"])
    tracks = [
        dict(r) for r in c.execute("SELECT name, kind FROM tracks ORDER BY rowid")
    ]
    out = {
        "seq": seq,
        "start": start,
        "end": end,
        "length": s["length"],
        "track_list": tracks,
    }
    if end - start + 1 > DETAIL_SPAN:
        b = int(c.execute("SELECT value FROM meta WHERE key = 'bin'").fetchone()[0])
        first, last = (start - 1) // b, (end - 1) // b
        k = max(1, -(-(last - first + 1) // MAX_BINS))  # merge k bins into one
        out["bin"] = b * k
        out["density"] = {}
        for r in c.execute(
            """SELECT track, bin / ? AS bb, sum(n) AS n FROM density
               WHERE seqid = ? AND bin BETWEEN ? AND ? GROUP BY track, bb""",
            (k, seq, first, last),
        ):
            out["density"].setdefault(r["track"], []).append([r["bb"], r["n"]])
        return out
    keys = ("start", "end", "strand", "type", "fid", "name", "parent")
    out["tracks"] = {t["name"]: [] for t in tracks}
    out["truncated"] = []
    genes = {}
    for r in _overlapping(c, s["num"], start, end):
        lst = out["tracks"].setdefault(r["track"], [])
        if len(lst) >= MAX_FEATURES:
            if r["track"] not in out["truncated"]:
                out["truncated"].append(r["track"])
            continue
        lst.append({k: r[k] for k in keys if r[k] is not None})
        if r["type"] == "gene" and r["fid"]:
            genes[r["fid"]] = lst[-1]
    ids = list(genes)
    for i in range(0, len(ids), 500):
        chunk = ids[i : i + 500]
        for g in c.execute(
            f"SELECT gene, defline FROM gene_info WHERE gene IN ({','.join('?' * len(chunk))})",
            chunk,
        ):
            genes[g["gene"]]["defline"] = g["defline"]
    if end - start + 1 <= SEQ_SPAN:
        out["bases"] = fetch_bases(slug, c, seq, start, end)
    return out


def _first_transcript_parts(c, gene_fid: str) -> tuple[list[dict], list]:
    tx = c.execute(
        """SELECT fid, start, "end" FROM feat WHERE parent = ? AND type IN ('mRNA','transcript')
           ORDER BY start, fid""",
        (gene_fid,),
    ).fetchall()
    out = []
    for t in tx:
        parts = c.execute(
            'SELECT type, start, "end" FROM feat WHERE parent = ? ORDER BY start',
            (t["fid"],),
        ).fetchall()
        out.append(
            {
                "fid": t["fid"],
                "start": t["start"],
                "end": t["end"],
                "parts": [dict(p) for p in parts],
            }
        )
    return out, (out[0]["parts"] if out else [])


def _classify(parts, pos: int) -> str | None:
    if not parts:
        return None
    inside = lambda t: any(
        p["type"] == t and p["start"] <= pos <= p["end"] for p in parts
    )  # noqa: E731
    return "CDS" if inside("CDS") else "exon (UTR)" if inside("exon") else "intron"


_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {
    a + b + c: _AA[i * 16 + j * 4 + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}


IUPAC = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
         "H": "ACT", "V": "ACG", "N": "ACGT"}


def _codon(c: str) -> str:
    if c in CODON:
        return CODON[c]
    aa = {CODON.get(x + y + z) for x in IUPAC.get(c[0], c[0]) for y in IUPAC.get(c[1], c[1])
          for z in IUPAC.get(c[2], c[2])}
    return aa.pop() if len(aa) == 1 else "X"  # GGN is still G


def translate(cds: str) -> str:
    """Standard code; an ambiguous codon is X unless every reading gives the same amino acid."""
    return "".join(_codon(cds[i : i + 3]) for i in range(0, len(cds) - 2, 3))


@app.get("/api/g/{slug}/transcript/{fid}")
def transcript(slug: str, fid: str):
    """One transcript spliced: mRNA, CDS and protein, and the exon/intron table, all derived from
    the exon and CDS features and the genome. A gene whose CDS hangs directly off it (organelle
    GFF3s) counts as its own transcript. Phase is recomputed from the CDS lengths, assuming the CDS
    starts on a full codon."""
    c = _db(slug)
    t = c.execute(
        """SELECT * FROM feat WHERE fid = ? AND type IN ('mRNA','transcript','gene')
           ORDER BY CASE type WHEN 'gene' THEN 1 ELSE 0 END, id LIMIT 1""",
        (fid,),
    ).fetchone()
    if not t:
        raise HTTPException(404, f"no transcript {fid}")
    kids = c.execute(
        'SELECT type, start, "end" FROM feat WHERE parent = ? AND track = ? ORDER BY start',
        (fid, t["track"]),
    ).fetchall()
    cds = [(k["start"], k["end"]) for k in kids if k["type"] == "CDS"]
    exons = [(k["start"], k["end"]) for k in kids if k["type"] == "exon"] or cds
    if not exons:
        raise HTTPException(404, f"{fid} has no exons or CDS")
    if t["end"] - t["start"] + 1 > MAX_FASTA:
        raise HTTPException(400, f"{fid} spans more than {MAX_FASTA:,} bp")
    minus = t["strand"] == "-"
    lo, hi = min(e[0] for e in exons), max(e[1] for e in exons)
    span = fetch_bases(slug, c, t["seqid"], lo, hi)

    def piece(a: int, b: int) -> str:
        s = span[a - lo : b - lo + 1]
        return revcomp(s) if minus else s

    order = (lambda xs: sorted(xs, reverse=True)) if minus else sorted
    exons, cds = order(exons), order(cds)
    mrna = "".join(piece(a, b) for a, b in exons)
    cds_seq = "".join(piece(a, b) for a, b in cds)
    c0, c1 = (min(a for a, _ in cds), max(b for _, b in cds)) if cds else (0, -1)
    phase, done = {}, 0
    for a, b in cds:
        phase[(a, b)] = (3 - done % 3) % 3
        done += b - a + 1
    rows = []
    for n, (a, b) in enumerate(exons, 1):
        up = max(0, min(b, c0 - 1) - a + 1) if cds else 0  # left of the CDS
        down = max(0, b - max(a, c1 + 1) + 1) if cds else 0  # right of the CDS
        utr5, utr3 = (down, up) if minus else (up, down)
        part = next((p for p in cds if p[0] <= b and p[1] >= a), None)
        rows.append(
            {
                "n": n,
                "start": a,
                "end": b,
                "length": b - a + 1,
                "utr5": utr5 if cds else None,
                "utr3": utr3 if cds else None,
                "cds": (b - a + 1 - utr5 - utr3) if cds else 0,
                "phase": phase[part] if part else None,
            }
        )
    introns = []
    for n, (x, y) in enumerate(zip(exons, exons[1:]), 1):
        a, b = (y[1] + 1, x[0] - 1) if minus else (x[1] + 1, y[0] - 1)
        if b < a:
            continue
        s = piece(a, b)
        introns.append(
            {"n": n, "start": a, "end": b, "length": b - a + 1, "splice": f"{s[:2]}…{s[-2:]}"}
        )
    protein = translate(cds_seq)
    stop = protein.endswith("*")
    return {
        "fid": fid,
        "gene": t["parent"] or (fid if t["type"] == "gene" else None),
        "seqid": t["seqid"],
        "strand": t["strand"],
        "start": t["start"],
        "end": t["end"],
        "exons": rows,
        "introns": introns,
        "mrna": mrna,
        "utr5": sum(r["utr5"] or 0 for r in rows),
        "utr3": sum(r["utr3"] or 0 for r in rows),
        "cds": cds_seq,
        "protein": protein[:-1] if stop else protein,
        "starts_atg": cds_seq[:3] == "ATG",
        "ends_stop": stop,
        "internal_stops": protein[:-1].count("*") if stop else protein.count("*"),
        "cds_whole_codons": len(cds_seq) % 3 == 0,
    }


@app.get("/api/g/{slug}/feature/{fid}")
def feature(slug: str, fid: str, track: str = ""):
    c = _db(slug)
    f = c.execute(
        f"""SELECT * FROM feat WHERE fid = ? {"AND track = ?" if track else ""}
            ORDER BY CASE type WHEN 'gene' THEN 0 ELSE 1 END, id LIMIT 1""",
        (fid, track) if track else (fid,),
    ).fetchone()
    if not f:
        raise HTTPException(404, f"no feature {fid}")
    out = {
        k: f[k]
        for k in (
            "track",
            "seqid",
            "start",
            "end",
            "strand",
            "type",
            "fid",
            "name",
            "parent",
            "attrs",
        )
    }
    kind = c.execute("SELECT kind FROM tracks WHERE name = ?", (f["track"],)).fetchone()
    out["kind"] = kind["kind"] if kind else "features"
    out["n_occurrences"] = c.execute(
        "SELECT count(*) FROM feat WHERE fid = ? AND track = ?", (fid, f["track"])
    ).fetchone()[0]
    g = c.execute("SELECT defline FROM gene_info WHERE gene = ?", (fid,)).fetchone()
    out["defline"] = g["defline"] if g else None
    if f["type"] != "gene":
        return out
    out["transcripts"], first = _first_transcript_parts(c, fid)
    out["family"] = _family(c, fid, out["defline"])
    num = _seq_row(c, f["seqid"])["num"]
    others = _overlapping(
        c,
        num,
        f["start"],
        f["end"],
        "AND f.track != ? ORDER BY f.track, f.start LIMIT 2000",
        (f["track"],),
    )
    out["overlaps"] = [
        {
            "track": r["track"],
            "fid": r["fid"],
            "name": r["name"],
            "type": r["type"],
            "start": r["start"],
            "end": r["end"],
            "strand": r["strand"],
            "where": _classify(first, r["start"]),
        }
        for r in others
    ]
    return out


FAMILY_RE = re.compile(r"\bPTHR\d+:SF\d+\b")


def _family(c, fid: str, defline: str | None) -> dict | None:
    """The genes that share this gene's PANTHER subfamily: its copies in this genome.
    Phytozome deflines read "(1 of 2) PTHR48010:SF104 - …"; the 2 counts the genes of that subfamily."""
    m = FAMILY_RE.search(defline or "")
    if not m:
        return None
    tok, genes, seen = m.group(0), [], set()
    for r in c.execute(
        """SELECT f.fid, f.name, f.track, f.seqid, f.start, f."end", f.strand, g.defline
           FROM gene_info g JOIN feat f ON f.fid = g.gene AND f.type = 'gene'
           WHERE g.defline LIKE ? ORDER BY f.seqid, f.start LIMIT 500""",
        (f"%{tok}%",),
    ):
        if r["fid"] in seen or tok not in FAMILY_RE.findall(r["defline"] or ""):
            continue  # PTHR1:SF10 must not match PTHR1:SF104
        seen.add(r["fid"])
        genes.append({k: r[k] for k in ("fid", "name", "track", "seqid", "start", "end", "strand")})
    return {"id": tok, "genes": genes}


@app.get("/api/g/{slug}/occurrences/{track}/{fid}")
def occurrences(slug: str, track: str, fid: str):
    """Every feature with this ID in one track (e.g. all insertions of one mutant line), each with
    the gene it falls in and where in that gene."""
    c = _db(slug)
    rows = c.execute(
        """SELECT f.seqid, f.start, f."end", f.strand, s.num FROM feat f JOIN seq s ON s.name = f.seqid
           WHERE f.track = ? AND f.fid = ? ORDER BY s.num, f.start LIMIT 1000""",
        (track, fid),
    ).fetchall()
    if not rows:
        raise HTTPException(404, f"no {fid} in track {track}")
    gene_tracks = tuple(
        r["name"] for r in c.execute("SELECT name FROM tracks WHERE kind = 'genes'")
    )
    out = []
    for r in rows:
        gene = where = defline = None
        if gene_tracks:
            hit = _overlapping(
                c,
                r["num"],
                r["start"],
                r["start"],
                f"AND f.type = 'gene' AND f.track IN ({','.join('?' * len(gene_tracks))}) LIMIT 1",
                gene_tracks,
            )
            if hit:
                gene = hit[0]["fid"]
                where = _classify(_first_transcript_parts(c, gene)[1], r["start"])
                d = c.execute(
                    "SELECT defline FROM gene_info WHERE gene = ?", (gene,)
                ).fetchone()
                defline = d["defline"] if d else None
        out.append(
            {
                "seqid": r["seqid"],
                "start": r["start"],
                "end": r["end"],
                "strand": r["strand"],
                "gene": gene,
                "where": where,
                "defline": defline,
            }
        )
    return {"track": track, "fid": fid, "occurrences": out}


@app.get("/api/g/{slug}/fasta", response_class=PlainTextResponse)
def fasta(slug: str, seq: str, start: int, end: int, strand: str = "+", fid: str = ""):
    if end - start + 1 > MAX_FASTA:
        raise HTTPException(400, f"at most {MAX_FASTA:,} bp per request")
    c = _db(slug)
    bases = fetch_bases(slug, c, seq, start, end)
    if strand == "-":
        bases = revcomp(bases)
    head = f">{fid + ' ' if fid else ''}{seq}:{start}-{end}({strand})"
    return (
        head
        + "\n"
        + "\n".join(bases[i : i + 60] for i in range(0, len(bases), 60))
        + "\n"
    )


# ── getting genomes in ────────────────────────────────────────────────────────


class NewUpload(BaseModel):
    filename: str
    size: int


class Source(BaseModel):
    upload: str | None = None  # an upload id
    path: str | None = None  # relative to GENOMES_IMPORT_DIR (admins)


class TrackIn(BaseModel):
    name: str
    kind: str = "features"
    source: Source


class GenomeIn(BaseModel):
    name: str
    description: str = ""
    fasta: Source
    tracks: list[TrackIn] = []
    deflines: Source | None = None


class TracksIn(BaseModel):
    tracks: list[TrackIn] = []
    deflines: Source | None = None


@app.post("/api/uploads")
def new_upload(body: NewUpload, request: Request):
    _may_add(request)
    try:
        return imports.create_upload(request.state.user, body.filename, body.size)
    except imports.ImportError_ as e:
        raise HTTPException(400, str(e)) from None


@app.get("/api/uploads/{uid}")
def upload_status(uid: str, request: Request):
    try:
        return imports.upload_status(uid, request.state.user)
    except KeyError:
        raise HTTPException(404, "no such upload") from None


@app.put("/api/uploads/{uid}")
async def upload_chunk(uid: str, offset: int, request: Request):
    data = bytearray()
    async for part in request.stream():
        data += part
        if len(data) > imports.MAX_CHUNK:
            raise HTTPException(413, "chunk too large")
    try:
        return imports.append_chunk(uid, request.state.user, offset, bytes(data))
    except KeyError:
        raise HTTPException(404, "no such upload") from None
    except imports.ImportError_ as e:
        raise HTTPException(409, str(e)) from None


@app.delete("/api/uploads/{uid}")
def cancel_upload(uid: str, request: Request):
    try:
        imports.upload_status(uid, request.state.user)
    except KeyError:
        raise HTTPException(404, "no such upload") from None
    imports.discard_upload(uid)
    return {"deleted": uid}


@app.get("/api/server-files")
def server_files(request: Request):
    if not request.state.admin:
        raise HTTPException(403, "admins only")
    return imports.server_files()


def _resolve(src: Source, request: Request, uploads: list[str]) -> str:
    if src.upload:
        uploads.append(src.upload)
        try:
            return str(imports.upload_file(src.upload, request.state.user))
        except KeyError:
            raise imports.ImportError_("no such upload") from None
    if src.path:
        if not request.state.admin:
            raise imports.ImportError_("only admins can import files from the server")
        return str(imports.server_file(src.path))
    raise imports.ImportError_("each file needs an upload or a server path")


def _tracks(tracks: list[TrackIn], request: Request, uploads: list[str]) -> list[dict]:
    out, seen = [], set()
    for t in tracks:
        if not TRACK_RE.match(t.name) or t.name in seen:
            raise imports.ImportError_(
                f"track name {t.name!r} is missing, invalid or used twice"
            )
        if t.kind not in ("genes", "features"):
            raise imports.ImportError_(
                f"track kind must be genes or features, not {t.kind!r}"
            )
        seen.add(t.name)
        out.append(
            {
                "name": t.name,
                "kind": t.kind,
                "path": _resolve(t.source, request, uploads),
            }
        )
    return out


@app.post("/api/genomes")
def add_genome(body: GenomeIn, request: Request):
    _may_add(request)
    try:
        slug = slugify(body.name)
        if store.genome_path(slug):
            raise imports.ImportError_(f"a genome called {slug!r} exists already")
        uploads: list[str] = []
        job = {
            "kind": "new_genome",
            "name": body.name.strip()[:80],
            "description": body.description[:500],
            "user": request.state.user,
            "fasta": _resolve(body.fasta, request, uploads),
            "tracks": _tracks(body.tracks, request, uploads),
            "deflines": _resolve(body.deflines, request, uploads)
            if body.deflines
            else None,
            "slug": slug,
            "uploads": uploads,
        }
    except (imports.ImportError_, BuildError) as e:
        raise HTTPException(400, str(e)) from None
    return imports.submit(job)


@app.post("/api/g/{slug}/tracks")
def add_tracks(slug: str, body: TracksIn, request: Request):
    _may_add(request)
    _db(slug).close()
    try:
        uploads: list[str] = []
        job = {
            "kind": "add_tracks",
            "slug": slug,
            "user": request.state.user,
            "tracks": _tracks(body.tracks, request, uploads),
            "deflines": _resolve(body.deflines, request, uploads)
            if body.deflines
            else None,
            "uploads": uploads,
        }
    except imports.ImportError_ as e:
        raise HTTPException(400, str(e)) from None
    if not job["tracks"] and not job["deflines"]:
        raise HTTPException(400, "nothing to add")
    return imports.submit(job)


@app.delete("/api/g/{slug}/tracks/{track}")
def drop_track(slug: str, track: str, request: Request):
    c = _db(slug)
    t = c.execute("SELECT added_by FROM tracks WHERE name = ?", (track,)).fetchone()
    if not t:
        raise HTTPException(404, f"no track {track}")
    _may_change(request, t["added_by"], store.meta(slug).get("added_by", ""))
    return imports.submit(
        {"kind": "drop_track", "slug": slug, "track": track, "user": request.state.user}
    )


@app.get("/api/jobs")
def jobs():
    return imports.list_jobs()


@app.get("/api/jobs/{jid}")
def job(jid: str):
    try:
        return imports.get_job(jid)
    except KeyError:
        raise HTTPException(404, "no such job") from None


def main():
    import uvicorn

    uvicorn.run(
        app,
        host=os.environ.get("GENOMES_HOST", "127.0.0.1"),
        port=int(os.environ.get("GENOMES_PORT", 8766)),
        proxy_headers=True,
    )


if __name__ == "__main__":
    main()
