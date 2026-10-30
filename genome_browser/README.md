# genome_browser — any genome, uploaded by the lab

A genome browser for whatever the lab works with. An admin adds a genome (FASTA plus gene models
and any extra tracks) by upload or from a folder on the server, and every signed-in lab member can
browse it once it is built. Set `GENOMES_UPLOADS=everyone` to let every member add genomes. It is
FastAPI over one SQLite file per genome. The page is vanilla JavaScript and SVG, with no JBrowse,
no htslib and no build step.

- **Add a genome** (＋ Add, admins): a FASTA (plain or `.gz`, any line length), optionally gene models
  (GFF3), extra tracks (GFF3 or BED: insertions, repeats, peaks, markers…) and a descriptions TSV.
  - **Uploads** go up in 16 MB chunks. That keeps them under the proxy's request limit, and a
    broken connection resumes. The server refuses an upload when the disk lacks room for it.
  - **Builds** run one at a time in a background process. The page shows the log.
  - **Large genomes:** chromosomes longer than 2 Gb (e.g. faba bean) work. Positions are stored as
    64-bit integers, and the spatial index is only a prefilter.
- **Admins** can also import from a folder already on the server (`GENOMES_IMPORT_DIR`, e.g. the
  lab share mounted read-only) without uploading anything.
- **Browse:** whole-chromosome density, re-binned to at most 2,000 bars. Up to 2 Mb shows the
  features; at 400 kb or less, transcripts with exons, CDS and UTRs; under 150 bp, both strands.
  Click a track name to hide it.
- **Search:** gene ID, gene name, any feature ID (also without its prefix), or words from the
  descriptions. Descriptions come from GFF3 `description`/`Note`/`product`, or the optional TSV.
- **Gene card:** location, transcripts, FASTA (gene or ±1 kb), and every feature from the other
  tracks inside the gene, classed as CDS, UTR exon or intron (from the first transcript). When the
  description carries a PANTHER subfamily (`PTHR48010:SF104`, as in Phytozome deflines), the card
  also lists the gene's copies: every gene of that subfamily in the genome, each a click away.
- **Transcript panel** (on the gene card, one transcript at a time): the spliced mRNA (UTRs in lower
  case), the CDS and the protein, each as a FASTA download, and an exon/intron table (5′ UTR, CDS,
  3′ UTR bases and phase per exon; splice dinucleotides per intron). Derived on request from the
  exon/CDS features and the genome (`GET /api/g/<genome>/transcript/<id>`); nothing new is stored.
  Phase is recomputed from the CDS lengths. Warnings for no ATG, no stop, in-frame stops.
- **Repeated IDs** (e.g. all insertions of one mutant line) get a card listing every location, with
  the gene and gene part each one hits.
- **Who may do what:** every signed-in member sees every genome. Admins add genomes and tracks
  (`GENOMES_UPLOADS=admins`, the default), so the disk holds one curated copy of each. The person
  who added a genome or track, or an admin, can remove it.
- **Deep links:** `#<genome>/<seq>:<start>-<end>` or `#<genome>/gene=<id>`.

## Sign-in

`GENOMES_AUTH=header` (in Docker): Caddy asks the plasmid bank's `/auth/check` before every request
(`forward_auth`) and passes the user and role on as `X-Auth-User` / `X-Auth-Role`. It drops any
such header a browser sends. So the same account works in both apps, and signing out of the bank
signs you out here too. `GENOMES_AUTH=none` is for a laptop: everyone is `local`, an admin.

## Data layout

`GENOMES_DIR/<slug>/` holds `genome.fa` (normalised: 60 bases per line), `genome.sqlite` and
`meta.json`. The schema is in the docstring of `genomebrowser/build.py`. Working space lives in
`.incoming/` (uploads), `.jobs/` (build queue and logs) and `.building/`. A genome appears only when
its build has finished. Adding or removing a track builds a copy of the database and swaps it in.

From the command line (e.g. a big genome already on the server):

```bash
python -m genomebrowser.build --genomes-dir DIR --name "Barley Morex V3" --fasta genome.fa.gz \
  --genes genes.gff3.gz --track repeats=repeats.gff3 --track peaks=peaks.bed [--deflines d.tsv]
```

**The worked example, *L. japonicus* Gifu v1.3:** the genome (assembly LjG1.1) and gene models (with
a LORE1 track) are distributed through [Phytozome](https://phytozome-next.jgi.doe.gov/) under JGI's
data-use terms. Download them there and upload them; they are not redistributed here. The demo
(`deploy/quickstart.sh up`) builds the *L. japonicus* chloroplast from `examples/genome/` on first
start.

## Run without Docker

```bash
../.venv/bin/python -m genomebrowser.server     # http://127.0.0.1:8766, GENOMES_AUTH=none
../.venv/bin/python -m pytest -q                # synthetic genomes, no data needed
```
