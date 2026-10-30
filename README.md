# lab_server_template

> 🚧 **Work in progress.** This started as the in-house server of a plant molecular biology lab and is
> being turned into a template that other labs can clone. It runs and its tests pass, but the
> lab-specific parts (the ID scheme, the folder layout of the share) are examples you will want
> to adapt. Issues and pull requests are welcome.

A self-hosted web server for a molecular biology lab, in one `docker compose` stack:

- **Plasmid map bank + bench tools.** Index the lab's folder of plasmid maps and search, view and
  work with them in the browser.
- **Genome browser for any genome.** Admins add a genome (FASTA) with gene models and any extra
  tracks (GFF3/BED), by upload or from a server folder, and the whole lab can browse it. Built for big
  plant genomes: chromosomes past 2 Gb work in the tests, but a real barley or faba bean build has
  not been timed yet.

Both sit behind one [Caddy](https://caddyserver.com/) on one hostname with HTTPS and share one
sign-in: the bank has local accounts (single sign-on through a header-setting login proxy is
optional), and Caddy asks the bank before letting anyone into the genome browser.

## Quick start (demo data)

Needs Docker and `python3`.

```bash
git clone <this repo> && cd lab_server_template
deploy/quickstart.sh up
```

The script prints a first `admin` password. Then open:

- **Plasmid bank:** <https://localhost:8443/> (self-signed certificate, so accept the browser warning)
- **Genome browser:** <https://localhost:8443/genomes/> (same sign-in)

The demo bank holds public plasmids from NCBI GenBank: pUC19, pBR322, pACYC184, pBluescript II
SK(+), pEGFP-N1, pCAMBIA1300 and pCAMBIA1301. It also has a list of universal sequencing primers.
The demo genome is the *Lotus japonicus* chloroplast (NC_002694.1, 150 kb, 127 genes), built on
first start. Add your own with **＋ Add** on the genome browser's page. Sources are
in [`examples/SOURCES.md`](examples/SOURCES.md).

`deploy/quickstart.sh user NAME` adds an account, and `deploy/quickstart.sh down` stops the stack.

## What the plasmid bank does

| Area | What |
|---|---|
| Index | Scans a folder tree (a network share, mounted read-only) into SQLite. It reads GenBank, FASTA and raw sequence, sniffs the format from content, flags CLC binaries, and writes a hygiene report (duplicates, ambiguous IDs, unreadable files) |
| Browse | Search by ID, name or path; circular/linear map view with features, ORFs, translation and restriction sites; find a pasted sequence on either strand, across the origin; Ctrl/⌘+C copies a selection |
| Sequence view | A Sequence tab with rows of the whole plasmid (the selection highlighted, one-strand selection, a translation row, recognition sequences and cuts, feature and primer arrows 5′→3′), and a FASTA tab with the text to copy |
| Primers | Which oligos in the lab collection bind a map, and where. Primer design with primer3, circular-aware. Site-directed mutagenesis primers. PCR product from two primers, with a cycling program for common polymerases (Q5, Phusion, SuperFi II, OneTaq, Taq) and the vendors' sources |
| Cloning | Restriction digest, ligation, Golden Gate (with every overhang checked base by base, and a `cloning module_<name>` feature for parts that carry no annotation), Gibson design and a virtual gel |
| Verify | Compare Sanger reads (`.ab1`, FASTA, FASTQ) with a map, several reads stacked under the reference, each read renamable; the result stays while the map is open |
| BLAST / align | BLAST+ against the whole bank ("which plasmids contain this?"), MAFFT alignments |
| Edit / clean | Edit features and bases, reverse-complement, join exons into one CDS (frame checked), your own feature types. Save as a new version, a copy in My workspace, or (after two confirmations) over the file on the share, with the old file kept in `_replaced/`. A feature of any length (one base for a point mutation) straight from a selection. The batch GenBank cleaner handles the mundane annotation cleanup and normalises other programs' GenBank dialects (CLC, Vector NTI, ApE: non-standard feature types, locale dates), and `python -m plasmidbank.cleanup_report` lists what to tidy in the share (empty files, unreadable files, names without an ID, copies) |
| People | My workspace (a folder tree of your own files: drag and drop, rename, duplicate, move; synced both ways with `Workspaces/<you>/` on the share), shared projects with draft maps, a common annotation library plus a private one per user, and an admin page with a Basic/Advanced setting for which tools everyone sees |
| Uploads | New maps go into the bank at once, one or several at a time. An optional add-only sync copies them to the share. A map without a lab ID goes to My workspace instead |
| Guides | A short guided tour per page, a "What's new" card per release (`plasmidbank/app/static/whatsnew.js`) and a `/changelog` page listing every release; light/dark theme on every page |
| Labels | QR and Code 128 PNGs for a lab database or label printer, behind a separate bearer token |
| Operations | A usage log and report (who uses what, where people wait, whether the server needs more), a load test, and a "starting" page plus restart countdown in Caddy. See [`deploy/README.md`](deploy/README.md) |

## Adapting it to your lab

- **Map IDs.** The example scheme is the source lab's: `PMC-01234` for plasmids and `EC12345` for
  parts. It lives in `plasmid_maps/plasmidbank/ids.py` (two regexes) and in
  `plasmid_maps/plasmidbank/uploads.py`, which sets the folder layout new uploads are synced into.
  A configurable ID scheme is on the to-do list.
- **Your maps.** Set `BANK_MOUNT` in `deploy/.env` to your folder, then run
  `docker compose -f deploy/docker-compose.yml --env-file deploy/.env --profile tools run --rm scan`.
- **Your oligos.** `OLIGO_FILE` is a tab-separated export with the columns of `examples/oligos.tab`.
- **A bigger common annotation list (optional).** `plasmidbank.annotation_import` loads a feature
  FASTA into the common list. pLannotate's SnapGene-derived feature set works well, but it is not
  shipped here because no licence is stated for that data. Download `plannotate-databases-v2.tar.gz`
  from the [pLannotate releases](https://github.com/mmcguffi/pLannotate/releases). Extract the FASTA
  with `blastdbcmd -db BLAST_dbs/snapgene -entry all -line_length 1000000 > snapgene.fa`. Export the
  `snapgene` table of `BLAST_dbs/snapgene.db` as `snapgene.csv`. Then run
  `python -m plasmidbank.annotation_import snapgene.fa --meta snapgene.csv --state /state --source "pLannotate (SnapGene features)"`
  inside the app container. That is a dry run; add `--write --by <name>` to import.
- **Your genomes.** Upload them on the genome browser's page, or let admins import from a server
  folder (`GENOMES_IMPORT_DIR`). See [`genome_browser/README.md`](genome_browser/README.md).
- **A real server.** Set `SITE_ADDRESS` to its DNS name, the ports to 80/443, and `TLS_DIRECTIVE`
  to your certificate. Schedule `scan` (and `sync`, if you use uploads) with cron or systemd timers.

## Adding a module

The dashboard shows two greyed-out tiles, **Safety quiz** and **Chemical register**: examples of lab
modules that are not part of this template but plug in the same way the bank's own extras do.

1. **Logic** in its own file under `plasmid_maps/plasmidbank/` (no FastAPI import there; raise
   `ValueError`/`KeyError` and turn them into HTTP errors in the route).
2. **Routes and page:** add them to the router in `plasmid_maps/plasmidbank/app/modules.py`, with the
   page as `plasmidbank/app/static/<module>.html`; build its store in `modules.setup()`, which the
   server calls once with the state folder. Sign-in already covers every route; use `_admin(request)`
   for admin-only ones.
3. **Dashboard tile:** in `plasmidbank/app/static/home.html`, turn the stub `<div class="tile off">`
   into `<a class="tile" href="/your-module">`.
4. **Admin card (optional):** `plasmidbank/app/static/admin.html` marks where module cards go, for
   uploads or settings admins manage.

## Layout

```
plasmid_maps/   plasmid bank: FastAPI app (plasmidbank/app), tools, scanner, tests
genome_browser/ genome browser: FastAPI + SQLite per genome, chunked uploads, build queue, tests
deploy/         docker-compose.yml, Caddyfile, .env.example, quickstart.sh
examples/       public demo data: bank/ (plasmids), oligos.tab, genome/ (chloroplast)
```

## Development

```bash
python3 -m venv .venv
.venv/bin/pip install -r plasmid_maps/requirements.txt -r genome_browser/requirements.txt pytest
(cd plasmid_maps && ../.venv/bin/python -m pytest -q)
(cd genome_browser && ../.venv/bin/python -m pytest -q)
(cd plasmid_maps && ../.venv/bin/python -m plasmidbank.app.server)    # http://127.0.0.1:8765, examples/bank
```

BLAST+ and MAFFT are optional outside Docker. Without them, BLAST against the bank and multiple
alignment are unavailable. Clone verification is pure Python and works either way.

## Built with

Python, FastAPI, Biopython, primer3-py, SQLite, vanilla JavaScript and SVG, Caddy and Docker. The
genome browser is modelled on Lotus Base, the former *L. japonicus* community portal
([lotusbase/lotus.au.dk](https://github.com/lotusbase/lotus.au.dk)). It is an independent implementation.

## License

[MIT](LICENSE). The demo data in `examples/` comes from NCBI GenBank/RefSeq (see
`examples/SOURCES.md`).
