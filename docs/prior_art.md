# Prior art: other genome browsers (surveyed 2026-10-08)

Why build one instead of deploying an existing browser? None of the open-source projects below
was confirmed to let a lab member add a genome (FASTA + GFF3) from a web page and have it built on
the server. The common pattern is a CLI or config file plus pre-built indexes, and that is the gap
this browser fills. Licences matter before porting any code: several of these are GPL.

Collected by a web-research pass. The licences come from the GitHub API, not from reading each
LICENSE file, so re-check one before reusing its code.

| Project | Licence | How data gets in | What we took from it |
|---|---|---|---|
| [JBrowse 2](https://github.com/GMOD/jbrowse-components) | Apache-2.0 | `jbrowse add-assembly` / `add-track` CLI; "Open track" for a URL or local file | Whole-genome views can run on sequence lengths alone ([docs](https://jbrowse.org/jb2/docs/config_guides/assemblies/)), so our density view needs only `seq` + `density` |
| [igv.js](https://github.com/igvteam/igv.js) | MIT | URLs to FASTA + `.fai` (+ `.gzi` for bgzip) | A large FASTA needs an index ([docs](https://igv.org/doc/igvjs/Reference-Genome/)). We normalise to fixed-width lines and keep our own `.fai` table |
| [genome-canvas](https://github.com/QuertermousLab/genome-canvas) | MIT | Configured server folders, no upload | The "import from a server folder" idea (`GENOMES_IMPORT_DIR`), and a track picker |
| [GenomeHubs](https://github.com/genomehubs/genomehubs) | MIT | A config file of paths or URLs, run in a container | One build step that takes FASTA + GFF and produces everything the site needs |
| [sequenceserver](https://github.com/wurmlab/sequenceserver) | AGPL-3.0 | FASTA files in a `db` folder (open source) | BLAST UX for a later step: detect the FASTA type and offer to build the database |
| [Lotus Base](https://github.com/lotusbase/lotus.au.dk) | MIT | — (PHP portal, last push 2023) | The original inspiration: the gene card and the LORE1 line view, generalised here into "features in this gene" and "every location of this ID" |
| [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks) | GPL-3.0 | Config file → static PDF/PNG | Not interactive; not used |

**Checked because of big plant genomes:**
- **htslib indexes:** TBI and BAI indexes stop at 512 Mb per chromosome (2^29), so CSI is needed
  above that ([tabix](https://www.htslib.org/doc/tabix.html)). Check your assembly's chromosome
  lengths against this before choosing tools for big plant genomes. We use no htslib index.
- **SQLite:** `rtree_i32` stores 32-bit integers, which top out at 2.1 Gb. We use the float `rtree`
  as a prefilter, then compare exactly on the 64-bit `feat.start/end`. A test covers positions
  around 2.5 Gb.
- **Uploads of many GB:** the [tus protocol](https://tus.io/protocols/resumable-upload) defines
  resumable uploads (an offset, then append). Ours is a small in-house version of the same idea:
  16 MB chunks, each accepted only at the current offset.
