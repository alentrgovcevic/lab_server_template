# Demo data sources

All files are public NCBI records, fetched 2026-10-08 with E-utilities (`efetch`). Map files are
renamed into the example ID scheme (`PMC-…`, `EC…`) so the bank's indexer has IDs to work with.
The numbers are made up.

| File in `bank/` | Accession | Record |
|---|---|---|
| `pmc_gb/PMC-03204.gbk`, `clc_gb_pairs/PMC-00147.gbk` | M77789 | Cloning vector pUC19 |
| `pmc_gb/PMC-03207 pBR322 CTRL.gbk` (+ an identical `.1` copy, to demo duplicate detection) | J01749 | Cloning vector pBR322 |
| `pmc_gb/PMC-03405 pACYC184 in PMC-02483.gbk` | X06403 | Cloning vector pACYC184 |
| `pmc_gb/PMC-03460-pEGFP-N1.gb` | U55762 | Cloning vector pEGFP-N1 |
| `pmc_gb/PMC-05000_pCAMBIA1300.gb` | AF234296 | Binary vector pCAMBIA-1300 |
| `pmc_gb/PMC-07001_pCAMBIA1301_binary_vector.gbk` | AF234297 | Binary vector pCAMBIA-1301 |
| `ec_gb/EC25111 pBluescript II SK(+).gb` | X52328 | pBluescript II SK(+) |
| `fasta/EC31492_EGFP_CDS.fa` | U55762 | EGFP CDS, extracted |
| `fasta/PMC-03254_pUC19_MCS_region.fa` | M77789 | pUC19 bases 141–600, extracted |
| `clc/*.clc`, `clc_gb_pairs/PMC-00147.clc` | — | Placeholders with the CLC magic bytes, to demo how unreadable CLC binaries are flagged |

`oligos.tab` holds universal sequencing primers (M13, T7, SP6, T3, CMV, EGFP, pBR322 ori, AmpR)
plus two made-up examples: a FAM-labelled primer and a tailed cloning primer.

`genome/`: the *Lotus japonicus* chloroplast genome, RefSeq **NC_002694.1**. The FASTA and the
GFF3 come from NCBI. `seed.json` tells the genome browser to build it on first start, when it has
no genomes yet; gene descriptions come from the GFF3's `product=` attributes.
