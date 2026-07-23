# Data availability

This repository intentionally does not track large raw or intermediate data files. This
section explains what is and isn't included, and how to regenerate what's excluded.

## Tracked in this repository

| Path | Contents | Approx. size |
|---|---|---|
| `annotation_barcodes/` | Per-cell-type barcode annotation tables | 11 MB |
| `annotation_enhancer/` | Enhancer annotation tables | <1 MB |
| `indexing/` | REF/ALT lookup tables, RNA/DNA barcode match tables, negative control list | 26 MB |
| `public_data/` | Public comparison datasets: Bond et al. MPRA results, GTEx/eQTL fine-mapping (`syn2580853`, `syn30863700`), and Cooper et al. *Science* supplementary data (`cooper_et_al/`) | 102 MB |

`public_data/syn2580853_finemapped_sign_eQTLOnly.csv` is ~91 MB — under GitHub's 100 MB
hard limit but close to its 50 MB soft-warning threshold; if this repo grows further,
consider moving it to Git LFS or an external mirror.

## Excluded from this repository

| Excluded data | Why | How to regenerate / where it will live |
|---|---|---|
| `mpra_fastq/` — raw paired-end sequencing reads (~35 GB) | Too large for git; raw sequencing data belongs in a dedicated repository | Will be deposited to SRA/GEO — accession number **TBD** |
| `outputs/` — result tables (allele effects, MPRAnalyze MAD/differential results, QC summaries, matched read counts) | Regenerable output, not source data; kept out of git to keep the repo lean | Regenerate via `scripts/` (fastq → count tables) followed by the MPRAnalyze R scripts, and/or the relevant notebooks in `notebooks/01_mpranalyze/`. Planned for deposit as a Zenodo record — DOI **TBD** |
| Other large intermediate result tables from the original working analysis (~2.7 GB: `allele_differences_*`, `enhancer_activities_*`, full-resolution `read_counts_R1R2`/`read_counts_R2only`, motif-scanning outputs, etc.) | Same as above | Same as above — planned for the same Zenodo record, DOI **TBD** |
| `for_collaborators/` (per-collaborator data slices) | Bespoke data shared privately with named collaborators, not intended for public release | N/A — not part of this repository |

## Regenerating an excluded intermediate

Each notebook loads specific input tables by path. If a notebook fails to load a file that
isn't listed under "Tracked in this repository" above, the file is a regenerable
intermediate: trace it back through the relevant `scripts/step01-09_*` step or the
corresponding `notebooks/01_mpranalyze/` notebook that produces it from raw counts, and
re-run that step against the raw data once obtained from SRA/GEO.
