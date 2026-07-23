# Analysis notebooks

These notebooks cover the downstream MPRA analysis: statistical modeling, annotation, and
figure generation. They are organized by pipeline stage with their own self-contained
numbering (`00`-`04`), independent of the numbering used inside `scripts/` — `scripts/`
covers the earlier raw-processing steps (fastq to count tables) as a separate pipeline; see
the top-level [`README.md`](../README.md) for how the two fit together.

| Folder | Stage |
|---|---|
| `00_annotation/` | SNP, enhancer, and motif annotation |
| `01_mpranalyze/` | MPRAnalyze differential activity and MAD-score statistical modeling |
| `02_plotting/` | Figure generation: heatmaps, volcano/Manhattan plots |
| `03_ml_motif_benchmark/` | ML logFC comparisons, SNP categorization, SpliceAI, motif enrichment |
| `04_risk_allele_gwas/` | GWAS risk-allele annotation, cell-type specificity |

Individual notebook filenames retain their original numbering from the source analysis
(e.g. `2.1...`, `6.1...`) — only the enclosing folder numbers were reindexed after some
notebooks were pruned, so file-level numbers no longer necessarily match their folder.

## Notes on scope and reproducibility

- These notebooks are kept as a record of the analysis, with their original output figures
  embedded — they were not required to be re-mapped one-to-one onto specific paper figures.
- Several notebooks load intermediate data tables that are not included in this repository
  (see [`../DATA.md`](../DATA.md) for what's tracked vs. excluded and why) — they will not
  re-execute end-to-end without first regenerating those tables from raw data via
  `scripts/` and the MPRAnalyze R scripts there.
