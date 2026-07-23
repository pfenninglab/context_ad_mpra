# Analysis notebooks

These notebooks cover the downstream MPRA analysis: statistical modeling, annotation, and
figure generation. They are organized by pipeline stage using the numeric prefix already
present in each notebook's original filename (category `3` was never used upstream and is
skipped intentionally, not a gap in this reorganization).

| Folder | Stage |
|---|---|
| `00_annotation/` | SNP, enhancer, and motif annotation |
| `01_count_processing/` | Count-table construction, tissue splitting, pseudobarcodes, dropout QC |
| `02_mpranalyze/` | MPRAnalyze differential activity and MAD-score statistical modeling |
| `04_differential_annotation/` | Annotating differential results, splitting by contributor |
| `05_plotting/` | Figure generation: heatmaps, volcano/Manhattan plots, Xenium probe selection |
| `06_ml_motif_benchmark/` | ML logFC comparisons, SNP categorization, SpliceAI, motif enrichment |
| `07_atac_overlap/` | ATAC accessibility vs. ML cluster overlap |
| `08_risk_allele_gwas/` | GWAS risk-allele annotation, cell-type specificity, comparisons to published MPRA datasets |
| `misc/` | Notebooks without a pipeline-stage prefix in the original analysis |

## Notes on scope and reproducibility

- These notebooks are kept as a record of the analysis, with their original output figures
  embedded — they were not required to be re-mapped one-to-one onto specific paper figures.
- A few near-duplicate filenames from the original analysis (e.g. two `2.0.2...` MAD-scoring
  notebooks, or `5.2.6.plot_manhattan.ipynb`) each had an explicit `(Copy)` sibling that was an
  older, superseded version (by file modification time) — those copies were dropped, keeping
  only the newer file.
- `6.3.motif_enrichment_snps_15bp_in_progress.ipynb` was left unfinished in the original analysis
  (15bp window); `6.3.motif_enrichment_snps_25bp.ipynb` is the completed 25bp analysis.
- Several notebooks load intermediate data tables that are not included in this repository
  (see [`../DATA.md`](../DATA.md) for what's tracked vs. excluded and why) — they will not
  re-execute end-to-end without first regenerating those tables from raw data via
  `scripts/step01-08_*` and the MPRAnalyze R scripts in `scripts/`.
