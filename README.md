# AD_MPRA

Massively parallel reporter assay (MPRA) analysis pipeline and downstream analysis code for
a study of Alzheimer's-disease-associated regulatory variants across multiple cell types and
tissues (HEK293T, HMC3 microglia, THP1 monocyte/macrophage under stimulation, iPSC-derived
neurons, brain, and gut).

This repository accompanies the manuscript (citation below, to be updated on acceptance).

## Repository layout

```
src/mpra/           Core Python package: read mapping, count-table construction,
                     RNA/DNA barcode matching, annotation expansion, MPRAnalyze R helpers
scripts/             Pipeline scripts (step01-step09) from raw reads to allele-level count
                     tables, plus R scripts that run MPRAnalyze differential/MAD testing
notebooks/           Downstream analysis notebooks, organized by stage (see notebooks/README.md)
annotation_barcodes/ Per-cell-type barcode-to-element annotation tables
annotation_enhancer/ Enhancer annotation tables
indexing/            REF/ALT lookup tables, RNA/DNA barcode match tables, negative controls
```

See [`DATA.md`](DATA.md) for exactly what is and isn't tracked in this repository, and why
(raw fastq and large intermediate result tables are excluded — see that file for where they
will be deposited).

## Environment setup

```bash
conda env create -f environment.yml
conda activate ad-mpra
pip install -e .
```

This installs the Python dependencies (numpy, pandas, scipy, matplotlib, Biopython, Jupyter)
plus R with the Bioconductor `MPRAnalyze` and `BiocParallel` packages used by the differential
activity and MAD-score analyses.

## Pipeline run order

1. **Read mapping & count tables** (`scripts/step01_map_reads_make_table*.py`,
   `scripts/step01_map_neuron_reads_make_table_streaming.py`) — map raw fastq reads to
   barcodes/enhancers and build raw count tables. Uses `src/mpra/map_reads.py` and
   `src/mpra/make_count_table.py`.
2. **Tissue separation & barcode matching** (`step02_separate_tissues_20231205.py`,
   `step03_make_rna_dna_matched_barcodes_with_neuron_v4.py`) — split by tissue/cell type and
   match RNA counts to DNA counts per barcode (`src/mpra/matching.py`).
3. **Count table conversion & pseudobarcodes** (`step04_convert_count_table.py`,
   `step04_make_neuron_pseudobarcodes.py`, `step05_convert_many_matched_barcodes_to_noSV40.py`).
4. **Annotation & allele-level counts** (`step06_make_annotations_split_counts.py`,
   `step07_make_alleleonly_counts.py`, `step08_convert_neuron_reshaped_to_altref.py`,
   `add_controls_to_altref.py`) — uses `src/mpra/expand_annotation.py` and
   `src/mpra/reshape_noSV40.py`.
5. **MAD-score QC** (`step09_analyze_mad_pseudobarcode_R1R2_nguyen_fdr.py`).
6. **MPRAnalyze differential/MAD testing** (R): `run_mpranalyze_mad_one.R` /
   `batch_run_mpranalyze_mad.sh` (batched via `mad_jobs.csv`) and
   `mpranalyze_comparative_analysis_alleleonly.R`, built on `src/mpra/mpranalyze_utils.R` and
   `src/mpra/mpranalyze_mad_utils.R`.
7. **Downstream analysis & figures**: notebooks in `notebooks/`, run in the stage order
   described in `notebooks/README.md` (annotation → MPRAnalyze → plotting → ML/motif
   benchmarking → risk allele/GWAS).

## Data availability

Raw sequencing data and large intermediate result tables are not stored in this repository —
see [`DATA.md`](DATA.md) for what's included, what's excluded, and where excluded data will be
deposited (SRA/GEO for raw reads; Zenodo for large intermediate tables — accessions to be
added).

## Citation

If you use this code, please cite:

> [Manuscript citation to be added upon publication.]

## License

MIT — see [`LICENSE`](LICENSE).
