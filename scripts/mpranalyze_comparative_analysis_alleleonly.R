#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library("MPRAnalyze")
  library("BiocParallel")
})

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

infer_repo_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", cmd, value = TRUE)
  if (length(f) == 1) {
    script_path <- sub("^--file=", "", f)
    return(normalizePath(file.path(dirname(script_path), "..")))
  }
  normalizePath(getwd())
}

BASE_DIR   <- get_arg("--base-dir", infer_repo_root())
COUNTS_DIR <- get_arg("--counts-dir", file.path(BASE_DIR, "outputs/read_counts_R1R2"))
ANNOT_DIR  <- get_arg("--annot-dir",  file.path(BASE_DIR, "annotation_barcodes"))
INDEX_DIR  <- get_arg("--index-dir",  file.path(BASE_DIR, "indexing"))
OUT_DIR    <- get_arg("--out-dir",    file.path(BASE_DIR, "outputs/mpranalyze_results"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(BASE_DIR, "src/mpra/mpranalyze_utils.R"))

register(MulticoreParam(24))
bpparam <- MulticoreParam(24, log = TRUE, stop.on.error = FALSE)


run_one <- function(tag, dna_csv, rna_csv, annot_csv, out_csv, rnaDesign, reducedDesign) {
  results <- process_datasets_alleleonly(
    dna_file = file.path(COUNTS_DIR, dna_csv),
    rna_file = file.path(COUNTS_DIR, rna_csv),
    annot_file = file.path(ANNOT_DIR, annot_csv),

  )

  obj <- MpraObject(
    dnaCounts = results$DNA,
    rnaCounts = results$RNA,
    colAnnot  = results$Annotation,
    BPPARAM   = bpparam
  )

  obj <- estimateDepthFactors(obj, lib.factor = c("Test"), which.lib = "both", depth.estimator = "uq")

  obj <- analyzeComparative(
    obj = obj,
    dnaDesign = ~ Barcode_Allele + Test,
    rnaDesign = rnaDesign,
    reducedDesign = reducedDesign
  )

  res <- testLrt(obj)
  out_path <- file.path(OUT_DIR, out_csv)
  write.csv(res, out_path)
  message(sprintf("[OK] %s -> %s", tag, out_path))
}

# -------------------------
# THP1Macrophage allele-only (your cell2)
# -------------------------
run_one(
  tag = "THP1Macrophage_alleleOnly",
  dna_csv = "THP1Macrophage_DNA_matched_barcodes_reshaped_allele.csv",
  rna_csv = "THP1Macrophage_RNA_matched_barcodes_reshaped_allele.csv",
  annot_csv = "mpra3_annot_THP1Macrophage_barcodes.csv",
  out_csv = "20240813_comparative_THP1Macrophage_alleleOnly.csv",
  rnaDesign = ~ Tissue + Allele_String,
  reducedDesign = ~ Tissue
)

# -------------------------
# THP1 single-condition alleleOnly (your cell3)
# -------------------------
run_one(
  tag = "THP1_Naive_alleleOnly",
  dna_csv = "THP1_Naive_DNA_matched_barcodes_reshaped_allele.csv",
  rna_csv = "THP1_Naive_RNA_matched_barcodes_reshaped_allele.csv",
  annot_csv = "mpra3_annot_THP1_Naive_barcodes.csv",
  out_csv = "20240616_comparative_THP1_Naive_alleleOnly.csv",
  rnaDesign = ~ Allele_String,
  reducedDesign = ~ 1
)

run_one(
  tag = "THP1_IFNB_alleleOnly",
  dna_csv = "THP1_IFNB_DNA_matched_barcodes_reshaped_allele.csv",
  rna_csv = "THP1_IFNB_RNA_matched_barcodes_reshaped_allele.csv",
  annot_csv = "mpra3_annot_THP1_IFNB_barcodes.csv",
  out_csv = "20240616_comparative_THP1_IFNB_alleleOnly.csv",
  rnaDesign = ~ Allele_String,
  reducedDesign = ~ 1
)

run_one(
  tag = "THP1_LPSIFNG_alleleOnly",
  dna_csv = "THP1_LPSIFNG_DNA_matched_barcodes_reshaped_allele.csv",
  rna_csv = "THP1_LPSIFNG_RNA_matched_barcodes_reshaped_allele.csv",
  annot_csv = "mpra3_annot_THP1_LPSIFNG_barcodes.csv",
  out_csv = "20240616_comparative_THP1_LPSIFNG_alleleOnly.csv",
  rnaDesign = ~ Allele_String,
  reducedDesign = ~ 1
)

run_one(
  tag = "THP1_IFNG_alleleOnly",
  dna_csv = "THP1_IFNG_DNA_matched_barcodes_reshaped_allele.csv",
  rna_csv = "THP1_IFNG_RNA_matched_barcodes_reshaped_allele.csv",
  annot_csv = "mpra3_annot_THP1_IFNG_barcodes.csv",
  out_csv = "20240616_comparative_THP1_IFNG_alleleOnly.csv",
  rnaDesign = ~ Allele_String,
  reducedDesign = ~ 1
)

# -------------------------
# THP1 interaction models (your cell3)
# -------------------------
run_one(
  tag = "THP1_IFNBvsNaive_interaction",
  dna_csv = "THP1_IFNBvsNaive_DNA_matched_barcodes_reshaped_allele.csv",
  rna_csv = "THP1_IFNBvsNaive_RNA_matched_barcodes_reshaped_allele.csv",
  annot_csv = "mpra3_annot_THP1_IFNBvsNaive_barcodes.csv",
  out_csv = "20240812_comparative_THP1_IFNBvsNaive_interaction.csv",
  rnaDesign = ~ Tissue + Allele_String + IFNB_interaction,
  reducedDesign = ~ Tissue + Allele_String
)

run_one(
  tag = "THP1_IFNGvsNaive_interaction",
  dna_csv = "THP1_IFNGvsNaive_DNA_matched_barcodes_reshaped_allele.csv",
  rna_csv = "THP1_IFNGvsNaive_RNA_matched_barcodes_reshaped_allele.csv",
  annot_csv = "mpra3_annot_THP1_IFNGvsNaive_barcodes.csv",
  out_csv = "20240812_comparative_THP1_IFNGvsNaive_interaction.csv",
  rnaDesign = ~ Tissue + Allele_String + IFNG_interaction,
  reducedDesign = ~ Tissue + Allele_String
)

run_one(
  tag = "THP1_LPSIFNGvsIFNG_interaction",
  dna_csv = "THP1_LPSIFNGvsIFNG_DNA_matched_barcodes_reshaped_allele.csv",
  rna_csv = "THP1_LPSIFNGvsIFNG_RNA_matched_barcodes_reshaped_allele.csv",
  annot_csv = "mpra3_annot_THP1_LPSIFNGvsIFNG_barcodes.csv",
  out_csv = "20240812_comparative_THP1_LPSIFNGvsIFNG_interaction.csv",
  rnaDesign = ~ Tissue + Allele_String + LPSIFNG_interaction,
  reducedDesign = ~ Tissue + Allele_String
)

run_one(
  tag = "THP1_LPSIFNGvsNaive_interaction",
  dna_csv = "THP1_LPSIFNGvsNaive_DNA_matched_barcodes_reshaped_allele.csv",
  rna_csv = "THP1_LPSIFNGvsNaive_RNA_matched_barcodes_reshaped_allele.csv",
  annot_csv = "mpra3_annot_THP1_LPSIFNGvsNaive_barcodes.csv",
  out_csv = "20240812_comparative_THP1_LPSIFNGvsNaive_interaction.csv",
  rnaDesign = ~ Tissue + Allele_String + LPSIFNG_interaction,
  reducedDesign = ~ Tissue + Allele_String
)

message("[DONE] allele-only MPRAnalyze finished.")

