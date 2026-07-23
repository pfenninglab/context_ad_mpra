# src/mpra/mpranalyze_mad_utils.R

suppressPackageStartupMessages({
  library("MPRAnalyze")
  library("BiocParallel")
})

# ---------- utils ----------
as_factor_df <- function(df) {
  for (i in colnames(df)) df[[i]] <- as.factor(df[[i]])
  df
}

select_first_k_of_group <- function(n_cols, group_size = 32, keep_k = 15) {
  idx <- seq_len(n_cols)
  idx[((idx - 1) %% group_size) < keep_k]
}

align_counts_to_annotation <- function(df_DNA, df_RNA, annot_DNA) {
  desired <- rownames(annot_DNA)
  if (is.null(desired) || length(desired) == 0) stop("Annotation has empty rownames.")

  keep <- desired[desired %in% colnames(df_DNA) & desired %in% colnames(df_RNA)]
  if (length(keep) == 0) {
    message("[DEBUG] Example DNA col: ", ifelse(ncol(df_DNA) > 0, colnames(df_DNA)[1], "NONE"))
    message("[DEBUG] Example RNA col: ", ifelse(ncol(df_RNA) > 0, colnames(df_RNA)[1], "NONE"))
    message("[DEBUG] Example annot row: ", ifelse(length(desired) > 0, desired[1], "NONE"))
    stop("No overlapping columns between DNA/RNA counts and annotation rownames.")
  }

  list(
    DNA = df_DNA[, keep, drop = FALSE],
    RNA = df_RNA[, keep, drop = FALSE],
    Annotation = annot_DNA[keep, , drop = FALSE]
  )
}

read_controls <- function(neg_path) {
  negative <- read.csv(neg_path, sep = "\t", header = FALSE)
  as.character(negative$V1)
}

# ---------- main runner ----------
run_mpranalyze_mad_quantification <- function(
  dna_csv,
  rna_csv,
  annot_csv,
  neg_path,
  out_csv,
  cores = 24,
  group_size = 32,
  keep_k = 15,
  drop_ref_in_annot = TRUE,
  bpparam_log = TRUE
) {
  register(MulticoreParam(cores))
  bpparam <- MulticoreParam(cores, log = bpparam_log, stop.on.error = FALSE)

  # Read inputs (as you did: row.names=1)
  df_DNA <- read.csv(dna_csv, header = TRUE, row.names = 1, check.names = FALSE)
  df_RNA <- read.csv(rna_csv, header = TRUE, row.names = 1, check.names = FALSE)
  annot_DNA <- read.csv(annot_csv, header = TRUE, row.names = 1, check.names = FALSE)

  df_DNA <- as.matrix(df_DNA)
  df_RNA <- as.matrix(df_RNA)

  # Annotation factors
  annot_DNA <- as_factor_df(annot_DNA)

  # Controls
  control <- read_controls(neg_path)

  # Drop REF rows in annotation (your original)
  if (drop_ref_in_annot) {
    annot_DNA <- annot_DNA[!grepl("REF", rownames(annot_DNA)), , drop = FALSE]
  }

  # Select columns (15 of every 32) by POSITION (your original)
  total_columns <- ncol(df_DNA)
  selected_columns <- select_first_k_of_group(total_columns, group_size = group_size, keep_k = keep_k)

  df_DNA <- df_DNA[, selected_columns, drop = FALSE]
  df_RNA <- df_RNA[, selected_columns, drop = FALSE]
  annot_DNA <- annot_DNA[selected_columns, , drop = FALSE]

  # Robust step: ensure counts columns match annotation rownames (drop extras + reorder)
  aligned <- align_counts_to_annotation(df_DNA, df_RNA, annot_DNA)
  df_DNA <- aligned$DNA
  df_RNA <- aligned$RNA
  annot_DNA <- aligned$Annotation

  # MPRAnalyze object
  obj <- MpraObject(
    dnaCounts = df_DNA,
    rnaCounts = df_RNA,
    colAnnot  = annot_DNA,
    control   = control,
    BPPARAM   = bpparam
  )

  obj <- estimateDepthFactors(obj, lib.factor = c("Test"), which.lib = "both", depth.estimator = "uq")

  obj <- analyzeQuantification(
    obj = obj,
    dnaDesign = ~ Barcode_Allele + Test,
    rnaDesign = ~ 1
  )

  test_results <- testEmpirical(obj = obj, twoSided = FALSE)
  write.csv(test_results, out_csv)

  invisible(test_results)
}
