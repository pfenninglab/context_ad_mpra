suppressPackageStartupMessages({
  library("MPRAnalyze")
  library("BiocParallel")
})

as_factor_df <- function(df) {
  for (i in colnames(df)) df[[i]] <- as.factor(df[[i]])
  df
}

select_first_k_of_group <- function(n_cols, group_size = 32, keep_k = 15) {
  idx <- seq_len(n_cols)
  idx[((idx - 1) %% group_size) < keep_k]
}

read_controls <- function(neg_path) {
  negative <- read.csv(neg_path, sep = "\t", header = FALSE)
  as.character(negative$V1)
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

# ---- main runner: analyzeQuantification + testEmpirical (MAD) ----
run_mpranalyze_quantification_mad <- function(
  dna_csv,
  rna_csv,
  annot_csv,
  neg_path,
  out_csv,
  cores = 24,
  group_size = 32,
  keep_k = 15,
  drop_ref_in_annot = TRUE,
  # designs passed in as strings, e.g. "~ Barcode_Allele + Test"
  dnaDesign_str = "~ Barcode_Allele + Test",
  rnaDesign_str = "~ 1",
  # test options
  twoSided = FALSE,
  bpparam_log = TRUE
) {
  # parse formulas
  dnaDesign <- as.formula(dnaDesign_str)
  rnaDesign <- as.formula(rnaDesign_str)

  register(MulticoreParam(cores))
  bpparam <- MulticoreParam(cores, log = bpparam_log, stop.on.error = FALSE)

  # read inputs (same as you did: row.names=1)
  df_DNA <- read.csv(dna_csv, header = TRUE, row.names = 1, check.names = FALSE)
  df_RNA <- read.csv(rna_csv, header = TRUE, row.names = 1, check.names = FALSE)
  annot_DNA <- read.csv(annot_csv, header = TRUE, row.names = 1, check.names = FALSE)

  df_DNA <- as.matrix(df_DNA)
  df_RNA <- as.matrix(df_RNA)

  annot_DNA <- as_factor_df(annot_DNA)
  control <- read_controls(neg_path)

  # your original: drop REF rows in annotation
  if (drop_ref_in_annot) {
    annot_DNA <- annot_DNA[!grepl("REF", rownames(annot_DNA)), , drop = FALSE]
  }

  # your original: keep first 15 of every 32 columns by POSITION
  total_columns <- ncol(df_DNA)
  selected_columns <- select_first_k_of_group(total_columns, group_size = group_size, keep_k = keep_k)

  df_DNA <- df_DNA[, selected_columns, drop = FALSE]
  df_RNA <- df_RNA[, selected_columns, drop = FALSE]
  annot_DNA <- annot_DNA[selected_columns, , drop = FALSE]

  # strongly recommended: align / filter / reorder by annotation rownames
  aligned <- align_counts_to_annotation(df_DNA, df_RNA, annot_DNA)
  df_DNA <- aligned$DNA
  df_RNA <- aligned$RNA
  annot_DNA <- aligned$Annotation

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
    dnaDesign = dnaDesign,
    rnaDesign = rnaDesign
  )

  test_results <- testEmpirical(obj = obj, twoSided = twoSided)
  write.csv(test_results, out_csv)

  invisible(test_results)
}
