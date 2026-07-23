#!/usr/bin/env Rscript

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

make_abs <- function(base_dir, p) {
  if (is.null(p)) return(NULL)
  if (grepl("^/", p)) return(p)
  file.path(base_dir, p)
}

BASE_DIR <- get_arg("--base-dir", infer_repo_root())
source(file.path(BASE_DIR, "src/mpra/mpranalyze_utils.R"))

dna_csv   <- get_arg("--dna", NULL)
rna_csv   <- get_arg("--rna", NULL)
annot_csv <- get_arg("--annot", NULL)
out_csv   <- get_arg("--out", NULL)

if (is.null(dna_csv) || is.null(rna_csv) || is.null(annot_csv) || is.null(out_csv)) {
  cat(
    "Usage:\n",
    "  Rscript scripts/run_mpranalyze_mad_one.R \\\n",
    "    --dna=PATH --rna=PATH --annot=PATH --out=PATH \\\n",
    "    [--neg=indexing/mpra3_negatives.csv] [--cores=24]\n",
    "    [--dnaDesign='~ Barcode_Allele + Test'] [--rnaDesign='~ 1']\n",
    "    [--twoSided=FALSE] [--group-size=32] [--keep-k=15]\n",
    sep = ""
  )
  quit(status = 2)
}

neg_path <- get_arg("--neg", "indexing/mpra3_negatives.csv")
cores    <- as.integer(get_arg("--cores", "24"))
group_sz <- as.integer(get_arg("--group-size", "32"))
keep_k   <- as.integer(get_arg("--keep-k", "15"))

dnaDesign_str <- get_arg("--dnaDesign", "~ Barcode_Allele + Test")
rnaDesign_str <- get_arg("--rnaDesign", "~ 1")
twoSided_str  <- get_arg("--twoSided", "FALSE")
twoSided      <- tolower(twoSided_str) %in% c("true", "t", "1", "yes", "y")

run_mpranalyze_quantification_mad(
  dna_csv = make_abs(BASE_DIR, dna_csv),
  rna_csv = make_abs(BASE_DIR, rna_csv),
  annot_csv = make_abs(BASE_DIR, annot_csv),
  neg_path = make_abs(BASE_DIR, neg_path),
  out_csv = make_abs(BASE_DIR, out_csv),
  cores = cores,
  group_size = group_sz,
  keep_k = keep_k,
  dnaDesign_str = dnaDesign_str,
  rnaDesign_str = rnaDesign_str,
  twoSided = twoSided
)

message("[DONE] wrote: ", make_abs(BASE_DIR, out_csv))
