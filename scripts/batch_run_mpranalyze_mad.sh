#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/media/zihengc/T7/ad_mpra_chen"
CORES="${CORES:-24}"
NEG="indexing/mpra3_negatives.csv"
JOBS_CSV="${BASE_DIR}/scripts/mad_jobs.csv"

mkdir -p "${BASE_DIR}/outputs/mpranalyze_results"

# Use python's csv parser to robustly read quoted CSV fields
python3 - <<'PY' "${JOBS_CSV}" "${BASE_DIR}" "${CORES}" "${NEG}"
import csv, sys, subprocess, os

jobs_csv, base_dir, cores, neg = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

with open(jobs_csv, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row["name"]
        dna = row["dna"]
        rna = row["rna"]
        annot = row["annot"]
        out = row["out"]
        dnaDesign = row.get("dnaDesign", "~ Barcode_Allele + Test")
        rnaDesign = row.get("rnaDesign", "~ 1")
        twoSided = row.get("twoSided", "FALSE")

        print(f"[RUN] {name}", flush=True)

        cmd = [
            "Rscript", os.path.join(base_dir, "scripts/run_mpranalyze_mad_one.R"),
            f"--base-dir={base_dir}",
            f"--dna={dna}",
            f"--rna={rna}",
            f"--annot={annot}",
            f"--neg={neg}",
            f"--out={out}",
            f"--cores={cores}",
            f"--dnaDesign={dnaDesign}",
            f"--rnaDesign={rnaDesign}",
            f"--twoSided={twoSided}",
        ]
        subprocess.check_call(cmd)

print("[DONE] All MAD jobs finished.", flush=True)
PY
