#!/bin/bash
set -euo pipefail

# Directory containing all .bed and .narrowPeak files
PEAKS_DIR="/ocean/projects/bio210062p/zihengc/ad_ml/expand_peaks_500bp_hg38"
# Output directory
OUTPUT_DIR="/jet/home/zihengc/ocean_zihengc/ad_ml/expand_peaks_500bp_hg38_enhancer_10kbfromTSS"
# Annotation file
ANNOT_FILE="/ocean/projects/bio210062p/zihengc/genome/gencode/hg38/gencode.v27.annotation.gff3.gz"
# Distance parameter
DISTANCE=10000

# Keep only canonical chromosomes (edit if needed)
# If you don't want chrM, remove it from the regex.
KEEP_CHR_REGEX='^(chr([1-9]|1[0-9]|2[0-2])|chrX|chrY|chrM)$'

# Create output directory if not exists
mkdir -p "$OUTPUT_DIR"

# Temp dir for filtered peak files
FILTERED_DIR="${OUTPUT_DIR}/_filtered_inputs_canonical_chr"
mkdir -p "$FILTERED_DIR"

echo "[INFO] Filtering peaks to canonical chromosomes into: $FILTERED_DIR"
echo "[INFO] KEEP_CHR_REGEX = $KEEP_CHR_REGEX"

# Find all bed and narrowPeak files (including gzipped versions)
mapfile -t RAW_FILES < <(find "$PEAKS_DIR" -type f \( -name "*.bed" -o -name "*.bed.gz" -o -name "*.narrowPeak" -o -name "*.narrowPeak.gz" \) | sort)

if [ "${#RAW_FILES[@]}" -eq 0 ]; then
  echo "No input files found in the directory."
  exit 1
fi

# Filter each file
for f in "${RAW_FILES[@]}"; do
  bn="$(basename "$f")"
  out="${FILTERED_DIR}/${bn}"

  # Skip if already filtered (optional)
  # if [ -s "$out" ]; then
  #   echo "[SKIP] exists: $out"
  #   continue
  # fi

  if [[ "$f" == *.gz ]]; then
    # Stream-decompress -> filter by chr -> recompress
    # Uses gzip -c to write .gz output
    zcat "$f" \
      | awk -v re="$KEEP_CHR_REGEX" 'BEGIN{FS=OFS="\t"} $1 ~ re {print}' \
      | gzip -c > "$out"
  else
    awk -v re="$KEEP_CHR_REGEX" 'BEGIN{FS=OFS="\t"} $1 ~ re {print}' "$f" > "$out"
  fi

  # If filtering produced empty file, warn (still keep it, but you can choose to delete)
  if [ ! -s "$out" ]; then
    echo "[WARN] Filtered file is empty: $out (from $f)"
  fi
done

# Build comma-separated list from filtered files
INPUT_FILES=$(find "$FILTERED_DIR" -type f \( -name "*.bed" -o -name "*.bed.gz" -o -name "*.narrowPeak" -o -name "*.narrowPeak.gz" \) | sort | paste -sd ',' -)

if [ -z "$INPUT_FILES" ]; then
  echo "No filtered input files found in: $FILTERED_DIR"
  exit 1
fi

echo "[INFO] Submitting sbatch with filtered inputs..."
sbatch /jet/home/zihengc/ocean_zihengc/ad_ml/scripts/filter_enhancers_from_peaks_with_gff.sh \
  -b "$INPUT_FILES" \
  -o "$OUTPUT_DIR" \
  -a "$ANNOT_FILE" \
  -d "$DISTANCE"