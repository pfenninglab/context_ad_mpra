#!/bin/bash

# Directory containing all .bed and .narrowPeak files
PEAKS_DIR="/ocean/projects/bio210062p/zihengc/ad_ml/expand_peaks_500bp_mm10"
# Output directory
OUTPUT_DIR="/ocean/projects/bio210062p/zihengc/ad_ml/expand_peaks_500bp_mm10_enhancer_5kbfromTSS"
# Annotation file
ANNOT_FILE="/ocean/projects/bio210062p/zihengc/genome/gencode/mm10/gencode.vM15.annotation.gff3.gz"
# Distance parameter
DISTANCE=5000

# Create output directory if not exists
#mkdir -p "$OUTPUT_DIR"

# Find all bed and narrowPeak files (including gzipped versions)
INPUT_FILES=$(find $PEAKS_DIR -type f \( -name "*.bed" -o -name "*.bed.gz" -o -name "*.narrowPeak" -o -name "*.narrowPeak.gz" \) | paste -sd ',' -)

# Check if input files are found
if [ -z "$INPUT_FILES" ]; then
    echo "No input files found in the directory."
    exit 1
fi

# Run the script with all input files as a comma-separated list
sbatch /jet/home/zihengc/ocean_zihengc/ad_ml/scripts/filter_enhancers_from_peaks_with_gff.sh -b "$INPUT_FILES" -o "$OUTPUT_DIR" -a "$ANNOT_FILE" -d "$DISTANCE"
