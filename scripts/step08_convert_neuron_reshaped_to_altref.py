#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_BASE_DIR = Path("/Volumes/T7/ad_mpra_chen")
DEFAULT_COUNTS_DIR = "outputs/read_counts_R1R2"
DEFAULT_LOOKUP = "indexing/ALT_REF_LookUpTable_amended_20231117.csv"


def convert_reshaped_to_altref(
    input_csv: Path,
    output_csv: Path,
    lookup_csv: Path,
    *,
    keep_controls: bool = False,
    controls_only: bool = False,
    drop_sv40: bool = True,
) -> tuple[int, int]:
    """
    Convert a reshaped REF/ALT count matrix to the notebook-style alt/ref table.

    This mirrors 1.3.3.convert_count_table_Enhancers_20240418.ipynb:
      - keep ALT columns directly
      - map REF columns through the ALT->REF lookup table
      - rename REF columns to ALT so both alleles share sample/barcode columns
      - concatenate ALT rows and mapped REF rows
      - by default keep only alt/ref/motifdisrupt enhancer rows
      - with controls_only, keep only Control rows with ALT-named columns
    """
    df = pd.read_csv(input_csv, index_col=0)

    alt_cols = df.columns[df.columns.astype(str).str.contains("ALT")]
    ref_cols = df.columns[df.columns.astype(str).str.contains("REF")]
    if len(alt_cols) == 0:
        raise ValueError(f"{input_csv} must contain ALT columns.")
    if not controls_only and len(ref_cols) == 0:
        raise ValueError(f"{input_csv} must contain REF columns.")

    alt = df.loc[:, alt_cols].copy()
    alt.index = alt.index.astype(str)

    if controls_only:
        controls = alt.loc[
            alt.index.str.contains("Control", case=False, na=False)
        ].copy()
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        controls.to_csv(output_csv)
        return df.shape[0], controls.shape[0]

    lookup = pd.read_csv(lookup_csv, header=None, names=["alt_id", "ref_id"])

    ref = df.loc[:, ref_cols].copy()
    ref = ref.merge(lookup.set_index("alt_id"), left_index=True, right_index=True)
    ref = ref.set_index("ref_id")
    ref = ref.drop_duplicates()
    ref = ref.rename(columns=lambda x: str(x).replace("REF", "ALT"))

    combined = pd.concat([alt, ref], axis=0)
    combined.index = combined.index.astype(str)

    if drop_sv40:
        combined = combined.loc[
            ~combined.index.str.contains("SV40", case=False, na=False)
        ].copy()

    if not keep_controls:
        row_mask = (
            combined.index.str.contains("ref", case=False, na=False)
            | combined.index.str.contains("alt", case=False, na=False)
            | combined.index.str.contains("motifdisrupt", case=False, na=False)
        )
        combined = combined.loc[row_mask].copy()

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_csv)
    return df.shape[0], combined.shape[0]


def default_jobs(counts_dir: Path) -> list[Path]:
    names = [
        "Neuron_RNA_matched_barcodes_reshaped.csv",
        "Neuron_DNA_matched_barcodes_reshaped.csv",
        "NeuronPseudobarcodes_RNA_matched_barcodes_reshaped.csv",
        "NeuronPseudobarcodes_DNA_matched_barcodes_reshaped.csv",
    ]
    return [counts_dir / name for name in names]


def output_name(input_csv: Path, suffix: str) -> Path:
    if input_csv.name.endswith("_matched_barcodes_reshaped_altref_withcontrol_noSV40.csv"):
        return input_csv.with_name(
            input_csv.name.replace(
                "_matched_barcodes_reshaped_altref_withcontrol_noSV40.csv", suffix
            )
        )
    if not input_csv.name.endswith("_matched_barcodes_reshaped.csv"):
        raise ValueError(f"Unexpected input file name: {input_csv.name}")
    return input_csv.with_name(
        input_csv.name.replace("_matched_barcodes_reshaped.csv", suffix)
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert Neuron reshaped count matrices to notebook-style alt/ref tables."
    )
    parser.add_argument("--base-dir", default=str(DEFAULT_BASE_DIR))
    parser.add_argument("--counts-dir", default=DEFAULT_COUNTS_DIR)
    parser.add_argument("--lookup", default=DEFAULT_LOOKUP)
    parser.add_argument(
        "--files",
        nargs="*",
        help="Optional input file names or paths. Defaults to all Neuron*_matched_barcodes_reshaped.csv files.",
    )
    parser.add_argument(
        "--controls-only",
        action="store_true",
        help="Write only Control rows with ALT-named columns. SV40 controls are kept in this mode.",
    )
    parser.add_argument(
        "--keep-controls",
        action="store_true",
        help="Keep Control rows. Default mimics the notebook and keeps only alt/ref/motifdisrupt rows.",
    )
    parser.add_argument(
        "--keep-sv40",
        action="store_true",
        help="Keep SV40 rows. Default drops rows containing SV40.",
    )
    parser.add_argument(
        "--output-suffix",
        default=None,
        help="Output suffix replacing _matched_barcodes_reshaped.csv.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    base_dir = Path(args.base_dir).expanduser().resolve()
    counts_dir = (base_dir / args.counts_dir).resolve()
    lookup_csv = (base_dir / args.lookup).resolve()

    if not lookup_csv.exists():
        raise FileNotFoundError(f"Lookup table not found: {lookup_csv}")
    if not counts_dir.exists():
        raise FileNotFoundError(f"Counts directory not found: {counts_dir}")

    if args.files:
        jobs = []
        for item in args.files:
            path = Path(item)
            if not path.is_absolute():
                path = counts_dir / path
            jobs.append(path.resolve())
    else:
        jobs = default_jobs(counts_dir)

    if not jobs:
        print(f"[WARN] No Neuron reshaped files found under {counts_dir}")
        return

    output_suffix = args.output_suffix
    if output_suffix is None:
        output_suffix = (
            "_matched_barcodes_reshaped_controls.csv"
            if args.controls_only
            else "_matched_barcodes_reshaped_altref.csv"
        )

    for input_csv in jobs:
        if not input_csv.exists():
            raise FileNotFoundError(f"Input file not found: {input_csv}")
        output_csv = output_name(input_csv, output_suffix)

        if args.dry_run:
            print(f"[DRY] {input_csv} -> {output_csv}")
            continue

        input_rows, output_rows = convert_reshaped_to_altref(
            input_csv,
            output_csv,
            lookup_csv,
            keep_controls=args.keep_controls,
            controls_only=args.controls_only,
            drop_sv40=False if args.controls_only else not args.keep_sv40,
        )
        print(
            f"[OK] {input_csv.name}: {input_rows} input rows -> "
            f"{output_rows} alt/ref rows, wrote {output_csv.name}"
        )

    print("[DONE] step08 Neuron alt/ref conversion complete.")


if __name__ == "__main__":
    main()
