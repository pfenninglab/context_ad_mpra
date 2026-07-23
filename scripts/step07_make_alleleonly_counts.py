#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

# ✅ Hard-coded default folder (as requested)
DEFAULT_SEARCH_DIR = Path("/media/zihengc/T7/ad_mpra_chen/outputs/read_counts_R1R2")


def out_name(in_path: Path) -> Path:
    """
    Convert:
      xxx_matched_barcodes_reshaped.csv
    to:
      xxx_matched_barcodes_reshaped_allele.csv
    """
    name = in_path.name
    if not name.endswith("_matched_barcodes_reshaped.csv"):
        raise ValueError(f"Unexpected filename: {name}")
    return in_path.with_name(
        name.replace(
            "_matched_barcodes_reshaped.csv",
            "_matched_barcodes_reshaped_allele.csv",
        )
    )


def main():
    parser = argparse.ArgumentParser(
        description="Extract rows with 'alt' in the index from *_matched_barcodes_reshaped.csv files."
    )
    parser.add_argument(
        "--search-dir",
        default=str(DEFAULT_SEARCH_DIR),
        help=f"Directory to search recursively (default: {DEFAULT_SEARCH_DIR}).",
    )
    parser.add_argument(
        "--pattern",
        default="*matched_barcodes_reshaped.csv",
        help="Glob pattern to match input files (searched recursively).",
    )
    parser.add_argument(
        "--index-col",
        default=0,
        type=int,
        help="Which column to use as index when reading CSV (default 0).",
    )
    parser.add_argument(
        "--contains",
        default="alt",
        help="Substring to match in index (case-insensitive). Default 'alt'.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="If set, do not write outputs; only print what would happen.",
    )
    args = parser.parse_args()

    search_dir = Path(args.search_dir).expanduser().resolve()
    if not search_dir.exists():
        raise FileNotFoundError(f"Search directory not found: {search_dir}")

    files = sorted(search_dir.rglob(args.pattern))
    if not files:
        print(f"[WARN] No files matched pattern '{args.pattern}' under {search_dir}")
        return

    needle = str(args.contains).lower()
    wrote = 0
    skipped = 0

    for fp in files:
        # Avoid re-processing already-generated allele outputs
        if fp.name.endswith("_matched_barcodes_reshaped_allele.csv"):
            continue
        # Only process intended inputs
        if not fp.name.endswith("_matched_barcodes_reshaped.csv"):
            continue

        try:
            df = pd.read_csv(fp, index_col=args.index_col)
        except Exception as e:
            print(f"[SKIP] Failed to read {fp}: {e}")
            skipped += 1
            continue

        idx = df.index.astype(str)
        mask = idx.str.lower().str.contains(needle, na=False)
        df_alt = df.loc[mask].copy()

        out_fp = out_name(fp)

        if args.dry_run:
            print(f"[DRY] {fp}: {mask.sum()} / {len(df)} rows -> {out_fp.name}")
            continue

        df_alt.to_csv(out_fp)
        print(f"[OK] {fp}: {mask.sum()} / {len(df)} rows -> {out_fp.name}")
        wrote += 1

    if not args.dry_run:
        print(f"[DONE] Wrote {wrote} file(s). Skipped {skipped} file(s).")


if __name__ == "__main__":
    main()
