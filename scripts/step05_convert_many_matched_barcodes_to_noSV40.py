#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

from mpra.reshape_noSV40 import load_drop_enhancers, reshape_matched_barcodes_to_nosv40


def guess_repo_root() -> Path:
    """
    Auto-detect repo root by searching upwards for pyproject.toml.
    Fallback: assume scripts/ is directly under repo root.
    """
    here = Path(__file__).resolve()
    for p in [here.parent, *here.parents]:
        if (p / "pyproject.toml").exists():
            return p
    return here.parents[1]


def default_files() -> list[str]:
    return [
        "HEK293T_DNA_matched_barcodes.csv",
        "HEK293T_RNA_matched_barcodes.csv",
        "HMC3_DNA_matched_barcodes.csv",
        "HMC3_RNA_matched_barcodes.csv",
        "THP1Macrophage_DNA_matched_barcodes.csv",
        "THP1Macrophage_RNA_matched_barcodes.csv",
        "THP1Monocyte_DNA_matched_barcodes.csv",
        "THP1Monocyte_RNA_matched_barcodes.csv",
    ]


def make_out_name(in_name: str) -> str:
    """
    Output name convention:
      <input_stem>_reshaped_altref_withcontrol_noSV40.csv
    """
    stem = in_name[:-4] if in_name.lower().endswith(".csv") else in_name
    return f"{stem}_reshaped_altref_withcontrol_noSV40.csv"


def main():
    parser = argparse.ArgumentParser(
        description="Reshape matched_barcodes CSVs to wide barcode-per-enhancer format, dropping SV40 and low-frequency enhancers."
    )

    # ✅ Defaults as you requested
    parser.add_argument(
        "--base-dir",
        default=None,
        help="Repo root. If omitted, auto-detect via pyproject.toml.",
    )
    parser.add_argument(
        "--io-dir",
        default="outputs/read_counts_R1R2",
        help="Directory containing input matched_barcodes CSVs (relative to base-dir).",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Where to write outputs. Default = same as io-dir.",
    )
    parser.add_argument(
        "--drop-table",
        default="indexing/low_frequency_enhancer_drop_table_20260126.csv",
        help="CSV with index enhancer_id and column Drop (Y to remove). Set empty string to disable.",
    )
    parser.add_argument(
        "--no-drop-sv40",
        action="store_true",
        help="If set, do NOT drop SV40.",
    )
    parser.add_argument(
        "--sv40-pattern",
        default="SV40",
        help="Pattern to identify SV40 enhancer_id.",
    )
    parser.add_argument(
        "--files",
        nargs="*",
        default=None,
        help="Input filenames. If omitted, uses default 8 files (HEK293T/HMC3/THP1Mac/THP1Mono DNA+RNA).",
    )

    args = parser.parse_args()

    repo_root = Path(args.base_dir).resolve() if args.base_dir else guess_repo_root()
    io_dir = (repo_root / args.io_dir).resolve()
    out_dir = (repo_root / (args.out_dir or args.io_dir)).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load drop table (Drop==Y)
    drop_enhancers = None
    if args.drop_table != "":
        drop_path = Path(args.drop_table)
        if not drop_path.is_absolute():
            drop_path = (repo_root / drop_path).resolve()
        drop_enhancers = load_drop_enhancers(drop_path)

    files = args.files if args.files else default_files()

    for fname in files:
        in_path = io_dir / fname
        if not in_path.exists():
            print(f"[SKIP] missing: {in_path}")
            continue

        df = pd.read_csv(in_path, index_col=0)

        reshaped = reshape_matched_barcodes_to_nosv40(
            df,
            drop_enhancers=drop_enhancers,
            drop_sv40=not args.no_drop_sv40,
            sv40_pattern=args.sv40_pattern,
        )

        out_path = out_dir / make_out_name(fname)
        reshaped.to_csv(out_path)
        print(f"[OK] {fname} -> {out_path.name}  shape={reshaped.shape}")


if __name__ == "__main__":
    main()

