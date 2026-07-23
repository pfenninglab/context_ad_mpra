# scripts/step01_separate_tissues_20231205.py

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

from mpra.sample_sets import (
    HEK293T_COLUMNS,
    THP1_COLUMNS,
    HMC3_COLUMNS,
    BRAIN_COLUMNS,  
)

def repo_root() -> Path:

    return Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="outputs/read_counts_R1R2/df_MPRA3_20231110_R1R2.csv",
        help="Combined master count table CSV.",
    )
    parser.add_argument(
        "--outdir",
        default="outputs/read_counts_R1R2",
        help="Output directory for separated tissue tables.",
    )
    parser.add_argument(
        "--index-col",
        default="ID",
        help="Index column name in the input CSV (use 'ID' if that's what the notebook used).",
    )
    args = parser.parse_args()

    ROOT = repo_root()
    input_path = (ROOT / args.input).resolve()
    outdir = (ROOT / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, index_col=args.index_col)

    # Same notebook logic: subset columns and write
    df[HEK293T_COLUMNS].to_csv(outdir / "HEK293T_DNA_RNA.csv")
    df[THP1_COLUMNS].to_csv(outdir / "THP1_DNA_RNA.csv")
    df[HMC3_COLUMNS].to_csv(outdir / "HMC3_DNA_RNA.csv")


    df[BRAIN_COLUMNS].to_csv(outdir / "Brain_DNA_RNA.csv")

if __name__ == "__main__":
    main()
