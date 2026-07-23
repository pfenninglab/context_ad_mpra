#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def distribute_values_order(df: pd.DataFrame, n_groups: int = 5) -> pd.DataFrame:
    """
    Match the Brain/Gut pseudobarcode notebook logic.

    Within each enhancer_id, rank barcodes by summed RNA counts and split them
    into n_groups pseudo-barcode groups.
    """
    df = df.copy()
    cols_to_sum = [col for col in df.columns if col != "enhancer_id"]
    df["sum"] = df[cols_to_sum].apply(pd.to_numeric, errors="coerce").sum(axis=1)
    df = df.sort_values(by=["enhancer_id", "sum"], ascending=[True, False], kind="mergesort")

    pieces: list[pd.DataFrame] = []
    for enhancer_id, group in df.groupby("enhancer_id", sort=False):
        group = group.copy()
        group["enhancer_id"] = enhancer_id
        group["Group"] = pd.qcut(
            group["sum"].rank(method="first"),
            n_groups,
            labels=list(range(1, n_groups + 1)),
        )
        pieces.append(group)

    return pd.concat(pieces, axis=0).drop(columns=["sum"])


def sum_barcode_groups(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Pseudo-barcode"] = df["enhancer_id"].astype(str) + "_Group_" + df["Group"].astype(str)
    return df.groupby("Pseudo-barcode", as_index=False).sum(numeric_only=True)


def make_pseudobarcode_table(
    count_df: pd.DataFrame,
    group_info: pd.DataFrame,
    sample_columns: list[str],
) -> pd.DataFrame:
    df = count_df[sample_columns].copy()
    df = pd.merge(
        df,
        group_info[["enhancer_id", "Group", "Pseudo-barcode"]],
        left_index=True,
        right_index=True,
        how="inner",
    )

    pseudo = sum_barcode_groups(df).set_index("Pseudo-barcode").sort_index()
    enhancer_map = (
        group_info[["enhancer_id", "Pseudo-barcode"]]
        .drop_duplicates("Pseudo-barcode")
        .set_index("Pseudo-barcode")
    )
    pseudo = pd.merge(pseudo, enhancer_map, left_index=True, right_index=True)
    return pseudo[[*sample_columns, "enhancer_id"]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make Neuron pseudo-barcode RNA/DNA matched tables using the same "
            "RNA-derived grouping logic as Brain_Gut_pseudobarcodes.ipynb."
        )
    )
    parser.add_argument(
        "--base-dir",
        default=Path(__file__).resolve().parents[1],
        type=Path,
        help="Repository root. Default: parent of scripts/.",
    )
    parser.add_argument(
        "--io-dir",
        default="outputs/read_counts_R1R2",
        help="Directory containing Neuron matched barcode inputs and receiving outputs.",
    )
    parser.add_argument("--cell-type", default="Neuron", help="Input file prefix. Default: Neuron.")
    parser.add_argument(
        "--output-prefix",
        default="NeuronPseudobarcodes",
        help="Output file prefix. Default: NeuronPseudobarcodes.",
    )
    parser.add_argument("--n-groups", default=5, type=int, help="Pseudo-barcode groups per enhancer.")
    parser.add_argument(
        "--enhancer-contains",
        default=None,
        help="Optional case-insensitive substring filter on enhancer_id before grouping.",
    )
    parser.add_argument(
        "--enhancer-excludes",
        default=None,
        help="Optional case-insensitive substring to exclude from enhancer_id before grouping.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base_dir = args.base_dir.resolve()
    io_dir = base_dir / args.io_dir

    rna_path = io_dir / f"{args.cell_type}_RNA_matched_barcodes.csv"
    dna_path = io_dir / f"{args.cell_type}_DNA_matched_barcodes.csv"
    rna_out = io_dir / f"{args.output_prefix}_RNA_matched_barcodes.csv"
    dna_out = io_dir / f"{args.output_prefix}_DNA_matched_barcodes.csv"
    group_out = io_dir / f"{args.output_prefix}_group_assignment.csv"

    rna = pd.read_csv(rna_path, index_col=0)
    dna = pd.read_csv(dna_path, index_col=0)

    sample_columns = [col for col in rna.columns if col != "enhancer_id"]
    missing_in_dna = [col for col in sample_columns if col not in dna.columns]
    if missing_in_dna:
        raise KeyError(f"DNA file is missing RNA sample columns: {missing_in_dna}")
    if "enhancer_id" not in rna.columns or "enhancer_id" not in dna.columns:
        raise KeyError("RNA and DNA inputs must both contain an enhancer_id column.")
    if not rna.index.equals(dna.index):
        raise ValueError("RNA and DNA barcode indexes differ; cannot reuse RNA grouping safely.")

    keep_mask = pd.Series(True, index=rna.index)
    if args.enhancer_contains:
        keep_mask &= rna["enhancer_id"].astype(str).str.contains(
            args.enhancer_contains, case=False, na=False
        )
    if args.enhancer_excludes:
        keep_mask &= ~rna["enhancer_id"].astype(str).str.contains(
            args.enhancer_excludes, case=False, na=False
        )
    rna = rna.loc[keep_mask].copy()
    dna = dna.loc[keep_mask].copy()
    if rna.empty:
        raise ValueError("No barcode rows remain after enhancer filters.")

    group_basis = rna[sample_columns + ["enhancer_id"]].copy()
    enhancer_sizes = group_basis.groupby("enhancer_id").size()
    too_small = enhancer_sizes[enhancer_sizes < args.n_groups]
    if not too_small.empty:
        raise ValueError(
            "Some enhancers have fewer barcode rows than --n-groups: "
            f"{too_small.to_dict()}"
        )
    grouped = distribute_values_order(group_basis, n_groups=args.n_groups)
    grouped["Pseudo-barcode"] = (
        grouped["enhancer_id"].astype(str) + "_Group_" + grouped["Group"].astype(str)
    )

    rna_pseudo = make_pseudobarcode_table(rna, grouped, sample_columns)
    dna_pseudo = make_pseudobarcode_table(dna, grouped, sample_columns)

    rna_pseudo.to_csv(rna_out)
    dna_pseudo.to_csv(dna_out)
    grouped[["enhancer_id", "Group", "Pseudo-barcode"]].to_csv(group_out)

    print(f"[INFO] base_dir: {base_dir}")
    print(f"[INFO] RNA input: {rna_path} {rna.shape}")
    print(f"[INFO] DNA input: {dna_path} {dna.shape}")
    print(f"[INFO] sample columns: {sample_columns}")
    print(f"[INFO] unique enhancers: {rna['enhancer_id'].nunique()}")
    print(f"[INFO] RNA pseudo: {rna_out} {rna_pseudo.shape}")
    print(f"[INFO] DNA pseudo: {dna_out} {dna_pseudo.shape}")
    print(f"[INFO] group assignment: {group_out} {grouped[['enhancer_id', 'Group', 'Pseudo-barcode']].shape}")
    print("[DONE] Neuron pseudobarcode tables finished.")


if __name__ == "__main__":
    main()
