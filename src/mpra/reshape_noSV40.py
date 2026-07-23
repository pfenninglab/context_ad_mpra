from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


def load_drop_enhancers(drop_table_csv: str | Path) -> set[str]:
    """
    Load enhancer IDs to drop from a drop-table CSV.

    Expected format:
      - index: enhancer_id
      - column: 'Drop' with 'Y' marking enhancers to remove
    """
    drop_table_csv = Path(drop_table_csv)
    df = pd.read_csv(drop_table_csv, index_col=0)

    if "Drop" not in df.columns:
        raise ValueError(f"Drop table must contain a 'Drop' column: {drop_table_csv}")

    drop_mask = df["Drop"].astype(str).str.upper() == "Y"
    return set(df.loc[drop_mask].index.astype(str))


def reshape_matched_barcodes_to_nosv40(
    df: pd.DataFrame,
    drop_enhancers: set[str] | None = None,
    drop_sv40: bool = True,
    sv40_pattern: str = "SV40",
    fill_value: int = 0,
    out_col_suffix: str = "_ALT_Barcode_",
) -> pd.DataFrame:
    """
    Convert a matched_barcodes table (barcode-level rows) into a reshaped wide table.

    Input df expectations:
      - index: seq_id (barcode-level rows)
      - columns: sample columns + 'enhancer_id'

    Output:
      - index: enhancer_id
      - columns: <sample>_ALT_Barcode_1..max_barcodes
        where max_barcodes is inferred as the maximum group size across enhancers
      - values: counts preserving original row order within each enhancer group;
        padded with fill_value
    """
    if "enhancer_id" not in df.columns:
        raise ValueError("Input CSV must contain a column named 'enhancer_id'.")

    df = df.copy()
    df["enhancer_id"] = df["enhancer_id"].astype(str)

    # Drop SV40 enhancers
    if drop_sv40 and sv40_pattern:
        df = df.loc[~df["enhancer_id"].str.contains(sv40_pattern, case=False, na=False)].copy()

    # Drop low-frequency enhancers
    if drop_enhancers:
        df = df.loc[~df["enhancer_id"].isin(drop_enhancers)].copy()

    sample_cols = [c for c in df.columns if c != "enhancer_id"]
    if not sample_cols:
        raise ValueError("No sample columns found (only 'enhancer_id' present).")

    # Force integer counts
    df[sample_cols] = df[sample_cols].fillna(0)
    df[sample_cols] = df[sample_cols].astype(np.int64)

    # Infer max barcodes per enhancer
    grp_sizes = df.groupby("enhancer_id", sort=True).size()
    if grp_sizes.empty:
        raise ValueError("No rows remain after filtering; check SV40/drop table settings.")
    max_barcodes = int(grp_sizes.max())

    # Build output columns
    out_cols: list[str] = []
    for col in sample_cols:
        out_cols.extend([f"{col}{out_col_suffix}{i}" for i in range(1, max_barcodes + 1)])

    out = pd.DataFrame(
        fill_value,
        index=grp_sizes.index.astype(str),
        columns=out_cols,
        dtype=np.int64,
    )

    # Fill each enhancer; preserve row order within group as in df
    for enh, g in df.groupby("enhancer_id", sort=True):
        enh = str(enh)
        for col_i, col in enumerate(sample_cols):
            vals = g[col].to_numpy()
            start = col_i * max_barcodes
            out.loc[enh, out_cols[start : start + len(vals)]] = vals

    return out

