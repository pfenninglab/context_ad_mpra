#!/usr/bin/env python3
"""Add control enhancers from a barcode-level table into an altref wide table.

Inputs
------
1) altref_wide_csv: a wide table like *_reshaped_altref*.csv
   - index column (first column) is enhancer_id-like identifiers
   - columns look like: <SAMPLE>_ALT_Barcode_1 ... <SAMPLE>_ALT_Barcode_N

2) long_barcode_csv: a barcode-level table like *_matched_barcodes.csv
   - must contain columns: ID, enhancer_id, plus one column per sample (e.g., Naive_HEK293T_ZC65_R)
   - each row corresponds to one barcode for a given enhancer_id

Controls
--------
Controls are detected as rows where either:
  - enhancer_id contains the substring 'Control', or
  - ID contains the substring 'Control'

For each control enhancer_id, we create ONE row in the wide table.
Barcode ordering is taken from the sorted order of the 'ID' column (stable and reproducible).

Notes
-----
- This script DOES NOT fabricate barcodes beyond the observed ones, except padding with zeros
  if a control has fewer barcodes than the wide table expects (rare).
- It will ERROR if sample names do not match between the two files.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import Dict, List, Tuple

import pandas as pd


ALT_PAT = re.compile(r"^(?P<sample>.+)_ALT_Barcode_(?P<idx>\d+)$")


def parse_alt_columns(cols: List[str]) -> Tuple[List[str], Dict[str, int]]:
    """Return ordered alt columns and max barcode index per sample."""
    sample_to_max: Dict[str, int] = {}
    alt_cols: List[Tuple[str, int, str]] = []  # (sample, idx, colname)

    for c in cols:
        m = ALT_PAT.match(c)
        if not m:
            continue
        s = m.group("sample")
        i = int(m.group("idx"))
        alt_cols.append((s, i, c))
        sample_to_max[s] = max(sample_to_max.get(s, 0), i)

    if not alt_cols:
        raise ValueError("No columns matching '<SAMPLE>_ALT_Barcode_<N>' were found in the wide table.")

    # Sort columns by sample then idx to get canonical order
    alt_cols_sorted = [t[2] for t in sorted(alt_cols, key=lambda x: (x[0], x[1]))]
    return alt_cols_sorted, sample_to_max


def infer_out_path(wide_path: str) -> str:
    base, ext = os.path.splitext(wide_path)
    return f"{base}_withcontrols{ext or '.csv'}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-w", "--wide", required=True, help="Wide altref CSV (no controls yet).")
    ap.add_argument("-l", "--long", required=True, help="Barcode-level CSV with controls.")
    ap.add_argument("-o", "--output", default=None, help="Output CSV path.")
    ap.add_argument("--id-col", default="ID", help="Column name for barcode ID in long table.")
    ap.add_argument("--enhancer-col", default="enhancer_id", help="Column name for enhancer ID in long table.")
    ap.add_argument(
        "--control-token",
        default="Control",
        help="Substring used to detect controls in ID/enhancer_id (default: 'Control').",
    )
    args = ap.parse_args()

    wide_path = args.wide
    long_path = args.long
    out_path = args.output or infer_out_path(wide_path)

    # Read wide table (index in first column)
    wide = pd.read_csv(wide_path, index_col=0)
    wide.index = wide.index.astype(str)

    _, sample_to_max = parse_alt_columns(list(wide.columns))
    samples_wide = sorted(sample_to_max.keys())

    # Read long table
    long_df = pd.read_csv(long_path)
    for col in (args.id_col, args.enhancer_col):
        if col not in long_df.columns:
            raise ValueError(f"'{col}' not found in long table. Columns: {list(long_df.columns)}")

    # Validate sample names match EXACTLY
    sample_cols_long = [c for c in long_df.columns if c not in (args.id_col, args.enhancer_col)]
    missing_in_long = [s for s in samples_wide if s not in sample_cols_long]
    extra_in_long = [c for c in sample_cols_long if c not in samples_wide]
    if missing_in_long:
        raise ValueError("Sample(s) in wide but missing in long: " + ", ".join(missing_in_long))
    if extra_in_long:
        raise ValueError("Sample(s) in long but not in wide: " + ", ".join(extra_in_long))

    # Identify control barcode rows
    token = args.control_token
    id_str = long_df[args.id_col].astype(str)
    enh_str = long_df[args.enhancer_col].astype(str)
    ctrl_mask = id_str.str.contains(token, na=False) | enh_str.str.contains(token, na=False)
    ctrl = long_df.loc[ctrl_mask].copy()
    if ctrl.empty:
        raise ValueError(f"No control rows detected using token: {token!r}")

    # Build ONE wide row per control enhancer_id
    control_rows = []
    control_index = []

    for enh_id, sub in ctrl.groupby(args.enhancer_col, sort=False):
        enh_id = str(enh_id)

        # Deterministic barcode numbering
        sub = sub.sort_values(args.id_col, kind="mergesort")
        n_obs = len(sub)

        row = {}
        for sample in samples_wide:
            max_n = sample_to_max[sample]
            vals = sub[sample].tolist()

            # Pad/truncate to match wide schema (rare; controls usually are exactly max_n)
            if n_obs < max_n:
                vals = vals + [0.0] * (max_n - n_obs)
            elif n_obs > max_n:
                vals = vals[:max_n]

            for i in range(1, max_n + 1):
                col = f"{sample}_ALT_Barcode_{i}"
                v = vals[i - 1]
                row[col] = float(v) if pd.notna(v) else 0.0

        control_rows.append(row)
        control_index.append(enh_id)

    controls_wide = pd.DataFrame(control_rows, index=pd.Index(control_index, name=wide.index.name))
    controls_wide = controls_wide.reindex(columns=wide.columns)

    # Refuse to overwrite anything in the wide table
    overlap = wide.index.intersection(controls_wide.index)
    if len(overlap) > 0:
        raise ValueError(f"Would overwrite {len(overlap)} existing enhancer_id(s); example: {overlap[0]}")

    merged = pd.concat([wide, controls_wide], axis=0)
    merged.to_csv(out_path)

    print(f"Wide (input) rows: {len(wide):,}")
    print(f"Controls added (unique enhancers): {len(controls_wide):,}")
    print(f"Merged rows: {len(merged):,}")
    print(f"Output: {out_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
