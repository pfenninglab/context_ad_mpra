"""
Utilities to create RNA-matched and DNA-matched barcode count tables.

This module is a script-friendly, path-robust refactor of notebook:
1.1.1.make_RNAmatch_DNAmatch_20231111.ipynb

Logic is preserved:
- change_names renames columns and separates RNA/DNA columns
- HEK293T: auto-match DNA to RNA by replacing trailing 'R' with 'D'
- THP1/HMC3: match via explicit lookup table (RNA<->DNA)
- enhancer_id is appended from the master table (ID -> enhancer_id)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd


def change_names(count_table_enhancer: pd.DataFrame) -> Tuple[pd.DataFrame, list[str], list[str]]:
    """
    Replicates the notebook's change_names() behavior:
    - rename columns based on parsing the original sample naming scheme
    - return (renamed_df, RNA_columns, DNA_columns)
    """
    column_names: list[str] = []
    for parts in count_table_enhancer.columns.str.split("-"):
        # Notebook logic
        if parts[0] in {"26", "106", "109"}:
            parts[0] = ""
        if parts[2] in {"26", "106", "109"}:
            parts[2] = ""
        if parts[-1] in {"26", "106", "109"}:
            parts[-1] = ""
        if parts[2][:2] == "AP":
            parts[2] = "Mouse" + parts[2]

        # last + "_" + middle + "_" + (ZC...) + "_" + (R/D)
        column_names.append(parts[-1] + "_" + parts[2] + "_" + parts[0][:-1] + "_" + parts[0][-1])

    for i in range(len(column_names)):
        column_names[i] = column_names[i].strip("_")

    count_table_enhancer = count_table_enhancer.copy()
    count_table_enhancer.columns = column_names
    count_table_enhancer = count_table_enhancer.sort_index(axis=1)

    RNA_columns: list[str] = []
    DNA_columns: list[str] = []

    for col in count_table_enhancer.columns:
        # Notebook logic: check last char of the first chunk
        # (kept as-is to preserve behavior)
        if col.split("-")[0][-1] == "R":
            RNA_columns.append(col)
        else:
            DNA_columns.append(col)

    return count_table_enhancer, RNA_columns, DNA_columns


def _read_enhancer_id_map(enhancer_id_file_path: Path) -> pd.Series:
    """
    Reads master table and returns Series indexed by barcode table index:
    ID -> enhancer_id
    """
    s = pd.read_csv(enhancer_id_file_path, index_col="ID")["enhancer_id"]
    # keep series name consistent
    s.name = "enhancer_id"
    return s


def _write_outputs(
    RNA_matched: pd.DataFrame,
    DNA_matched: pd.DataFrame,
    output_directory: Path,
    cell_type: str,
) -> Tuple[Path, Path]:
    output_directory.mkdir(parents=True, exist_ok=True)
    rna_path = output_directory / f"{cell_type}_RNA_matched_barcodes.csv"
    dna_path = output_directory / f"{cell_type}_DNA_matched_barcodes.csv"
    RNA_matched.to_csv(rna_path)
    DNA_matched.to_csv(dna_path)
    return rna_path, dna_path


def match_auto_by_suffix(
    dna_rna_file_path: Path,
    enhancer_id_file_path: Path,
    output_directory: Path,
    cell_type: str,
    rename_dict: Optional[Dict[str, str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Notebook cell for HEK293T:
    - optional renaming of some columns
    - change_names()
    - DNA matched to RNA columns by mapping <RNAcol without last char> + "D"
    - append enhancer_id from master table
    - save outputs
    """
    count_table_barcode = pd.read_csv(dna_rna_file_path, index_col=0)
    if rename_dict is not None:
        count_table_barcode = count_table_barcode.rename(rename_dict, axis=1)

    matched_barcodes, RNA_columns, DNA_columns = change_names(count_table_barcode)
    DNA = matched_barcodes[DNA_columns]
    RNA = matched_barcodes[RNA_columns]

    DNA_matched = pd.DataFrame(index=RNA.index)
    for col in RNA.columns:
        DNA_matched[col] = DNA[col[:-1] + "D"]

    RNA_matched = RNA

    enhancer_id = _read_enhancer_id_map(enhancer_id_file_path)

    RNA_matched = pd.concat([RNA_matched, enhancer_id], axis=1)
    DNA_matched = pd.concat([DNA_matched, enhancer_id], axis=1)

    _write_outputs(RNA_matched, DNA_matched, output_directory, cell_type)
    return RNA_matched, DNA_matched


def match_with_lookup_table(
    dna_rna_file_path: Path,
    enhancer_id_file_path: Path,
    output_directory: Path,
    cell_type: str,
    lookup_table_path: Path,
    rename_dict: Optional[Dict[str, str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Notebook cells for THP1 and HMC3:
    - optional renaming of columns
    - change_names()
    - use explicit lookup table with columns ["DNA", "RNA"]
    - append enhancer_id and save
    """
    count_table_barcode = pd.read_csv(dna_rna_file_path, index_col=0)
    if rename_dict is not None:
        count_table_barcode = count_table_barcode.rename(rename_dict, axis=1)

    matched_barcodes, RNA_columns, DNA_columns = change_names(count_table_barcode)
    DNA = matched_barcodes[DNA_columns]
    RNA = matched_barcodes[RNA_columns]

    DNA_RNA_lookup = pd.read_csv(lookup_table_path)
    DNA_matched = DNA[DNA_RNA_lookup["DNA"]]
    DNA_matched.columns = DNA_RNA_lookup["RNA"]

    RNA_matched = RNA[DNA_RNA_lookup["RNA"]]

    enhancer_id = _read_enhancer_id_map(enhancer_id_file_path)

    RNA_matched = pd.concat([RNA_matched, enhancer_id], axis=1)
    DNA_matched = pd.concat([DNA_matched, enhancer_id], axis=1)

    _write_outputs(RNA_matched, DNA_matched, output_directory, cell_type)
    return RNA_matched, DNA_matched
