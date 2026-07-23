#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable


def bootstrap_local_repo_imports() -> None:
    """
    Make the local repo importable when this script is run from scripts/.

    This fixes: ModuleNotFoundError: No module named 'mpra'

    Supports both common layouts:
      - repo/mpra
      - repo/src/mpra

    For your repository layout, the important path is:
      /Volumes/T7/ad_mpra_chen/src

    The function searches upward from both this script and the current working
    directory, then adds any directory that can expose the mpra package.
    """
    here = Path(__file__).resolve()
    cwd = Path.cwd().resolve()

    roots: list[Path] = []
    for start in (here.parent, cwd):
        roots.append(start)
        roots.extend(start.parents)

    candidates: list[Path] = []
    for root in roots:
        # Flat package layout: root/mpra
        candidates.append(root)
        # src package layout: root/src/mpra
        candidates.append(root / "src")

    seen: set[str] = set()
    for candidate in candidates:
        candidate = candidate.resolve()
        candidate_str = str(candidate)
        if candidate_str in seen:
            continue
        seen.add(candidate_str)

        if (candidate / "mpra").is_dir():
            if candidate_str not in sys.path:
                sys.path.insert(0, candidate_str)



bootstrap_local_repo_imports()

import pandas as pd

try:
    from mpra.matching import match_auto_by_suffix, match_with_lookup_table
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "Could not import 'mpra.matching'. Run this script from the repo root or scripts/ "
        "inside /Volumes/T7/ad_mpra_chen. For src layout, run with: "
        "PYTHONPATH=/Volumes/T7/ad_mpra_chen/src python <script> ... "
        "or install the local package with: "
        "python -m pip install -e /Volumes/T7/ad_mpra_chen"
    ) from exc


VALID_CELLS = {"HEK293T", "THP1Macrophage", "THP1Monocyte", "HMC3", "Neuron"}


def guess_repo_root() -> Path:
    """Find repo/data root by walking upward until pyproject.toml is found."""
    here = Path(__file__).resolve()
    for p in [here.parent, *here.parents]:
        if (p / "pyproject.toml").exists():
            return p
    return here.parents[1]


def resolve_path(repo_root: Path, path_like: str | Path) -> Path:
    """Resolve absolute paths as-is; resolve relative paths under repo_root."""
    p = Path(path_like)
    if p.is_absolute():
        return p.resolve()
    return (repo_root / p).resolve()


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")


def make_enhancer_map_from_input(
    dna_rna_file_path: Path,
    out_dir: Path,
    cell_type: str,
) -> Path:
    """
    Create a minimal enhancer map from a DNA/RNA count table that already has
    ID and enhancer_id columns. This is useful for Neuron_DNA_RNA.csv, which
    already contains both columns.
    """
    header = pd.read_csv(dna_rna_file_path, nrows=0).columns.tolist()
    required = {"ID", "enhancer_id"}
    if not required.issubset(header):
        missing = sorted(required - set(header))
        raise FileNotFoundError(
            "Master enhancer map was not found, and the input table cannot be "
            f"used to build one because it is missing columns: {missing}"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    fallback_path = out_dir / f"{cell_type}_enhancer_map_from_input.csv"
    enhancer_map = pd.read_csv(dna_rna_file_path, usecols=["ID", "enhancer_id"])
    enhancer_map = enhancer_map.drop_duplicates(subset=["ID"]).set_index("ID")
    enhancer_map.to_csv(fallback_path)
    print(f"[INFO] Created fallback enhancer map: {fallback_path}")
    return fallback_path


def resolve_enhancer_map(
    dna_rna_file_path: Path,
    master_map: Path,
    out_dir: Path,
    cell_type: str,
    allow_fallback_from_input: bool,
) -> Path:
    """Use master_map if it exists; optionally build one from the input table."""
    if master_map.exists():
        return master_map
    if allow_fallback_from_input:
        return make_enhancer_map_from_input(dna_rna_file_path, out_dir, cell_type)
    raise FileNotFoundError(f"Missing master enhancer map: {master_map}")


def get_r_or_d_suffix(column_name: str) -> tuple[str, str] | None:
    """
    Return (kind, suffix) for simple ZC sample names.

    Example:
      ZC130R-MPRA3-iPSCNeuron-26-Rep1 -> ("R", "MPRA3-iPSCNeuron-26-Rep1")
      ZC126D-MPRA3-iPSCNeuron-26-Rep1 -> ("D", "MPRA3-iPSCNeuron-26-Rep1")
    """
    if "-" not in column_name:
        return None
    first_token, suffix = column_name.split("-", 1)
    match = re.match(r"^ZC\d+([RD])$", first_token)
    if not match:
        return None
    return match.group(1), suffix


def validate_auto_suffix_pairs(dna_rna_file_path: Path, cell_type: str) -> None:
    """
    Lightweight validation for tables that should be matched by shared suffix.
    It checks that each RNA suffix has a corresponding DNA suffix.
    """
    columns = pd.read_csv(dna_rna_file_path, nrows=0).columns.tolist()
    r_by_suffix: dict[str, list[str]] = {}
    d_by_suffix: dict[str, list[str]] = {}

    for col in columns:
        parsed = get_r_or_d_suffix(col)
        if parsed is None:
            continue
        kind, suffix = parsed
        if kind == "R":
            r_by_suffix.setdefault(suffix, []).append(col)
        elif kind == "D":
            d_by_suffix.setdefault(suffix, []).append(col)

    paired_suffixes = sorted(set(r_by_suffix) & set(d_by_suffix))
    if not paired_suffixes:
        raise ValueError(
            f"No RNA/DNA suffix pairs detected in {dna_rna_file_path}. "
            "Use --single-lookup-table if this table requires a lookup table."
        )

    missing_dna = sorted(set(r_by_suffix) - set(d_by_suffix))
    missing_rna = sorted(set(d_by_suffix) - set(r_by_suffix))
    if missing_dna or missing_rna:
        print(f"[WARN] {cell_type}: incomplete RNA/DNA suffix pairs detected.")
        if missing_dna:
            print(f"[WARN] RNA suffixes without DNA: {missing_dna}")
        if missing_rna:
            print(f"[WARN] DNA suffixes without RNA: {missing_rna}")

    print(f"[INFO] {cell_type}: detected {len(paired_suffixes)} RNA/DNA suffix pair(s).")
    for suffix in paired_suffixes:
        print(f"[INFO]   suffix={suffix} | RNA={r_by_suffix[suffix]} | DNA={d_by_suffix[suffix]}")


def build_suffix_match_rename_dict(dna_rna_file_path: Path, cell_type: str) -> dict[str, str]:
    """
    Rename DNA columns so legacy auto matching can pair by RNA column stem.

    mpra.matching.match_auto_by_suffix() expects each RNA column to have a DNA
    column with the same sample stem and only the final R changed to D. Neuron
    uses shared suffixes instead, e.g. ZC130R-...-Rep1 pairs with
    ZC126D-...-Rep1. For those cases, temporarily rename the DNA column to
    ZC130D-...-Rep1 before handing off to the legacy matcher.
    """
    columns = pd.read_csv(dna_rna_file_path, nrows=0).columns.tolist()
    r_by_suffix: dict[str, list[str]] = {}
    d_by_suffix: dict[str, list[str]] = {}

    for col in columns:
        parsed = get_r_or_d_suffix(col)
        if parsed is None:
            continue
        kind, suffix = parsed
        if kind == "R":
            r_by_suffix.setdefault(suffix, []).append(col)
        elif kind == "D":
            d_by_suffix.setdefault(suffix, []).append(col)

    rename_dict: dict[str, str] = {}
    for suffix in sorted(set(r_by_suffix) & set(d_by_suffix)):
        r_cols = r_by_suffix[suffix]
        d_cols = d_by_suffix[suffix]
        if len(r_cols) != 1 or len(d_cols) != 1:
            raise ValueError(
                f"{cell_type}: expected exactly one RNA and one DNA column for suffix "
                f"{suffix!r}; found RNA={r_cols}, DNA={d_cols}"
            )

        r_first_token = r_cols[0].split("-", 1)[0]
        d_first_token = d_cols[0].split("-", 1)[0]
        expected_d_first_token = r_first_token[:-1] + "D"
        if d_first_token != expected_d_first_token:
            rename_dict[d_cols[0]] = expected_d_first_token + "-" + suffix

    if rename_dict:
        print(f"[INFO] {cell_type}: DNA column rename map for suffix matching: {rename_dict}")
    return rename_dict


def get_auto_match_sample_columns(dna_rna_file_path: Path) -> list[str]:
    """
    Return ZC RNA/DNA sample count columns that can be used by match_auto_by_suffix.

    The Neuron_DNA_RNA.csv file includes metadata columns such as ID,
    enhancer_id, barcode_id, barcode_seq, and Barcode_RC_* columns. The current
    mpra.matching.change_names() implementation assumes every non-index column is
    a hyphen-delimited sample name. Passing metadata columns to it causes:

        IndexError: list index out of range

    Therefore, for Neuron we build a clean temporary count matrix whose first
    column is ID and whose remaining columns are only ZC...R / ZC...D sample
    count columns.
    """
    columns = pd.read_csv(dna_rna_file_path, nrows=0).columns.tolist()
    sample_columns = [col for col in columns if get_r_or_d_suffix(col) is not None]
    if not sample_columns:
        raise ValueError(f"No ZC RNA/DNA sample columns found in {dna_rna_file_path}")
    return sample_columns


def make_clean_auto_match_input(
    *,
    dna_rna_file_path: Path,
    out_dir: Path,
    cell_type: str,
    id_column: str = "ID",
) -> Path:
    """
    Create a temporary input CSV compatible with mpra.matching.match_auto_by_suffix.

    The compatible format is:
      first column: ID
      remaining columns: RNA/DNA count sample columns only

    This avoids two common failure modes for richer count tables:
      1. the original first column is not ID, so index_col=0 would be wrong;
      2. metadata columns are passed into mpra.matching.change_names().
    """
    header = pd.read_csv(dna_rna_file_path, nrows=0).columns.tolist()
    if id_column not in header:
        raise KeyError(
            f"Cannot build clean auto-match input for {cell_type}; "
            f"missing required ID column: {id_column}"
        )

    sample_columns = get_auto_match_sample_columns(dna_rna_file_path)
    clean_columns = [id_column, *sample_columns]

    out_dir.mkdir(parents=True, exist_ok=True)
    clean_path = out_dir / f"_{cell_type}_DNA_RNA.clean_for_auto_match.csv"

    df = pd.read_csv(dna_rna_file_path, usecols=clean_columns)
    df = df[clean_columns]

    if df[id_column].isna().any():
        n_missing = int(df[id_column].isna().sum())
        raise ValueError(f"{cell_type}: {id_column} contains {n_missing} missing value(s).")

    duplicated_ids = int(df[id_column].duplicated().sum())
    if duplicated_ids:
        print(f"[WARN] {cell_type}: {duplicated_ids} duplicate {id_column} value(s) detected.")

    df.to_csv(clean_path, index=False)
    print(f"[INFO] {cell_type}: wrote clean auto-match input: {clean_path}")
    print(f"[INFO] {cell_type}: clean input columns: {clean_columns}")
    return clean_path


def run_clean_auto_match(
    *,
    dna_rna_file_path: Path,
    master_map: Path,
    out_dir: Path,
    cell_type: str,
    rename_dict: dict[str, str] | None = None,
    skip_missing: bool = False,
    validate_suffix_pairs: bool = True,
    allow_enhancer_map_fallback: bool = True,
) -> None:
    """
    Safe auto-match path for rich count tables like Neuron_DNA_RNA.csv.

    It first strips metadata columns and reorders the table so ID is the first
    column, then delegates to the existing mpra.matching.match_auto_by_suffix().
    """
    if not dna_rna_file_path.exists():
        message = f"Missing input table for {cell_type}: {dna_rna_file_path}"
        if skip_missing:
            print(f"[SKIP] {message}")
            return
        raise FileNotFoundError(message)

    if validate_suffix_pairs:
        validate_auto_suffix_pairs(dna_rna_file_path, cell_type)

    suffix_rename_dict = build_suffix_match_rename_dict(dna_rna_file_path, cell_type)
    if rename_dict:
        suffix_rename_dict.update(rename_dict)

    clean_input = make_clean_auto_match_input(
        dna_rna_file_path=dna_rna_file_path,
        out_dir=out_dir,
        cell_type=cell_type,
    )

    enhancer_map = resolve_enhancer_map(
        dna_rna_file_path=dna_rna_file_path,
        master_map=master_map,
        out_dir=out_dir,
        cell_type=cell_type,
        allow_fallback_from_input=allow_enhancer_map_fallback,
    )

    print(f"[RUN] clean match_auto_by_suffix: {cell_type}")
    match_auto_by_suffix(
        dna_rna_file_path=clean_input,
        enhancer_id_file_path=enhancer_map,
        output_directory=out_dir,
        cell_type=cell_type,
        rename_dict=suffix_rename_dict or None,
    )


def run_auto_match(
    *,
    dna_rna_file_path: Path,
    master_map: Path,
    out_dir: Path,
    cell_type: str,
    rename_dict: dict[str, str] | None = None,
    skip_missing: bool = False,
    validate_suffix_pairs: bool = False,
    allow_enhancer_map_fallback: bool = True,
) -> None:
    if not dna_rna_file_path.exists():
        message = f"Missing input table for {cell_type}: {dna_rna_file_path}"
        if skip_missing:
            print(f"[SKIP] {message}")
            return
        raise FileNotFoundError(message)

    if validate_suffix_pairs:
        validate_auto_suffix_pairs(dna_rna_file_path, cell_type)

    enhancer_map = resolve_enhancer_map(
        dna_rna_file_path=dna_rna_file_path,
        master_map=master_map,
        out_dir=out_dir,
        cell_type=cell_type,
        allow_fallback_from_input=allow_enhancer_map_fallback,
    )

    print(f"[RUN] match_auto_by_suffix: {cell_type}")
    match_auto_by_suffix(
        dna_rna_file_path=dna_rna_file_path,
        enhancer_id_file_path=enhancer_map,
        output_directory=out_dir,
        cell_type=cell_type,
        rename_dict=rename_dict,
    )


def run_lookup_match(
    *,
    dna_rna_file_path: Path,
    master_map: Path,
    out_dir: Path,
    cell_type: str,
    lookup_table_path: Path,
    rename_dict: dict[str, str] | None = None,
    skip_missing: bool = False,
    allow_enhancer_map_fallback: bool = True,
) -> None:
    missing = []
    if not dna_rna_file_path.exists():
        missing.append(f"input table: {dna_rna_file_path}")
    if not lookup_table_path.exists():
        missing.append(f"lookup table: {lookup_table_path}")
    if missing:
        message = f"Missing file(s) for {cell_type}: " + "; ".join(missing)
        if skip_missing:
            print(f"[SKIP] {message}")
            return
        raise FileNotFoundError(message)

    enhancer_map = resolve_enhancer_map(
        dna_rna_file_path=dna_rna_file_path,
        master_map=master_map,
        out_dir=out_dir,
        cell_type=cell_type,
        allow_fallback_from_input=allow_enhancer_map_fallback,
    )

    print(f"[RUN] match_with_lookup_table: {cell_type}")
    match_with_lookup_table(
        dna_rna_file_path=dna_rna_file_path,
        enhancer_id_file_path=enhancer_map,
        output_directory=out_dir,
        cell_type=cell_type,
        lookup_table_path=lookup_table_path,
        rename_dict=rename_dict,
    )


def patch_hek293t_zc65_inplace(
    out_dir: Path,
    filename: str = "HEK293T_DNA_matched_barcodes.csv",
) -> None:
    """
    Patch HEK293T DNA matched table by adding/replacing Naive_HEK293T_ZC65_R,
    then overwrite the same file in-place. No backup is created.

    Logic:
      zc65 = round((df[selected] / df[selected].sum()).mean(axis=1) * 5,000,000)
    """
    path = out_dir / filename
    require_file(path, "HEK293T DNA matched table for patching")

    df = pd.read_csv(path, index_col=0)

    selected_columns = [
        "Naive_HEK293T_ZC66_R",
        "Naive_HEK293T_ZC71_R",
        "Naive_HEK293T_ZC76_R",
        "Naive_HEK293T_ZC81_R",
    ]
    missing = [c for c in selected_columns if c not in df.columns]
    if missing:
        raise KeyError(f"Cannot patch HEK293T ZC65; missing required columns: {missing}")

    zc65 = round((df[selected_columns] / df[selected_columns].sum()).mean(axis=1) * 5_000_000)
    df["Naive_HEK293T_ZC65_R"] = zc65
    df.to_csv(path)
    print(f"[INFO] Patched HEK293T ZC65 in-place: {path}")


def normalize_cells(cells: Iterable[str]) -> list[str]:
    requested = list(cells)
    if not requested:
        return ["Neuron"]
    lowered = {c.lower() for c in requested}
    if "all" in lowered:
        return ["HEK293T", "THP1Macrophage", "THP1Monocyte", "HMC3", "Neuron"]

    invalid = [c for c in requested if c not in VALID_CELLS]
    if invalid:
        raise ValueError(
            f"Invalid --cells value(s): {invalid}. Valid values: all, "
            + ", ".join(sorted(VALID_CELLS))
        )
    return requested


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Make RNA/DNA matched barcode tables. This version adds Neuron_DNA_RNA.csv "
            "support, cleans rich Neuron tables before matching, and allows running Neuron only without patching HEK293T."
        )
    )

    parser.add_argument("--base-dir", default=None, help="Repo/data root. Defaults to auto-detected repo root.")
    parser.add_argument("--counts-dir", default="outputs/read_counts_R1R2", help="Directory containing <CELL>_DNA_RNA.csv")
    parser.add_argument("--indexing-dir", default="indexing", help="Directory containing *_RNA_DNA_MatchTable.csv")
    parser.add_argument(
        "--master-enhancer-map",
        default="outputs/read_counts_R1R2/df_MPRA3_20231110_R1R2.csv",
        help=(
            "CSV with index col 'ID' and column 'enhancer_id'. If this file is missing "
            "and the input table has ID/enhancer_id columns, a fallback map is created."
        ),
    )
    parser.add_argument("--out-dir", default="outputs/read_counts_R1R2", help="Where to write matched outputs.")

    parser.add_argument(
        "--cells",
        nargs="+",
        default=["Neuron"],
        help=(
            "Which built-in cell types to process. Use 'all' for HEK293T, THP1Macrophage, "
            "THP1Monocyte, HMC3, and Neuron. Default: Neuron."
        ),
    )
    parser.add_argument(
        "--neuron-file",
        default="Neuron_DNA_RNA.csv",
        help="Neuron input file name under --counts-dir, or an absolute path. Default: Neuron_DNA_RNA.csv",
    )
    parser.add_argument(
        "--neuron-cell-type",
        default="Neuron",
        help="cell_type value passed to mpra.matching for the Neuron table. Default: Neuron",
    )
    parser.add_argument(
        "--neuron-no-clean-input",
        action="store_true",
        help=(
            "Disable the Neuron-specific clean-input preprocessing step. "
            "Normally you should not use this for Neuron_DNA_RNA.csv because it contains metadata columns."
        ),
    )
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="Skip a requested built-in cell type when its required input files are missing.",
    )
    parser.add_argument(
        "--no-enhancer-map-fallback",
        action="store_true",
        help="Do not create a fallback enhancer map from input ID/enhancer_id columns when master map is missing.",
    )

    parser.add_argument(
        "--patch-hek293t",
        action="store_true",
        help="After matching, patch HEK293T_DNA_matched_barcodes.csv in-place to add Naive_HEK293T_ZC65_R.",
    )
    parser.add_argument(
        "--hek-patch-file",
        default="HEK293T_DNA_matched_barcodes.csv",
        help="Which HEK293T DNA matched file under --out-dir to patch when --patch-hek293t is used.",
    )

    parser.add_argument(
        "--single-dna-rna-file",
        default=None,
        help="Optional single input CSV path for a custom cell type. If set, built-in --cells are ignored.",
    )
    parser.add_argument(
        "--single-cell-type",
        default=None,
        help="cell_type value for --single-dna-rna-file.",
    )
    parser.add_argument(
        "--single-lookup-table",
        default=None,
        help="Optional lookup table path for --single-dna-rna-file. If omitted, auto suffix matching is used.",
    )

    args = parser.parse_args()

    repo_root = Path(args.base_dir).resolve() if args.base_dir else guess_repo_root()
    counts_dir = resolve_path(repo_root, args.counts_dir)
    indexing_dir = resolve_path(repo_root, args.indexing_dir)
    master_map = resolve_path(repo_root, args.master_enhancer_map)
    out_dir = resolve_path(repo_root, args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    allow_fallback = not args.no_enhancer_map_fallback

    print(f"[INFO] repo_root: {repo_root}")
    print(f"[INFO] counts_dir: {counts_dir}")
    print(f"[INFO] indexing_dir: {indexing_dir}")
    print(f"[INFO] master_map: {master_map}")
    print(f"[INFO] out_dir: {out_dir}")

    if args.single_dna_rna_file is not None:
        if args.single_cell_type is None:
            raise ValueError("--single-cell-type is required when --single-dna-rna-file is set.")

        single_file = resolve_path(repo_root, args.single_dna_rna_file)
        if args.single_lookup_table is None:
            run_auto_match(
                dna_rna_file_path=single_file,
                master_map=master_map,
                out_dir=out_dir,
                cell_type=args.single_cell_type,
                rename_dict=None,
                skip_missing=args.skip_missing,
                validate_suffix_pairs=True,
                allow_enhancer_map_fallback=allow_fallback,
            )
        else:
            run_lookup_match(
                dna_rna_file_path=single_file,
                master_map=master_map,
                out_dir=out_dir,
                cell_type=args.single_cell_type,
                lookup_table_path=resolve_path(repo_root, args.single_lookup_table),
                rename_dict=None,
                skip_missing=args.skip_missing,
                allow_enhancer_map_fallback=allow_fallback,
            )
    else:
        cells = normalize_cells(args.cells)

        for cell in cells:
            if cell == "HEK293T":
                hek_rename = {"ZC76RCol-MPRA3-HEK293T-26-Naive": "ZC76R-MPRA3-HEK293T-26-Naive"}
                run_auto_match(
                    dna_rna_file_path=counts_dir / "HEK293T_DNA_RNA.csv",
                    master_map=master_map,
                    out_dir=out_dir,
                    cell_type="HEK293T",
                    rename_dict=hek_rename,
                    skip_missing=args.skip_missing,
                    validate_suffix_pairs=False,
                    allow_enhancer_map_fallback=allow_fallback,
                )

            elif cell == "THP1Macrophage":
                run_lookup_match(
                    dna_rna_file_path=counts_dir / "THP1_DNA_RNA.csv",
                    master_map=master_map,
                    out_dir=out_dir,
                    cell_type="THP1Macrophage",
                    lookup_table_path=indexing_dir / "THP1mac_RNA_DNA_MatchTable.csv",
                    rename_dict=None,
                    skip_missing=args.skip_missing,
                    allow_enhancer_map_fallback=allow_fallback,
                )

            elif cell == "THP1Monocyte":
                run_lookup_match(
                    dna_rna_file_path=counts_dir / "THP1_DNA_RNA.csv",
                    master_map=master_map,
                    out_dir=out_dir,
                    cell_type="THP1Monocyte",
                    lookup_table_path=indexing_dir / "THP1mono_RNA_DNA_MatchTable.csv",
                    rename_dict=None,
                    skip_missing=args.skip_missing,
                    allow_enhancer_map_fallback=allow_fallback,
                )

            elif cell == "HMC3":
                hmc3_rename = {
                    "ZC22-24D-MPRA3-HMC3-26-Stim": "ZC22D-MPRA3-HMC3-26-Naive",
                    "ZC3-5D-MPRA3-HMC3-26-Stim": "ZC3D-MPRA3-HMC3-26-Naive",
                    "ZC72RAco-MPRA3-HMC3-26-Naive": "ZC72R-MPRA3-HMC3Aco-26-Naive",
                    "ZC72RCol-MPRA3-HMC3-26-Naive": "ZC72R-MPRA3-HMC3Col-26-Naive",
                    "ZC73RAco-MPRA3-HMC3-26-IFNG": "ZC73R-MPRA3-HMC3Aco-26-IFNG",
                    "ZC73RCol-MPRA3-HMC3-26-IFNG": "ZC73R-MPRA3-HMC3Col-26-IFNG",
                    "ZC74RAco-MPRA3-HMC3-26-IFNB": "ZC74R-MPRA3-HMC3Aco-26-IFNB",
                    "ZC74RCol-MPRA3-HMC3-26-IFNB": "ZC74R-MPRA3-HMC3Col-26-IFNB",
                    "ZC75RAco-MPRA3-HMC3-26-LPSIFNG": "ZC75R-MPRA3-HMC3Aco-26-LPSIFNG",
                    "ZC75RCol-MPRA3-HMC3-26-LPSIFNG": "ZC75R-MPRA3-HMC3Col-26-LPSIFNG",
                }
                run_lookup_match(
                    dna_rna_file_path=counts_dir / "HMC3_DNA_RNA.csv",
                    master_map=master_map,
                    out_dir=out_dir,
                    cell_type="HMC3",
                    lookup_table_path=indexing_dir / "HMC3_RNA_DNA_MatchTable.csv",
                    rename_dict=hmc3_rename,
                    skip_missing=args.skip_missing,
                    allow_enhancer_map_fallback=allow_fallback,
                )

            elif cell == "Neuron":
                neuron_path = resolve_path(counts_dir, args.neuron_file)
                if args.neuron_no_clean_input:
                    run_auto_match(
                        dna_rna_file_path=neuron_path,
                        master_map=master_map,
                        out_dir=out_dir,
                        cell_type=args.neuron_cell_type,
                        rename_dict=None,
                        skip_missing=args.skip_missing,
                        validate_suffix_pairs=True,
                        allow_enhancer_map_fallback=allow_fallback,
                    )
                else:
                    run_clean_auto_match(
                        dna_rna_file_path=neuron_path,
                        master_map=master_map,
                        out_dir=out_dir,
                        cell_type=args.neuron_cell_type,
                        rename_dict=None,
                        skip_missing=args.skip_missing,
                        validate_suffix_pairs=True,
                        allow_enhancer_map_fallback=allow_fallback,
                    )

    if args.patch_hek293t:
        patch_hek293t_zc65_inplace(
            out_dir=out_dir,
            filename=args.hek_patch_file,
        )

    print("[DONE] Matching pipeline finished.")


if __name__ == "__main__":
    main()
