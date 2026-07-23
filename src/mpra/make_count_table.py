import pandas as pd
import numpy as np
from pathlib import Path


class BarcodeProcessor:
    """
    Same logic as your original implementation, but paths are now robust:
    - Accepts str or Path everywhere.
    - Uses a configurable base_dir (repo/data root) so relative paths don't depend on CWD.
    - Keeps your output/lookup conventions (read_counts_R1R2/<experiment>_<type>_...).
    """

    def __init__(self, lookup_table_path, data_drop_table_path, base_dir=None):
        # base_dir: where relative paths (like "read_counts_R1R2/...") are anchored
        self.base_dir = Path(base_dir).resolve() if base_dir is not None else Path.cwd().resolve()

        self.lookup_table = self.load_lookup_table(self._resolve_path(lookup_table_path))
        self.filter_lookup_table(self._resolve_path(data_drop_table_path))  # lookup_table filtering

    def _resolve_path(self, p):
        """
        Resolve a path flexibly:
        - If p is absolute: use it.
        - If relative: interpret relative to base_dir.
        """
        p = Path(p)
        return p if p.is_absolute() else (self.base_dir / p)

    def load_lookup_table(self, file_path):
        """Load lookup table from CSV file."""
        return pd.read_csv(file_path, header=None)

    def filter_lookup_table(self, file_path):
        """Filter lookup table based on criteria defined in another CSV."""
        df_to_drop = pd.read_csv(file_path, index_col=0)
        strings_to_drop = df_to_drop[df_to_drop["Drop"] == "Y"].index
        for string_to_drop in strings_to_drop:
            self.lookup_table = self.lookup_table[
                ~self.lookup_table.apply(
                    lambda row: row.astype(str).str.contains(string_to_drop).any(), axis=1
                )
            ]

    def process_columns(self, count_table):
        """Process columns to separate and calculate 'REF' and 'ALT' data."""
        count_table_new = pd.DataFrame()

        # Exclude 'seq_id' column explicitly
        columns_to_process = [col for col in count_table.columns if col != "seq_id"]

        for col in columns_to_process:
            df_tmp = count_table[[col, "seq_id"]]
            ref = df_tmp.loc[self.lookup_table[0]].set_index("seq_id")
            alt = df_tmp.loc[self.lookup_table[1]].set_index("seq_id")
            count_table_new[col + "_ALT"] = ref[col]
            count_table_new[col + "_REF"] = alt[col].tolist()

        return count_table_new

    def reshape_group(self, group, samples, max_barcodes):
        """Reshape the grouped data for each enhancer."""
        n_samples = len(samples)
        n_barcodes = group.shape[0]  # Number of rows in the group
        reshaped_data = {}

        for i in range(n_samples):
            for j in range(max_barcodes):
                column_name = f"{samples[i]}_Barcode_{j+1}"
                reshaped_data[column_name] = group.iloc[j, i + 1] if j < n_barcodes else np.nan
        return pd.Series(reshaped_data)

    def convert_REFALT_Barcodes(
        self,
        count_table_path,
        experiment_name,
        nucleic_acid_type,
        io_dir="read_counts_R1R2",
        output_dir=None,
    ):
        """
        Convert REF/ALT barcodes for a given nucleic acid type and experiment.

        Flexibility additions (logic unchanged):
        - count_table_path can be absolute or relative to base_dir.
        - lookup/output are still constructed from experiment_name/nucleic_acid_type,
          but io_dir/output_dir can be set to move data around.
        """
        io_dir = self._resolve_path(io_dir)
        output_dir = self._resolve_path(output_dir) if output_dir is not None else io_dir

        # Construct file paths (same naming convention as your original)
        lookup_file_path = io_dir / f"{experiment_name}_{nucleic_acid_type}_matched_barcodes.csv"
        output_file_path = output_dir / f"{experiment_name}_{nucleic_acid_type}_matched_barcodes_reshaped.csv"

        # Load count table and set index
        count_table_path = self._resolve_path(count_table_path)
        count_table = pd.read_csv(count_table_path, index_col=0)
        count_table["seq_id"] = count_table.index
        count_table.set_index("enhancer_id", inplace=True)

        # Process columns
        processed_count_table = self.process_columns(count_table)

        # Load lookup table, rename axis, and merge with processed count table
        lookup_table = pd.read_csv(lookup_file_path, index_col=0)["enhancer_id"]
        lookup_table = lookup_table.rename_axis("seq_id").reset_index()
        annotated_df = processed_count_table.merge(lookup_table, on="seq_id", how="left")

        # Reshape data
        max_barcodes = annotated_df.groupby("enhancer_id").size().max()
        columns_to_keep = [col for col in annotated_df.columns if col not in ["seq_id", "enhancer_id"]]
        reshaped_df = annotated_df.groupby("enhancer_id").apply(
            lambda group: self.reshape_group(group, columns_to_keep, max_barcodes)
        )

        # Save reshaped dataframe
        output_file_path.parent.mkdir(parents=True, exist_ok=True)
        reshaped_df.fillna(0).to_csv(output_file_path)
        return reshaped_df

