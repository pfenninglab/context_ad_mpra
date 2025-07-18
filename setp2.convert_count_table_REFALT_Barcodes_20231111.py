import pandas as pd
import numpy as np

class BarcodeProcessor:
    def __init__(self, lookup_table_path,data_drop_table_path):
        self.lookup_table = self.load_lookup_table(lookup_table_path) 
        self.filter_lookup_table(data_drop_table_path) # lookup_table filtering

    def load_lookup_table(self, file_path):
        """
        Load lookup table from CSV file.
        """
        return pd.read_csv(file_path, header=None)

    def filter_lookup_table(self, file_path):
        """
        Filter lookup table based on criteria defined in another CSV.
        """
        df_to_drop = pd.read_csv(file_path, index_col=0)
        strings_to_drop = df_to_drop[df_to_drop["Drop"] == "Y"].index
        for string_to_drop in strings_to_drop:
            self.lookup_table = self.lookup_table[~self.lookup_table.apply(lambda row: row.astype(str).str.contains(string_to_drop).any(), axis=1)]


    def process_columns(self, count_table):
        """
        Process columns to separate and calculate 'REF' and 'ALT' data.
        """
        count_table_new = pd.DataFrame()

        # Exclude 'seq_id' column explicitly
        columns_to_process = [col for col in count_table.columns if col != 'seq_id']

        for col in columns_to_process:
            df_tmp = count_table[[col, "seq_id"]]
            ref = df_tmp.loc[self.lookup_table[0]].set_index("seq_id")
            alt = df_tmp.loc[self.lookup_table[1]].set_index("seq_id")
            count_table_new[col + "_ALT"] = ref[col]
            count_table_new[col + "_REF"] = alt[col].tolist()
        
        return count_table_new

    def reshape_group(self, group, samples, max_barcodes):
        """
        Reshape the grouped data for each enhancer.
        """
        n_samples = len(samples)
        n_barcodes = group.shape[0]  # Number of rows in the group
        reshaped_data = {}
        
        for i in range(n_samples):
            for j in range(max_barcodes):
                column_name = f'{samples[i]}_Barcode_{j+1}'
                reshaped_data[column_name] = group.iloc[j, i + 1] if j < n_barcodes else np.nan
        return pd.Series(reshaped_data)

    def convert_REFALT_Barcodes(self, count_table_path, experiment_name, nucleic_acid_type):
        """
        Convert REF/ALT barcodes for a given nucleic acid type and experiment.
        """
        # Construct file paths
        lookup_file_path = f"read_counts_R1R2/{experiment_name}_{nucleic_acid_type}_matched_barcodes.csv"
        output_file_path = f"read_counts_R1R2/{experiment_name}_{nucleic_acid_type}_matched_barcodes_reshaped.csv"

        # Load count table and set index
        count_table = pd.read_csv(count_table_path, index_col=0)
        count_table['seq_id'] = count_table.index
        count_table.set_index('enhancer_id', inplace=True)

        # Process columns
        processed_count_table = self.process_columns(count_table)

        # Load lookup table, rename axis, and merge with processed count table
        lookup_table = pd.read_csv(lookup_file_path, index_col=0)['enhancer_id']
        lookup_table = lookup_table.rename_axis('seq_id').reset_index()
        annotated_df = processed_count_table.merge(lookup_table, on="seq_id", how='left')

        # Reshape data
        max_barcodes = annotated_df.groupby('enhancer_id').size().max()
        columns_to_keep = [col for col in annotated_df.columns if col not in ['seq_id', 'enhancer_id']]
        reshaped_df = annotated_df.groupby('enhancer_id').apply(lambda group: self.reshape_group(group, columns_to_keep, max_barcodes))

        # Save reshaped dataframe
        reshaped_df.fillna(0).to_csv(output_file_path)
        return reshaped_df


# Example usage of the class
# Initialize the BarcodeProcessor
processor = BarcodeProcessor(lookup_table_path="indexing/ALT_REF_LookUpTable_20231111.csv",data_drop_table_path = "read_counts_R1R2/low_frequency_enhancer_drop_table_20231111.csv")

# Convert REF/ALT barcodes
reshaped_df = processor.convert_REFALT_Barcodes(count_table_path="read_counts_R1R2/THP1_RNA_matched_barcodes.csv", experiment_name="THP1", nucleic_acid_type="DNA")