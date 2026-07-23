import pandas as pd
def expand_annotation_32barcodes_from_enhancer_annotation(enhancer_annotation_path, save_path):
    df = pd.read_csv(enhancer_annotation_path,index_col=0)

    # Assuming 'df' is your original DataFrame
    # Duplicate the DataFrame
    duplicated_df = pd.concat([df, df]).sort_index(kind='merge')

    # Create an alternating index
    new_index = [f"{index}_ALT" if i % 2 == 0 else f"{index}_REF" for i, index in enumerate(duplicated_df.index)]

    # Assign the new index to the DataFrame
    duplicated_df.index = new_index

    # Now duplicated_df has each row duplicated with alternating "_ALT" and "_REF" suffixes in the index
    df = duplicated_df

    # Expand dataframe to have 32 rows for each original row
    expanded_df = pd.DataFrame()

    for _, row in df.iterrows():
        for i in range(1, 33):
            new_row = row.copy()
            new_row['Barcode_number'] = i
            expanded_df = pd.concat([expanded_df, new_row], axis=1)
    expanded_df = expanded_df.T
    expanded_df["sample_barcodes"] = (expanded_df.index + "_Barcode_" + expanded_df["Barcode_number"].astype(str)).tolist()

    expanded_df["Barcode_Allele"] = expanded_df["Barcode_number"].astype(str) + expanded_df.index.str.slice(-4,)
    expanded_df["Allele_String"] = expanded_df.index.str.slice(-3,)
    expanded_df["Test_Allele"] = expanded_df["Test"]+"_"+expanded_df["Allele_String"]
    expanded_df["Animal_Allele"] = expanded_df["Animal"]+"_"+expanded_df["Allele_String"]
    expanded_df["Tissue_Allele"] = expanded_df["Tissue"]+"_"+expanded_df["Allele_String"]
    #expanded_df = expanded_df.reset_index()
    expanded_df = expanded_df.set_index("sample_barcodes")
    expanded_df.rename(columns={"index": "sample_allele"}, inplace=True)
    expanded_df.index.name = None
    if save_path:
        expanded_df.to_csv(save_path)
    return expanded_df

import pandas as pd
def expand_annotation_5barcodes_from_enhancer_annotation(enhancer_annotation_path, save_path):
    df = pd.read_csv(enhancer_annotation_path,index_col=0)

    # Assuming 'df' is your original DataFrame
    # Duplicate the DataFrame
    duplicated_df = pd.concat([df, df]).sort_index(kind='merge')

    # Create an alternating index
    new_index = [f"{index}_ALT" if i % 2 == 0 else f"{index}_REF" for i, index in enumerate(duplicated_df.index)]

    # Assign the new index to the DataFrame
    duplicated_df.index = new_index

    # Now duplicated_df has each row duplicated with alternating "_ALT" and "_REF" suffixes in the index
    df = duplicated_df

    # Expand dataframe to have 32 rows for each original row
    expanded_df = pd.DataFrame()

    for _, row in df.iterrows():
        for i in range(1, 6):
            new_row = row.copy()
            new_row['Barcode_number'] = i
            expanded_df = pd.concat([expanded_df, new_row], axis=1)
    expanded_df = expanded_df.T
    expanded_df["sample_barcodes"] = (expanded_df.index + "_Barcode_" + expanded_df["Barcode_number"].astype(str)).tolist()

    expanded_df["Barcode_Allele"] = expanded_df["Barcode_number"].astype(str) + expanded_df.index.str.slice(-4,)
    expanded_df["Allele_String"] = expanded_df.index.str.slice(-3,)
    expanded_df["Test_Allele"] = expanded_df["Test"]+"_"+expanded_df["Allele_String"]
    expanded_df["Animal_Allele"] = expanded_df["Animal"]+"_"+expanded_df["Allele_String"]
    expanded_df["Tissue_Allele"] = expanded_df["Tissue"]+"_"+expanded_df["Allele_String"]
    #expanded_df = expanded_df.reset_index()
    expanded_df = expanded_df.set_index("sample_barcodes")
    expanded_df.rename(columns={"index": "sample_allele"}, inplace=True)
    expanded_df.index.name = None
    if save_path:
        expanded_df.to_csv(save_path)
    return expanded_df

import pandas as pd
def expand_annotation_3barcodes_from_enhancer_annotation(enhancer_annotation_path, save_path):
    df = pd.read_csv(enhancer_annotation_path,index_col=0)

    # Assuming 'df' is your original DataFrame
    # Duplicate the DataFrame
    duplicated_df = pd.concat([df, df]).sort_index(kind='merge')

    # Create an alternating index
    new_index = [f"{index}_ALT" if i % 2 == 0 else f"{index}_REF" for i, index in enumerate(duplicated_df.index)]

    # Assign the new index to the DataFrame
    duplicated_df.index = new_index

    # Now duplicated_df has each row duplicated with alternating "_ALT" and "_REF" suffixes in the index
    df = duplicated_df

    # Expand dataframe to have 32 rows for each original row
    expanded_df = pd.DataFrame()

    for _, row in df.iterrows():
        for i in range(1, 4):
            new_row = row.copy()
            new_row['Barcode_number'] = i
            expanded_df = pd.concat([expanded_df, new_row], axis=1)
    expanded_df = expanded_df.T
    expanded_df["sample_barcodes"] = (expanded_df.index + "_Barcode_" + expanded_df["Barcode_number"].astype(str)).tolist()

    expanded_df["Barcode_Allele"] = expanded_df["Barcode_number"].astype(str) + expanded_df.index.str.slice(-4,)
    expanded_df["Allele_String"] = expanded_df.index.str.slice(-3,)
    expanded_df["Test_Allele"] = expanded_df["Test"]+"_"+expanded_df["Allele_String"]
    expanded_df["Animal_Allele"] = expanded_df["Animal"]+"_"+expanded_df["Allele_String"]
    expanded_df["Tissue_Allele"] = expanded_df["Tissue"]+"_"+expanded_df["Allele_String"]
    #expanded_df = expanded_df.reset_index()
    expanded_df = expanded_df.set_index("sample_barcodes")
    expanded_df.rename(columns={"index": "sample_allele"}, inplace=True)
    expanded_df.index.name = None
    if save_path:
        expanded_df.to_csv(save_path)
    return expanded_df
