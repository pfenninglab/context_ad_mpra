from pathlib import Path
import sys


def bootstrap_local_repo_imports() -> None:
    here = Path(__file__).resolve()
    repo_root = here.parents[1]
    src_dir = repo_root / "src"
    for candidate in (src_dir, repo_root):
        candidate_str = str(candidate)
        if candidate_str not in sys.path and (candidate / "mpra").exists():
            sys.path.insert(0, candidate_str)


bootstrap_local_repo_imports()

import pandas as pd

from mpra.expand_annotation import expand_annotation_32barcodes_from_enhancer_annotation, expand_annotation_5barcodes_from_enhancer_annotation
from mpra.make_count_table import BarcodeProcessor


NEURON_SAMPLES = [
    ("Rep1_iPSCNeuron_ZC130_R", 1),
    ("Rep2_iPSCNeuron_ZC131_R", 2),
    ("Rep3_iPSCNeuron_ZC132_R", 3),
    ("Rep4_iPSCNeuron_ZC133_R", 4),
]


def make_neuron_pseudobarcode_enhancer_annotation(save_path: str) -> pd.DataFrame:
    rows = []
    for sample, rep in NEURON_SAMPLES:
        rows.append(
            {
                "sample": sample,
                "Tissue": "iPSCNeuron",
                "Delivery": "MPRA",
                "Date": pd.NA,
                "Animal": f"Rep{rep}",
                "Sex": pd.NA,
                "Test": sample,
                "Function": "Neuron",
                "Parent": "NeuronPseudobarcodes",
                "Neuron": 1,
                "Replicate": rep,
                "RT_reaction": sample,
            }
        )
    annotation = pd.DataFrame(rows).set_index("sample")
    annotation.index.name = None
    annotation.to_csv(save_path)
    return annotation


def make_neuron_enhancer_annotation(save_path: str) -> pd.DataFrame:
    rows = []
    for sample, rep in NEURON_SAMPLES:
        rows.append(
            {
                "sample": sample,
                "Tissue": "iPSCNeuron",
                "Delivery": "MPRA",
                "Date": pd.NA,
                "Animal": f"Rep{rep}",
                "Sex": pd.NA,
                "Test": sample,
                "Function": "Neuron",
                "Parent": "Neuron",
                "Neuron": 1,
                "Replicate": rep,
                "RT_reaction": sample,
            }
        )
    annotation = pd.DataFrame(rows).set_index("sample")
    annotation.index.name = None
    annotation.to_csv(save_path)
    return annotation


def write_altref_barcode_annotation(
    barcode_annotation_path: str,
    count_table_path: str,
    save_path: str,
) -> pd.DataFrame:
    """
    Write the colAnnot that matches *_reshaped_altref.csv count matrices.

    The alt/ref count matrices have only ALT-named columns, because REF rows are
    mapped into enhancer rows while keeping the ALT column namespace. Therefore
    colAnnot must also be ALT-only and ordered exactly like the count columns.
    """
    annotation = pd.read_csv(barcode_annotation_path, index_col=0)
    count_table = pd.read_csv(count_table_path, index_col=0, nrows=0)
    wanted = list(count_table.columns)

    missing = [col for col in wanted if col not in annotation.index]
    if missing:
        raise KeyError(
            f"Missing altref count columns in annotation {barcode_annotation_path}: {missing}"
        )

    altref_annotation = annotation.loc[wanted].copy()
    altref_annotation.to_csv(save_path)
    return altref_annotation


def run_neuron_pseudobarcodes_step06() -> None:
    io_dir = "outputs/read_counts_R1R2"

    processor_neuron = BarcodeProcessor(
        lookup_table_path="indexing/ALT_REF_LookUpTable_amended_pseudobarcode_20240404.csv",
        data_drop_table_path="indexing/low_frequency_enhancer_drop_table_20260126.csv",
        base_dir=Path.cwd(),
    )

    for na in ("RNA", "DNA"):
        processor_neuron.convert_REFALT_Barcodes(
            count_table_path=f"{io_dir}/NeuronPseudobarcodes_{na}_matched_barcodes.csv",
            experiment_name="NeuronPseudobarcodes",
            nucleic_acid_type=na,
            io_dir=io_dir,
            output_dir=io_dir,
        )

    make_neuron_pseudobarcode_enhancer_annotation(
        "annotation_enhancer/mpra3_annot_NeuronPseudobarcodes.csv"
    )
    expand_annotation_5barcodes_from_enhancer_annotation(
        enhancer_annotation_path="annotation_enhancer/mpra3_annot_NeuronPseudobarcodes.csv",
        save_path="annotation_barcodes/mpra3_annot_NeuronPseudobarcodes_barcodes.csv",
    )
    write_altref_barcode_annotation(
        barcode_annotation_path="annotation_barcodes/mpra3_annot_NeuronPseudobarcodes_barcodes.csv",
        count_table_path=f"{io_dir}/NeuronPseudobarcodes_RNA_matched_barcodes_reshaped_altref.csv",
        save_path="annotation_barcodes/mpra3_annot_NeuronPseudobarcodes_altref_barcodes.csv",
    )


def run_neuron_normal_step06() -> None:
    enhancer_annotation_path = "annotation_enhancer/mpra3_annot_Neuron.csv"
    barcode_annotation_path = "annotation_barcodes/mpra3_annot_Neuron_barcodes.csv"

    make_neuron_enhancer_annotation(enhancer_annotation_path)
    expand_annotation_32barcodes_from_enhancer_annotation(
        enhancer_annotation_path=enhancer_annotation_path,
        save_path=barcode_annotation_path,
    )
    write_altref_barcode_annotation(
        barcode_annotation_path=barcode_annotation_path,
        count_table_path="outputs/read_counts_R1R2/Neuron_RNA_matched_barcodes_reshaped_altref.csv",
        save_path="annotation_barcodes/mpra3_annot_Neuron_altref_barcodes.csv",
    )

    neuron_annotation = pd.read_csv(barcode_annotation_path, index_col=0)

    neuron_columns = neuron_annotation.index
    for na in ("RNA", "DNA"):
        count_table = pd.read_csv(
            f"outputs/read_counts_R1R2/Neuron_{na}_matched_barcodes_reshaped.csv",
            index_col=0,
        )
        missing = [col for col in neuron_columns if col not in count_table.columns]
        if missing:
            raise KeyError(f"Missing Neuron annotation columns in {na} reshaped count table: {missing}")
        count_table[neuron_columns].to_csv(
            f"outputs/read_counts_R1R2/Neuron_{na}_matched_barcodes_reshaped.csv"
        )

def main():

    REPO_ROOT = Path(__file__).resolve().parents[1]  # script in scripts/

    ##############################################################HEK293T###########################################################
    expand_annotation_32barcodes_from_enhancer_annotation(enhancer_annotation_path='annotation_enhancer/mpra3_annot_HEK293T.csv', 
        save_path='annotation_barcodes/mpra3_annot_HEK293T_barcodes.csv')

    ##############################################################THP1###########################################################
    #Expand 32 annotation
    expand_annotation_32barcodes_from_enhancer_annotation(enhancer_annotation_path='annotation_enhancer/mpra3_annot_THP1.csv', 
        save_path='annotation_barcodes/mpra3_annot_THP1_barcodes.csv')

    #THP1 Interaction
    import pandas as pd
    THP1_annotation = pd.read_csv('annotation_barcodes/mpra3_annot_THP1_barcodes.csv',index_col=0)

    THP1_annotation['Allele'] = THP1_annotation['Allele_String'].map({'ALT': 1, 'REF': 0})
    THP1_annotation["IFNB_interaction"] = THP1_annotation["IFNB"] * THP1_annotation["Allele"]
    THP1_annotation["IFNG_interaction"] = THP1_annotation["IFNG"] * THP1_annotation["Allele"]
    THP1_annotation["LPSIFNG_interaction"] = THP1_annotation["LPSIFNG"] * THP1_annotation["Allele"]
    THP1_annotation["Naive_interaction"] = THP1_annotation["Naive_Macrophage"] * THP1_annotation["Allele"]
    THP1_annotation.to_csv('annotation_barcodes/mpra3_annot_THP1_barcodes_interaction.csv')

    #Split annotation
    THP1 = pd.read_csv('annotation_barcodes/mpra3_annot_THP1_barcodes_interaction.csv',index_col=0)
    THP1[THP1["IFNB"]==1].to_csv("annotation_barcodes/mpra3_annot_THP1_IFNB_barcodes.csv")
    THP1[THP1["IFNG"]==1].to_csv("annotation_barcodes/mpra3_annot_THP1_IFNG_barcodes.csv")
    THP1[THP1["LPSIFNG"]==1].to_csv("annotation_barcodes/mpra3_annot_THP1_LPSIFNG_barcodes.csv")
    THP1[THP1["Naive_Macrophage"]==1].to_csv("annotation_barcodes/mpra3_annot_THP1_Naive_barcodes.csv")
    THP1[THP1["Tissue"]=="THP1_Monocyte"].to_csv("annotation_barcodes/mpra3_annot_THP1_Monocyte_barcodes.csv")

    THP1 = pd.read_csv('annotation_barcodes/mpra3_annot_THP1_barcodes_interaction.csv',index_col=0)
    THP1[(THP1["IFNB"]==1) | (THP1["Naive_Macrophage"]==1)].to_csv("annotation_barcodes/mpra3_annot_THP1_IFNBvsNaive_barcodes.csv")
    THP1[(THP1["IFNG"]==1) | (THP1["Naive_Macrophage"]==1)].to_csv("annotation_barcodes/mpra3_annot_THP1_IFNGvsNaive_barcodes.csv")
    THP1[(THP1["LPSIFNG"]==1) | (THP1["Naive_Macrophage"]==1)].to_csv("annotation_barcodes/mpra3_annot_THP1_LPSIFNGvsNaive_barcodes.csv")
    THP1[(THP1["LPSIFNG"]==1) | (THP1["IFNG"]==1)].to_csv("annotation_barcodes/mpra3_annot_THP1_LPSIFNGvsIFNG_barcodes.csv")
    THP1[(THP1["Tissue"]=="THP1_Monocyte") | (THP1["Tissue"]=="THP1_Naive")].to_csv("annotation_barcodes/mpra3_annot_THP1_MacrophagevsMonocyte_barcodes.csv")

    #Split barcodes according to annotation
    import pandas as pd
    THP1Macrophage_barcodes = pd.read_csv("outputs/read_counts_R1R2/THP1Macrophage_DNA_matched_barcodes_reshaped.csv",index_col=0)
    THP1Monocyte_barcodes = pd.read_csv("outputs/read_counts_R1R2/THP1Monocyte_DNA_matched_barcodes_reshaped.csv",index_col=0)
    THP1_barcodes = pd.merge(THP1Macrophage_barcodes, THP1Monocyte_barcodes,left_index=True,right_index=True)
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_IFNB_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_IFNB_DNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_IFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_IFNG_DNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_LPSIFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_LPSIFNG_DNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_Naive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_Naive_DNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_Monocyte_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_Monocyte_DNA_matched_barcodes_reshaped.csv")

    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_IFNBvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_IFNBvsNaive_DNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_IFNGvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_IFNGvsNaive_DNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_LPSIFNGvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_LPSIFNGvsNaive_DNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_LPSIFNGvsIFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_LPSIFNGvsIFNG_DNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_MacrophagevsMonocyte_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_MacrophagevsMonocyte_DNA_matched_barcodes_reshaped.csv")

###################################################################################################################################################
    THP1Macrophage_barcodes = pd.read_csv("outputs/read_counts_R1R2/THP1Macrophage_RNA_matched_barcodes_reshaped.csv",index_col=0)
    THP1Monocyte_barcodes = pd.read_csv("outputs/read_counts_R1R2/THP1Monocyte_RNA_matched_barcodes_reshaped.csv",index_col=0)
    THP1_barcodes = pd.merge(THP1Macrophage_barcodes, THP1Monocyte_barcodes,left_index=True,right_index=True)
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_IFNB_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_IFNB_RNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_IFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_IFNG_RNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_LPSIFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_LPSIFNG_RNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_Naive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_Naive_RNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_Monocyte_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_Monocyte_RNA_matched_barcodes_reshaped.csv")


    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_IFNBvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_IFNBvsNaive_RNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_IFNGvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_IFNGvsNaive_RNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_LPSIFNGvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_LPSIFNGvsNaive_RNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_LPSIFNGvsIFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_LPSIFNGvsIFNG_RNA_matched_barcodes_reshaped.csv")
    THP1_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_THP1_MacrophagevsMonocyte_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/THP1_MacrophagevsMonocyte_RNA_matched_barcodes_reshaped.csv")

    ##############################################################HMC3###########################################################
    #Expand 32 annotation
    expand_annotation_32barcodes_from_enhancer_annotation(enhancer_annotation_path='annotation_enhancer/mpra3_annot_HMC3.csv', 
        save_path='annotation_barcodes/mpra3_annot_HMC3_barcodes.csv')

    #Split annotation
    HMC3 = pd.read_csv('annotation_barcodes/mpra3_annot_HMC3_barcodes.csv',index_col=0)
    HMC3[HMC3["IFNB"]==1].to_csv("annotation_barcodes/mpra3_annot_HMC3_IFNB_barcodes.csv")
    HMC3[HMC3["IFNG"]==1].to_csv("annotation_barcodes/mpra3_annot_HMC3_IFNG_barcodes.csv")
    HMC3[HMC3["LPSIFNG"]==1].to_csv("annotation_barcodes/mpra3_annot_HMC3_LPSIFNG_barcodes.csv")
    HMC3[HMC3["Naive_Macrophage"]==1].to_csv("annotation_barcodes/mpra3_annot_HMC3_Naive_barcodes.csv")

    #Split barcodes according to annotation
    import pandas as pd
    HMC3_barcodes = pd.read_csv("outputs/read_counts_R1R2/HMC3_DNA_matched_barcodes_reshaped.csv",index_col=0)
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_IFNB_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_IFNB_DNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_IFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_IFNG_DNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_LPSIFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_LPSIFNG_DNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_Naive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_Naive_DNA_matched_barcodes_reshaped.csv")

    HMC3_barcodes = pd.read_csv("outputs/read_counts_R1R2/HMC3_RNA_matched_barcodes_reshaped.csv",index_col=0)
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_IFNB_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_IFNB_RNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_IFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_IFNG_RNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_LPSIFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_LPSIFNG_RNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_Naive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_Naive_RNA_matched_barcodes_reshaped.csv")

    #HMC3 Interaction
    import pandas as pd
    HMC3_annotation = pd.read_csv('annotation_barcodes/mpra3_annot_HMC3_barcodes.csv',index_col=0)

    HMC3_annotation['Allele'] = HMC3_annotation['Allele_String'].map({'ALT': 1, 'REF': 0})

    HMC3_annotation["IFNB_interaction"] = HMC3_annotation["IFNB"] * HMC3_annotation["Allele"]
    HMC3_annotation["IFNG_interaction"] = HMC3_annotation["IFNG"] * HMC3_annotation["Allele"]
    HMC3_annotation["LPSIFNG_interaction"] = HMC3_annotation["LPSIFNG"] * HMC3_annotation["Allele"]
    HMC3_annotation["Naive_interaction"] = HMC3_annotation["Naive_Macrophage"] * HMC3_annotation["Allele"]
    HMC3_annotation.to_csv('annotation_barcodes/mpra3_annot_HMC3_barcodes_interaction.csv')

    HMC3 = pd.read_csv('annotation_barcodes/mpra3_annot_HMC3_barcodes_interaction.csv',index_col=0)
    HMC3[(HMC3["IFNB"]==1) | (HMC3["Naive_Macrophage"]==1)].to_csv("annotation_barcodes/mpra3_annot_HMC3_IFNBvsNaive_barcodes.csv")
    HMC3[(HMC3["IFNG"]==1) | (HMC3["Naive_Macrophage"]==1)].to_csv("annotation_barcodes/mpra3_annot_HMC3_IFNGvsNaive_barcodes.csv")
    HMC3[(HMC3["LPSIFNG"]==1) | (HMC3["Naive_Macrophage"]==1)].to_csv("annotation_barcodes/mpra3_annot_HMC3_LPSIFNGvsNaive_barcodes.csv")
    HMC3[(HMC3["LPSIFNG"]==1) | (HMC3["IFNG"]==1)].to_csv("annotation_barcodes/mpra3_annot_HMC3_LPSIFNGvsIFNG_barcodes.csv")

    #Split barcodes according to annotation
    import pandas as pd
    HMC3_barcodes = pd.read_csv("outputs/read_counts_R1R2/HMC3_DNA_matched_barcodes_reshaped.csv",index_col=0)
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_IFNBvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_IFNBvsNaive_DNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_IFNGvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_IFNGvsNaive_DNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_LPSIFNGvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_LPSIFNGvsNaive_DNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_LPSIFNGvsIFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_LPSIFNGvsIFNG_DNA_matched_barcodes_reshaped.csv")


    #Split barcodes according to annotation
    import pandas as pd
    HMC3_barcodes = pd.read_csv("outputs/read_counts_R1R2/HMC3_RNA_matched_barcodes_reshaped.csv",index_col=0)
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_IFNBvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_IFNBvsNaive_RNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_IFNGvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_IFNGvsNaive_RNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_LPSIFNGvsNaive_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_LPSIFNGvsNaive_RNA_matched_barcodes_reshaped.csv")
    HMC3_barcodes[pd.read_csv("annotation_barcodes/mpra3_annot_HMC3_LPSIFNGvsIFNG_barcodes.csv",index_col=0).index].to_csv("outputs/read_counts_R1R2/HMC3_LPSIFNGvsIFNG_RNA_matched_barcodes_reshaped.csv")

    ##############################################################Neuron normal barcode version###########################################################
    run_neuron_normal_step06()

    ##############################################################Neuron pseudobarcodes###########################################################
    run_neuron_pseudobarcodes_step06()

    


if __name__ == "__main__":
    main()
