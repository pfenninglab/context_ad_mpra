from pathlib import Path
from mpra.make_count_table import BarcodeProcessor

def main():
    REPO_ROOT = Path(__file__).resolve().parents[1]  # script in scripts/

    # 统一用一个 io_dir（相对 REPO_ROOT）
    IO_DIR = "outputs/read_counts_R1R2"

    # ---- BrainR1R2merged20240404 ----
    processor_brain = BarcodeProcessor(
        lookup_table_path="indexing/ALT_REF_LookUpTable_amended_pseudobarcode_20240404.csv",
        data_drop_table_path="indexing/low_frequency_enhancer_drop_table_20260126.csv",
        base_dir=REPO_ROOT,
    )

    processor_brain.convert_REFALT_Barcodes(
        count_table_path=f"{IO_DIR}/BrainR1R2merged20240404_RNA_matched_barcodes.csv",
        experiment_name="BrainR1R2merged20240404",
        nucleic_acid_type="RNA",
        io_dir=IO_DIR,
        output_dir=IO_DIR,
    )

    # ---- HEK/THP1/HMC3 ----
    processor_cells = BarcodeProcessor(
        lookup_table_path="indexing/ALT_REF_LookUpTable_amended_20231117.csv",
        data_drop_table_path="indexing/low_frequency_enhancer_drop_table_20260126.csv",
        base_dir=REPO_ROOT,
    )

    jobs = [
        ("HEK293T", "RNA"),
        ("HEK293T", "DNA"),
        ("THP1Macrophage",    "RNA"),
        ("THP1Macrophage",    "DNA"),
        ("THP1Monocyte",    "RNA"),
        ("THP1Monocyte",    "DNA"),
        ("HMC3",    "RNA"),
        ("HMC3",    "DNA"),
    ]

    for exp, na in jobs:
        processor_cells.convert_REFALT_Barcodes(
            count_table_path=f"{IO_DIR}/{exp}_{na}_matched_barcodes.csv",
            experiment_name=exp,
            nucleic_acid_type=na,
            io_dir=IO_DIR,
            output_dir=IO_DIR,
        )

if __name__ == "__main__":
    main()

