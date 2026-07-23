# R1/R2 MPRA read mapping with selectable mapping mode
import os
import pandas as pd
import gzip
from Bio import SeqIO
from collections import Counter
from typing import List, Dict, Tuple, Optional
import numpy as np


class MPRA:
    def __init__(self, test_library: str, flank_bp: int = 10, match_mode: str = "both") -> None:
        """
        Constructor to initialize MPRA object with test library.

        Parameters:
        test_library (str): Path to the test library CSV file.
        flank_bp (int): Number of enhancer bases to use for mapping.
                        Original behavior is 10. Use 5 for 5 bp mapping.
        match_mode (str): Mapping key mode.
                          "both": barcode + R2 flank + R1 flank. Original behavior.
                          "r1":   barcode + R1 flank only.
                          "r2":   barcode + R2 flank only.
        """
        if not isinstance(flank_bp, int) or flank_bp <= 0:
            raise ValueError("flank_bp must be a positive integer.")

        allowed_modes = {"both", "r1", "r2"}
        if match_mode not in allowed_modes:
            raise ValueError(f"match_mode must be one of {sorted(allowed_modes)}.")

        self.test_library = pd.read_csv(test_library)
        self.flank_bp = flank_bp
        self.match_mode = match_mode
        self.match_col = self._get_match_col()
        self._ensure_match_column()

    def _get_match_col(self) -> str:
        """
        Return the library/count-table column name used for mapping.
        """
        if self.match_mode == "both":
            return f"Barcode_RC_Enhancer{self.flank_bp}_R2_R1"
        if self.match_mode == "r1":
            return f"Barcode_RC_Enhancer{self.flank_bp}_R1"
        if self.match_mode == "r2":
            return f"Barcode_RC_Enhancer{self.flank_bp}_R2"
        raise ValueError(f"Unsupported match_mode: {self.match_mode}")

    def _ensure_match_column(self) -> None:
        """
        Ensure the library has the key column used for read mapping.

        If the requested key column does not exist, derive it from
        Barcode_RC_Enhancer10_R2_R1 by truncating the R2 and/or R1 flanks.
        """
        if self.match_col in self.test_library.columns:
            self._warn_if_duplicate_keys()
            return

        source_col = "Barcode_RC_Enhancer10_R2_R1"

        if source_col not in self.test_library.columns:
            raise KeyError(
                f"Missing '{self.match_col}' and cannot derive it because "
                f"'{source_col}' is not in the library table."
            )

        if self.flank_bp > 10:
            raise ValueError(
                f"Cannot derive {self.match_col} from {source_col}; "
                "the source column only contains 10 bp flanks."
            )

        parts = self.test_library[source_col].astype(str).str.split("_", n=2, expand=True)

        if parts.shape[1] != 3:
            raise ValueError(
                f"Column '{source_col}' should have format barcode_R2flank_R1flank."
            )

        barcode = parts[0]
        r2_flank = parts[1].str[:self.flank_bp]
        r1_flank = parts[2].str[:self.flank_bp]

        if self.match_mode == "both":
            self.test_library[self.match_col] = barcode + "_" + r2_flank + "_" + r1_flank
        elif self.match_mode == "r1":
            self.test_library[self.match_col] = barcode + "_" + r1_flank
        elif self.match_mode == "r2":
            self.test_library[self.match_col] = barcode + "_" + r2_flank
        else:
            raise ValueError(f"Unsupported match_mode: {self.match_mode}")

        self._warn_if_duplicate_keys()

    def _warn_if_duplicate_keys(self) -> None:
        """
        Warn if the selected mapping key is not unique in the library.
        """
        n_duplicate_keys = self.test_library[self.match_col].duplicated().sum()

        if n_duplicate_keys > 0:
            print(
                f"Warning: {n_duplicate_keys} duplicated keys found in {self.match_col}. "
                "Using fewer bases or only one flank can collapse distinct constructs into the same mapping key."
            )

    @staticmethod
    def _count_frequency(my_list: List[str]) -> Dict[str, int]:
        """
        Count the frequency of each item in a list.
        """
        return dict(Counter(my_list))

    @staticmethod
    def _parse_fastq_to_list(file: str) -> List[str]:
        """
        Parse a gzip-compressed FASTQ file and return a list of sequences.

        This function checks the gzip magic bytes first, so if a macOS
        AppleDouble file such as ._sample_R1_001.fastq.gz is accidentally read,
        the error message will identify the exact bad file.
        """
        with open(file, "rb") as test_handle:
            magic = test_handle.read(2)

        if magic != b"\x1f\x8b":
            raise gzip.BadGzipFile(
                f"Not a valid gzip FASTQ file: {file}. "
                f"First two bytes are {magic!r}, expected b'\\x1f\\x8b'. "
                "This is often caused by macOS hidden files such as ._*.fastq.gz."
            )

        with gzip.open(file, "rt") as handle:
            return [str(record.seq) for record in SeqIO.parse(handle, "fastq")]

    def _find_barcodes(
        self,
        i: str,
        locate_seq_ad: str,
        locate_seq_control: str,
    ) -> Tuple[Optional[str], Optional[str]]:
        """
        Find barcode and R1 enhancer flank sequences in a read1 string.

        Returns:
        Tuple[Optional[str], Optional[str]]: barcode and R1 enhancer flank.
        """
        for j in range(39, 42):
            if i[j:j + 5] in {locate_seq_ad, locate_seq_control}:
                barcode = i[j + 5:j + 21]
                enhancer_seq_R1 = i[j + 33:j + 33 + self.flank_bp]
                return str(barcode), str(enhancer_seq_R1)
        return None, None

    def _find_read2_enhancer(
        self,
        read2: str,
        barcode: str,
        enhancer_seq_R1: str,
    ) -> Optional[str]:
        """
        Find R2 enhancer flank sequence in read2 and build the mapping key.

        Output key format depends on self.match_mode:
        both: barcode_R2flank_R1flank
        r2:   barcode_R2flank
        """
        for j in [0, -1, 1]:
            if read2[46 + j:51 + j] == "TGCTA":
                read2_enhancer_seq = read2[j + 51:j + 51 + self.flank_bp]

                if self.match_mode == "both":
                    return barcode + "_" + read2_enhancer_seq + "_" + enhancer_seq_R1

                if self.match_mode == "r2":
                    return barcode + "_" + read2_enhancer_seq

                raise ValueError("_find_read2_enhancer should not be called for match_mode='r1'.")

        return None

    def _build_r1_mapping_key(self, barcode: str, enhancer_seq_R1: str) -> str:
        """
        Build an R1-only mapping key.

        Output key format:
        barcode_R1flank
        """
        return barcode + "_" + enhancer_seq_R1

    def _process_reads(
        self,
        read1: List[str],
        read2: Optional[List[str]] = None,
    ) -> Tuple[List[str], int]:
        """
        Process reads and extract mapping keys.

        For match_mode='r1', only read1 is required.
        For match_mode='both' or 'r2', read2 is required.

        Returns:
        Tuple[List[str], int]: mapping keys and total read1 count.
        """
        barcode_list = []

        if self.match_mode != "r1" and read2 is None:
            raise ValueError("read2 is required when match_mode is 'both' or 'r2'.")

        if read2 is not None and len(read1) != len(read2):
            print("read1 != read2")

        for i in range(len(read1)):
            current_read1 = read1[i]
            barcode, enhancer_seq_R1 = self._find_barcodes(
                current_read1,
                locate_seq_ad="TCTAA",
                locate_seq_control="TCTGA",
            )

            if not barcode:
                continue

            if self.match_mode == "r1":
                barcode_list.append(self._build_r1_mapping_key(barcode, enhancer_seq_R1))
                continue

            current_read2 = read2[i]
            read2_barcode_enhancer = self._find_read2_enhancer(
                current_read2,
                barcode,
                enhancer_seq_R1,
            )
            if read2_barcode_enhancer:
                barcode_list.append(read2_barcode_enhancer)

        return barcode_list, len(read1)

    def _frequency_to_table(self, barcode_list: List[str]) -> pd.DataFrame:
        """
        Convert a list of mapping keys to a DataFrame with frequencies.
        """
        count_dict = self._count_frequency(barcode_list)
        count_table = pd.DataFrame.from_dict(count_dict, orient="index", columns=["read_count"])
        count_table[self.match_col] = count_table.index
        return count_table

    def _annotate_AD_barcodes(
        self,
        sample_name: str,
        read1: List[str],
        read2: Optional[List[str]],
        df_AD: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Annotate AD/control barcodes from one FASTQ sample.
        """
        barcode_list, total_read = self._process_reads(read1, read2)
        count_table = self._frequency_to_table(barcode_list)

        df_AD = pd.merge(
            count_table,
            df_AD,
            how="right",
            on=self.match_col,
        ).rename(columns={"read_count": sample_name})

        mapped_reads = df_AD[sample_name].sum()
        mapped_fraction = round(mapped_reads / total_read, 2) if total_read else 0

        print(
            sample_name,
            "Percentage of mapped reads",
            mapped_fraction,
            "total read number",
            total_read / 1000000,
        )

        return df_AD

    @staticmethod
    def _delete_files(file_path: str) -> None:
        """
        Delete a file if it exists.
        """
        if os.path.exists(file_path):
            os.remove(file_path)
        else:
            print(f"{file_path} does not exist")

    def _get_read1_read2(self, dir_folder: str) -> Dict[str, Dict[str, Optional[str]]]:
        """
        Get read1 and read2 files from a directory.

        Skips macOS hidden files such as ._*.fastq.gz and .DS_Store.
        """
        if not os.path.exists(dir_folder):
            print(f"The directory {dir_folder} does not exist")
            return {}

        r1_suffix = "_R1_001.fastq.gz"
        r2_suffix = "_R2_001.fastq.gz"
        paired_files = {}

        for file in sorted(os.listdir(dir_folder)):
            full_path = os.path.join(dir_folder, file)

            if not os.path.isfile(full_path):
                continue

            if file.startswith(".") or file.startswith("._"):
                continue

            if not file.endswith(".fastq.gz"):
                continue

            if file.endswith(r1_suffix):
                sample_name = file[:-len(r1_suffix)]
                read_key = "R1"
            elif file.endswith(r2_suffix):
                sample_name = file[:-len(r2_suffix)]
                read_key = "R2"
            else:
                print(f"Skipping unrecognized FASTQ filename: {file}")
                continue

            if sample_name not in paired_files:
                paired_files[sample_name] = {"R1": None, "R2": None}

            paired_files[sample_name][read_key] = file

        return paired_files

    def read_files_for_sample(
        self,
        sample_name: str,
        files: Dict[str, Optional[str]],
        dir_folder: str,
    ) -> Tuple[Optional[List[str]], Optional[List[str]]]:
        """
        Read FASTQ files for a sample and return read1/read2 sequence lists.
        For match_mode='r1', read2 is optional and will not be loaded.
        """
        file_read1 = os.path.join(dir_folder, files["R1"]) if files["R1"] else None
        file_read2 = os.path.join(dir_folder, files["R2"]) if files["R2"] else None

        if not file_read1:
            print(f"Missing R1 file for sample {sample_name}")
            return None, None

        if self.match_mode != "r1" and not file_read2:
            print(f"Missing R2 file for sample {sample_name}")
            return None, None

        read1 = self._parse_fastq_to_list(file_read1)
        read2 = None

        if self.match_mode != "r1":
            read2 = self._parse_fastq_to_list(file_read2)

        if not read1:
            print("read1 is empty")

        if self.match_mode != "r1" and not read2:
            print("read2 is empty")

        return read1, read2

    def loop_through_samples_and_build_count_tables(self, dir_folder: str) -> None:
        """
        Loop through all samples in a directory, build count tables, and save them as CSV files.
        """
        paired_files = self._get_read1_read2(dir_folder)
        df = self.test_library

        if not paired_files:
            print(f"No FASTQ files found in {dir_folder}")
            return

        print(f"Using flank_bp = {self.flank_bp}")
        print(f"Using match_mode = {self.match_mode}")
        print(f"Using mapping column = {self.match_col}")

        for n, (sample_name, files) in enumerate(paired_files.items(), start=1):
            read1, read2 = self.read_files_for_sample(sample_name, files, dir_folder)

            if read1 is None:
                print(sample_name, "Read1 is None")
                continue

            if self.match_mode != "r1" and read2 is None:
                print(sample_name, "Read2 is None")
                continue

            df = self._annotate_AD_barcodes(sample_name, read1, read2, df)

            if "enhancer_id" in df.columns:
                df = df.sort_values(by="enhancer_id", ascending=False)

            df.to_csv(f"df_MPRA3_{n}.csv", index=False)

            if n > 1:
                self._delete_files(f"df_MPRA3_{n - 1}.csv")
