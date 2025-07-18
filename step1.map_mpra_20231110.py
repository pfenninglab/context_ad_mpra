#R1+R2

import os
import pandas as pd
import gzip
from Bio import SeqIO
from collections import Counter
from typing import List, Dict, Tuple
import numpy as np

class MPRA:
    def __init__(self, test_library: str) -> None:
        """
        Constructor to initialize MPRA object with test library.
        
        Parameters:
        test_library (str): Path to the test library CSV file.
        """
        self.test_library = pd.read_csv(test_library)
        self.Barcode_Enhancer10_R2 = self.test_library["Barcode_RC_Enhancer10_R2"].values

    @staticmethod
    def _count_frequency(my_list: List[str]) -> Dict[str, int]:
        """
        Count the frequency of each item in a list.
        
        Parameters:
        my_list (List[str]): List of strings whose frequencies need to be counted.
        
        Returns:
        Dict[str, int]: Dictionary with items as keys and their frequencies as values.
        """
        return dict(Counter(my_list))

    @staticmethod
    def _parse_fastq_to_list(file: str) -> List[str]:
        """
        Parse a gzip compressed fastq file and return a list of sequences.
        
        Parameters:
        file (str): Path to the gzip compressed fastq file.
        
        Returns:
        List[str]: List of sequence strings from the fastq file.
        """
        with gzip.open(file, "rt") as handle:
            return [str(record.seq) for record in SeqIO.parse(handle, "fastq")]

    def _find_barcodes(self, i: str, locate_seq_ad: str, locate_seq_control:str) -> str:
        """
        Find barcode sequences in a string.
        
        Parameters:
        i (str): Input string where barcodes need to be located.
        locate_seq_ad (str): Sequence used for locating AD barcodes.
        locate_seq_control (str): Sequence used for locating control barcodes.
        
        Returns:
        str: Found barcode sequence.
        """
        for j in range(39, 42):
            if i[j:j+5] in {locate_seq_ad, locate_seq_control}:
                barcode = i[j+5:j+21]
                enhancer_seq_R1 = i[j+33:j+43]
                return str(barcode), str(enhancer_seq_R1)
        return None, None
    
    def _find_read2_enhancer(self, read2: str,barcode:str, enhancer_seq_R1:str) -> str:
        """
        Find enhancer sequence in read2 and concatenate with barcode.
        
        Parameters:
        read2 (str): Read2 sequence string.
        barcode (str): Barcode sequence string.
        
        Returns:
        str: Concatenated barcode and enhancer sequence.
        """
        for j in [0,-1,1]:
            if read2[46+j:51+j]=="TGCTA":
                read2_enhancer_seq = read2[j+51:j+61]
                read2_barcode_enhancer = barcode + "_" + read2_enhancer_seq+"_"+enhancer_seq_R1
                return read2_barcode_enhancer
        return None
    
    def _process_reads(self, read1, read2):
        """
        Process list of reads and classify barcodes.
        
        Parameters:
        read1 (List[str]): List of read1 sequences.
        read2 (List[str]): List of read2 sequences.
        
        Returns:
        List[str]: List of processed barcode sequences.
        """
        barcode_list = []

        if len(read1) != len(read2):
            print("read1 != read2")

        for i in range(len(read1)):
            current_read1 = read1[i]
            barcode, enhancer_seq_R1 = self._find_barcodes(current_read1, locate_seq_ad='TCTAA', locate_seq_control='TCTGA')
            if barcode:
                current_read2 = read2[i]
                read2_barcode_enhancer = self._find_read2_enhancer(current_read2, barcode,enhancer_seq_R1)
                barcode_list.append(read2_barcode_enhancer)
        return barcode_list,len(read1)

    def _frequency_to_table(self, barcode_list: List[str]) -> pd.DataFrame:
        """
        Convert a list of barcodes to a DataFrame with their frequencies.
        
        Parameters:
        barcode_list (List[str]): List of barcode sequences.
        
        Returns:
        pd.DataFrame: DataFrame with barcode frequencies.
        """
        count_dict = self._count_frequency(barcode_list)
        count_table = pd.DataFrame.from_dict(count_dict, orient="index", columns=["read_count"])
        count_table["Barcode_RC_Enhancer10_R2_R1"] = count_table.index
        return count_table

    def _annotate_AD_barcodes(self, sample_name: str, read1,read2, df_AD: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, str]]:
        """
        Annotate AD and control barcodes from a given fastq file.
        
        Parameters:
        sample_name (str): Name of the sample being processed.
        read1 (List[str]): List of read1 sequences.
        read2 (List[str]): List of read2 sequences.
        df_AD (pd.DataFrame): DataFrame containing AD barcode information.
        
        Returns:
        Tuple: Tuple containing DataFrame of annotated AD barcodes, DataFrame of annotated control barcodes,
               and a dictionary mapping barcode sequences to their corresponding AD and control IDs.
        """
        barcode_list,total_read = self._process_reads(read1,read2)
        count_table = self._frequency_to_table(barcode_list)
        df_AD = pd.merge(count_table, df_AD, how="right", on="Barcode_RC_Enhancer10_R2_R1").rename(columns={"read_count": sample_name})
        print(sample_name, "Percentage of mapped reads", round(df_AD[sample_name].sum()/total_read,2),"total read number",total_read/1000000 )
        return df_AD
    
    @staticmethod
    def _delete_files(file_path: str) -> None:
        """
        Delete a file if it exists.
        
        Parameters:
        file_path (str): Path to the file to be deleted.
        """
        if os.path.exists(file_path):
            os.remove(file_path)
            #print(f"{file_path} deleted successfully")
        else:
            print(f"{file_path} does not exist")

    def _get_read1_read2(self,dir_folder):
        """
        Get read1 and read2 files from a directory.
        
        Parameters:
        dir_folder (str): Path to the directory containing fastq files.
        
        Returns:
        Dict[str, Dict[str, str]]: Dictionary mapping sample names to their corresponding read1 and read2 file names.
        """
        if not os.path.exists(dir_folder):
            print(f"The directory {dir_folder} does not exist")
            return
        
        fastq_files = [f for f in os.listdir(dir_folder) if f.endswith('.fastq.gz')]
        paired_files = {}
        for file in fastq_files:
            sample_name = file[:-16]
            if sample_name not in paired_files:
                paired_files[sample_name] = {'R1': None, 'R2': None}
            
            if file.endswith('_R1_001.fastq.gz'):
                paired_files[sample_name]['R1'] = file
            elif file.endswith('_R2_001.fastq.gz'):
                paired_files[sample_name]['R2'] = file
        return paired_files


    def read_files_for_sample(self, sample_name, files, dir_folder):
        """
        Read files for a specific sample and return sequences from read1 and read2.

        Parameters:
        sample_name (str): Name of the sample being processed.
        files (Dict[str, str]): Dictionary with 'R1' and 'R2' as keys and their respective file names as values.
        dir_folder (str): Directory where the files are located.

        Returns:
        Tuple[List[str], List[str]]: Tuple containing lists of sequences from read1 and read2 respectively.
        """
        file_read1 = os.path.join(dir_folder, files['R1']) if files['R1'] else None
        file_read2 = os.path.join(dir_folder, files['R2']) if files['R2'] else None

        if not file_read1 or not file_read2:
            print(f"Missing R1 or R2 file for sample {sample_name}")
            return None, None

        read1 = self._parse_fastq_to_list(file_read1)
        read2 = self._parse_fastq_to_list(file_read2)
        
        if not read1:
            print("read1 is empty")
            # Consider handling empty read1 here if necessary

        if not read2:
            print("read2 is empty")
            # Consider handling empty read2 here if necessary

        return read1, read2

    def loop_through_samples_and_build_count_tables(self, dir_folder):
        """
        Loop through all samples in a directory, build count tables, and save them as CSV files.

        Parameters:
        dir_folder (str): Directory where the fastq files are located.
        """
        paired_files = self._get_read1_read2(dir_folder)
        df = self.test_library

        for n, (sample_name, files) in enumerate(paired_files.items(), start=1):
            read1, read2 = self.read_files_for_sample(sample_name, files, dir_folder)

            if read1 is None or read2 is None:
                print(sample_name, "Read1 or Read2 is None")
                continue  # Skip to the next iteration if read1 or read2 is None

            #print(sample_name, "Read1 and Read2 Received")

            # Annotate and update DataFrame, then save it as CSV
            df = self._annotate_AD_barcodes(sample_name, read1, read2, df).sort_values(by="enhancer_id", ascending=False)
            df.to_csv(f"df_MPRA3_{n}.csv",index=False)

            # Delete previous CSV file to save space
            if n > 1:
                self._delete_files(f"df_MPRA3_{n-1}.csv")

if __name__ == "__main__":
    mpra = MPRA(test_library= 'MPRA3_Contributor_20231108_unique_GeneName_BarcodeEnhancerPair.csv')
    mpra.loop_through_samples_and_build_count_tables(dir_folder = './mpra_fastq')

