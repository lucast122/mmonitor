import sys
import gzip
import logging
import multiprocessing
import os
import subprocess
import time
from contextlib import redirect_stdout

import pandas as pd

from build_mmonitor_pyinstaller import ROOT
from lib import emu

class EmuRunner:

    def __init__(self, custom_db_path=None):
        self.concat_file_name = None
        self.logger = logging.getLogger('timestamp')
        self.emu_path = os.path.join(ROOT, "lib", "emu.py")
        self.check_emu()
        self.emu_out = ""

    @staticmethod
    def unpack_fastq_list(ls):
        """
        Gets a list of paths and outputs a comma seperated string containing the paths used as input for emu.py
        """
        if len(ls) == 1:
            return f"{ls[0]}"
        elif len(ls) > 1:
            return ",".join(ls)

    def check_emu(self):
        try:
            result = subprocess.run([sys.executable, self.emu_path, '-h'], capture_output=True, text=True, check=True)
            output = result.stdout.strip()
            if output.startswith("usage: emu.py [-h]"):
                print("Emu is installed and accessible.")
            else:
                print("Unexpected output from Emu. Please check your installation.")
                print(f"Output: {output}")
        except subprocess.CalledProcessError as e:
            print(f"Error running Emu: {e}")
            print(f"Error output: {e.stderr}")
        except FileNotFoundError:
            print("Emu not found. Please make sure Emu is installed and accessible.")

    def run_emu(self, sequence_list, sample_name, min_abundance):
        print(f"Running emu with min abundance of {min_abundance}")
        self.emu_out = f"{ROOT}/src/resources/pipeline_out/{sample_name}/"

        if not sequence_list:
            print(f"Error: No input files provided for sample {sample_name}")
            return False

        #remove concatenated files from sequence list to avoid concatenating twice
        sequence_list = [s for s in sequence_list if "concatenated" not in s]
        print(f"Input files: {sequence_list}")
        
        if not sequence_list:
            print(f"Error: No valid input files found for sample {sample_name}")
            return False

        concat_file_name = f"{os.path.dirname(sequence_list[0])}/{sample_name}_concatenated.fastq.gz"
        self.concat_file_name = concat_file_name
        print(f"concat_file_name: {concat_file_name}")
        if not os.path.exists(concat_file_name):
            self.concatenate_fastq_files(sequence_list, concat_file_name)

        if not os.path.exists(concat_file_name):
            print(f"Error: Failed to create concatenated file {concat_file_name}")
            return False

        emu_db = f"{ROOT}/src/resources/emu_db/"
        if not os.path.exists(emu_db):
            print(f"Error: Emu database not found at {emu_db}")
            return False

        try:
            df_taxonomy = pd.read_csv(os.path.join(emu_db, "taxonomy.tsv"), sep='\t',
                                      index_col='tax_id', dtype=str)
        except Exception as e:
            print(f"Error reading taxonomy file: {e}")
            return False

        db_species_tids = df_taxonomy.index
        print(f"Emu out: {self.emu_out}")
        if not os.path.exists(self.emu_out):
            os.makedirs(self.emu_out)

        out_file_base = self.emu_out
        sam_out = f"{out_file_base}/emu_alignments.sam"
        tsv_out = f"{out_file_base}/{sample_name}_rel-abundance"
        print(f"min abundance: {min_abundance}")

        cmd = [
            sys.executable,
            self.emu_path,
            "abundance",
            "-i", concat_file_name,
            "-d", emu_db,
            "-o", self.emu_out,
            "--min-abundance", str(min_abundance)
        ]

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(result.stdout)
            print("Emu analysis completed successfully")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Emu analysis failed with return code {e.returncode}")
            print(f"Error output: {e.stderr}")
            return False

    def get_files_from_folder(self, folder_path):
        """
        Gets a path to a folder, checks if path contains sequencing files with specified endings and returns list
        containing paths to sequencing files.
        """
        files = []
        found = False
        try:
            for file in os.listdir(folder_path):
                # print(file)
                if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fasta") or file.endswith(
                        ".fastq.gz"):
                    files.append(f"{folder_path}/{file}")
                    found = True
            if not found:
                self.logger.error(f"No sequencing files (.fastq, .fq found at {folder_path}")
            return files
        except FileNotFoundError:
            self.logger.error(f"Invalid folder path")

    def concatenate_fastq_files(self, input_files, output_file):
        with gzip.open(output_file, 'wt') as output:
            for input_file in input_files:
                is_gzipped = input_file.endswith(".gz")
                open_func = gzip.open if is_gzipped else open
                mode = 'rt' if is_gzipped else 'r'

                with open_func(input_file, mode) as input2:
                    for line in input2:
                        output.write(line)
