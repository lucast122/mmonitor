import gzip
import logging
import multiprocessing
import os
import subprocess

import pandas as pd

from build_mmonitor_pyinstaller import ROOT
from lib import emu


class EmuRunner:

    def __init__(self, custom_db_path=None):
        self.concat_file_name = None
        self.logger = logging.getLogger('timestamp')
        self.check_emu()
        self.emu_out = ""
        self.emu_path = os.path.join(ROOT, "lib", "emu.py")
        self.default_db_path = os.path.join(ROOT, "src", "resources", "emu_db")
        self.db_path = custom_db_path if custom_db_path else self.default_db_path

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
            subprocess.run(["emu", '-h'], stdout=open(os.devnull, 'w'),
                            stderr=subprocess.STDOUT)
        except FileNotFoundError:
            self.logger.error(
                "Make sure that emu is installed and on the system path. For more info visit https://gitlab.com/treangenlab/emu")

    def run_emu(self, input_files, sample_name, min_abundance=0.01):
        output_dir = os.path.join(ROOT, "src", "resources", "pipeline_out", sample_name)
        os.makedirs(output_dir, exist_ok=True)

        print(f"Running emu with min abundance of {min_abundance}")
        print(f"Emu out: {output_dir}")
        print(f"min abundance: {min_abundance}")
        print(f"Using database: {self.db_path}")

        cmd = [
            self.emu_path, "abundance",
            "--type", "map-ont",
            "--keep-counts",
            "--output-dir", output_dir,
            "--min-abundance", str(min_abundance),
            "--threads", "12",
            "--output-basename", f"{sample_name}_rel-abundance",
            "--db", self.db_path
        ] + input_files

        try:
            subprocess.run(cmd, check=True)
            print(f"Emu analysis completed for sample {sample_name}.")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error running Emu: {e}")
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
