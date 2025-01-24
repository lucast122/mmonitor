import os
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
import gzip
import logging
import sys
from ..paths import SRC_DIR, LIB_DIR, RESOURCES_DIR

# Add local-build/emu to Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
emu_path = os.path.join(project_root, 'local-build', 'emu')
if emu_path not in sys.path:
    sys.path.insert(0, emu_path)

try:
    from emu import emu
except ImportError:
    print(f"Warning: Could not import emu module from {emu_path}")
    print("Make sure the emu module is installed in local-build/emu/")

class EmuRunner:

    def __init__(self, custom_db_path=None):
        self.concat_file_name = None
        self.logger = logging.getLogger('timestamp')
        self.emu_path = os.path.join(emu_path, "emu")
        self.check_emu()
        self.emu_out = ""
        self.custom_db_path = custom_db_path or os.environ.get('EMU_DATABASE_DIR')

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
        except subprocess.CalledProcessError:
            print("Error running Emu. Please check your installation.")
        except FileNotFoundError:
            print(f"Emu script not found at {self.emu_path}")

    def run_emu(self, input_file, output_dir, db_dir, threads, N=50, K="500M", minimap_type="map-ont"):
        """Run EMU analysis with proper parameter handling"""
        self.emu_out = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Get sample name from output directory
        sample_name = os.path.basename(output_dir)
        
        # Convert K to numeric format if given as string with M/G suffix
        if isinstance(K, str):
            if K.endswith('M'):
                K = str(int(K[:-1]) * 1000000)
            elif K.endswith('G'):
                K = str(int(K[:-1]) * 1000000000)
        
        cmd = [
            sys.executable,
            self.emu_path,
            "abundance",
            "--db", str(db_dir),
            "--threads", str(threads),
            "--output-dir", str(output_dir),
            "--type", minimap_type,
            "--min-abundance", "0.0001",
            "--keep-files",
            "--N", str(N),
            "--K", str(K),
            "--output-basename", sample_name,
            str(input_file)
        ]

        try:
            print(f"Running EMU analysis:")
            print(f"Input file: {input_file}")
            print(f"Output directory: {output_dir}")
            print(f"Database: {db_dir}")
            print(f"Threads: {threads}")
            print("Full EMU command:", ' '.join(cmd))
            
            result = subprocess.run(
                cmd, 
                check=True, 
                capture_output=True, 
                text=True
            )
            
            print("\nEMU stdout:")
            print(result.stdout)
            
            if result.stderr:
                print("\nEMU stderr:")
                print(result.stderr)
            
            print("EMU analysis completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"EMU analysis failed with return code {e.returncode}")
            print(f"Error output: {e.stderr}")
            return False
        except Exception as e:
            print(f"Unexpected error running EMU: {e}")
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
                if file.endswith((".fastq", ".fq", ".fasta", ".fastq.gz")):
                    files.append(os.path.join(folder_path, file))
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
