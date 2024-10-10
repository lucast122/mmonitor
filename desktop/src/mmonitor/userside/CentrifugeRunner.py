import gzip
import logging
import multiprocessing
import os
import subprocess
from build_mmonitor_pyinstaller import ROOT
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import concurrent
from threading import Lock
import gzip
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager
import shutil
import platform
import subprocess
import os
def process_file(file, temp_dir):
    is_gzipped = file.endswith('.gz')
    temp_file = os.path.join(temp_dir, os.path.basename(file) + ".tmp")

    open_func = gzip.open if is_gzipped else open

    with open_func(file, 'rt') as infile, open(temp_file, 'wt') as outfile:
        shutil.copyfileobj(infile, outfile)

    return temp_file

class CentrifugeRunner:

    def __init__(self):
        self.centrifuge_path = "centrifuge"  # Assuming centrifuge is in the system PATH
        self.logger = logging.getLogger('timestamp')
        self.cent_out = ""
        self.concat_file_name = ""
        self.pipeline_out = os.path.join(ROOT, "src", "resources", "pipeline_out")

    @staticmethod
    def concatenate_files_in_memory(file_paths):
        """
        Create a file-like object that virtually concatenates multiple files.
        """
        concatenated_content = StringIO()
        for file_path in file_paths:
            with open(file_path, 'r') as file:
                concatenated_content.write(file.read())
        concatenated_content.seek(0)
        return concatenated_content

    @staticmethod
    def concatenate_gzipped_files(input_files, output_file):
        print("Concatenating gzipped fastq files to prepare for analysis...")
        with open(output_file, 'wb') as outfile:
            for file in input_files:
                with open(file, 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)

    @staticmethod
    def concatenate_fastq_fast(input_files, output_file, is_gzipped=False):
        print("Concatenating fastq files to prepare for analysis...")
        open_func = gzip.open if is_gzipped else open
        with open_func(output_file, 'wb') as outfile:
            for file in input_files:
                with open_func(file, 'rb') as infile:
                    # Copy the contents of each input file to the output file
                    shutil.copyfileobj(infile, outfile)

    @staticmethod
    def concatenate_fastq_parallel(files, output_file):
        files_to_process = [file for file in files if "concatenated" not in file]
        print("files to process:")
        print(files_to_process)
        print(f"Concatenating {len(files_to_process)} files to {output_file} ...")

        temp_dir = "temp_fastq"
        os.makedirs(temp_dir, exist_ok=True)

        # Use process_file directly without lambda
        with ProcessPoolExecutor(max_workers=12) as executor:
            temp_files = list(executor.map(process_file, files_to_process, [temp_dir] * len(files_to_process)))

        with gzip.open(output_file, 'wt') as outfile:
            for temp_file in temp_files:
                with open(temp_file, 'rt') as infile:
                    shutil.copyfileobj(infile, outfile)
                os.remove(temp_file)

        shutil.rmtree(temp_dir)
        print(f"Concatenation completed: {output_file}")

    @staticmethod
    def concatenate_to_tempfile(file_paths):
        """
        Concatenate multiple files into a temporary file.
        Returns the path of the temporary file.
        """
        with tempfile.NamedTemporaryFile(delete=False, mode='w+') as temp_file:
            for file_path in file_paths:
                with open(file_path, 'r') as file:
                    temp_file.write(file.read())
            return temp_file.name


    @staticmethod
    def make_kraken_report(centrifuge_index_path, cent_out):
        cmd = f"centrifuge-kreport -x {centrifuge_index_path} {cent_out} > {cent_out.replace('cent_out', 'kraken_out')}"
        print(cmd)
        os.system(cmd)

    @staticmethod
    def process_line(line, centrifuge_index_path):
        columns = line.strip().split('\t')
        if len(columns) < 5:
            return
        col_path, output_file = columns[3], columns[4]
        cmd = ["centrifuge-kreport", "-x", centrifuge_index_path, col_path]
        with open(output_file, 'w') as output:
            subprocess.run(cmd, stdout=output, stderr=subprocess.DEVNULL)

    @staticmethod
    def make_kraken_report_from_tsv(file_path, centrifuge_index_path, max_workers=64):
        with open(file_path, 'r') as file, concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(CentrifugeRunner.process_line, line, centrifuge_index_path) for line in file]
            concurrent.futures.wait(futures)

    @staticmethod
    def run_centrifuge_multi_sample(centrifuge_tsv_path, database_path):
        cmd = f'centrifuge -x "{database_path}" --sample-sheet {centrifuge_tsv_path} -p {multiprocessing.cpu_count()}'
        os.system(cmd)

    @staticmethod
    def create_centrifuge_input_file(sample_names, concat_file_names, output_tsv_file):
        """
        Creates a TSV file for running centrifuge in multi sample mode.

        :param sample_names: List of sample names.
        :param concat_file_name: The concatenated file name for read-file1.
        :param output_tsv_file: The path to the output TSV file.
        """
        with open(output_tsv_file, 'w') as file:
            # no header
            # Write the data rows
            for idx, sample_name in enumerate(sample_names):
                classification_result_output = f"{ROOT}/src/resources/pipeline_out/{sample_name}_cent_out"
                report_file = f"{ROOT}/src/resources/pipeline_out/{sample_name}_kraken_out"
                file.write(
                    f"1\t{concat_file_names[idx]}\tread-file2\t{classification_result_output}\t{report_file}\n")

    @staticmethod
    def unpack_fastq_list(ls):
        """
        Gets a list of paths and outputs a comma seperated string containing the paths used as input for centrifuge
        """
        if len(ls) == 1:
            return f"{ls[0]}"
        elif len(ls) > 1:
            return ",".join(ls)

    def check_centrifuge(self):
        try:
            subprocess.call(['centrifuge', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except FileNotFoundError:
            self.logger.error(
                "Make sure that centrifuge is installed and on the sytem path. For more info visit http://www.ccb.jhu.edu/software/centrifuge/manual.shtml")

    def run_centrifuge(self, concat_file_name, sample_name, cent_db_path):
        os.makedirs(self.pipeline_out, exist_ok=True)
        cmd = f"centrifuge -x {cent_db_path} -U {concat_file_name} --report-file {self.pipeline_out}/{sample_name}_centrifuge_report.tsv -S {self.pipeline_out}/{sample_name}_centrifuge_out.txt"
        print(cmd)
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            print("Centrifuge analysis completed successfully")
        except subprocess.CalledProcessError as e:
            print(f"Error running Centrifuge: {e}")
            print(f"Centrifuge stderr: {e.stderr}")
            if "command not found" in e.stderr:
                print("Centrifuge is not installed or not in the system PATH")
        except FileNotFoundError:
            print("Centrifuge is not installed or not in the system PATH")

    def get_files_from_folder(self, folder_path):
        """
        Gets a path to a folder, checks if path contains sequencing files with specified endings and returns list
        containing paths to sequencing files.
        """
        files = []
        for file in os.listdir(folder_path):
            if file.endswith((".fastq", ".fq", ".fasta", ".fastq.gz")):
                files.append(os.path.join(folder_path, file))
        return files

    def run(self, input_file, output_dir, db_path, threads=1):
        output_file = os.path.join(output_dir, "centrifuge_output.tsv")
        report_file = os.path.join(output_dir, "centrifuge_report.tsv")

        cmd = [
            self.centrifuge_path,
            "-x", db_path,
            "-U", input_file,
            "--threads", str(threads),
            "-S", output_file,
            "--report-file", report_file
        ]

        try:
            subprocess.run(cmd, check=True)
            print(f"Centrifuge analysis completed successfully. Results saved in {output_dir}")
        except subprocess.CalledProcessError as e:
            print(f"Error running Centrifuge: {e}")
            raise
