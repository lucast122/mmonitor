import argparse
import csv
import gzip
import json
import os
import sys
import numpy as np
import logging

from build_mmonitor_pyinstaller import ROOT
from mmonitor.userside.FastqStatistics import FastqStatistics

from src.mmonitor.database.django_db_interface import DjangoDBInterface
from src.mmonitor.userside.CentrifugeRunner import CentrifugeRunner
from src.mmonitor.userside.FunctionalRunner import FunctionalRunner
from src.mmonitor.userside.EmuRunner import EmuRunner
from Bio import SeqIO
import gzip
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from datetime import date, datetime
import argparse
import csv
import gzip
import json
import os
import sys
import numpy as np

from datetime import date, datetime

class NumpyEncoder(json.JSONEncoder):
    """ Custom encoder for numpy data types """
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyEncoder, self).default(obj)


class MMonitorCMD:
    @staticmethod
    def get_files_from_folder(folder_path, recursive=False):
        """
        Gets a path to a folder, checks if path contains sequencing files with specified endings and returns list
        containing paths to sequencing files.
        """
        files = []

        if recursive:
            for dirpath, _, filenames in os.walk(folder_path):
                for file in filenames:
                    if file.endswith((".fastq", ".fq", ".fasta", ".fastq.gz")) and "concatenated" not in file:
                        print(file)
                        files.append(os.path.join(dirpath, file))
        else:
            for file in os.listdir(folder_path):
                full_path = os.path.join(folder_path, file)
                if os.path.isfile(full_path) and file.endswith(
                        (".fastq", ".fq", ".fasta", ".fastq.gz")) and "concatenated" not in file:
                    files.append(full_path)

        if not files:
            print(f"No sequencing files (.fastq, .fq, .fasta, .fastq.gz) found at {folder_path}")

        return files

    def __init__(self):
        self.use_multiplexing = None
        self.multi_sample_input = None
        self.emu_runner = EmuRunner()
        self.centrifuge_runner = CentrifugeRunner()
        self.functional_runner = FunctionalRunner()
        # self.args = self.parse_arguments()
        self.db_config = {}
        # self.django_db = DjangoDBInterface(self.args.config)
        self.pipeline_out = os.path.join(ROOT, "src", "resources", "pipeline_out")

    def parse_arguments(self):
        parser = argparse.ArgumentParser(description='MMonitor command line tool for various genomic analyses.')

        # Main analysis type
        parser.add_argument('-a', '--analysis', required=True, choices=['taxonomy-wgs', 'taxonomy-16s', 'assembly', 'functional', 'stats'],
                            help='Type of analysis to perform. Choices are taxonomy-wgs, taxonomy-16s, assembly, functional and stats.'
                                 'Functional will run the functional analysis pipeline including assembly, correction, binning, annotation and KEGG analysis while assembly will only run assembly, correction and binning.')

        # Configuration file
        parser.add_argument('-c', '--config', required=True, type=str,
                            help='Path to JSON config file. Ensure the file is accessible.')

        # Input options: Multi CSV or single input folder
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-m', '--multicsv', type=self.valid_file,
                           help='Path to CSV containing information for multiple samples.')
        group.add_argument('-i', '--input', type=self.valid_directory,
                           help='Path to folder containing sequencing data.')

        # Additional parameters
        parser.add_argument('-s', '--sample', type=str, help='Sample name.')
        parser.add_argument('-d', '--date', type=self.valid_date, help='Sample date in YYYY-MM-DD format.')
        parser.add_argument('-p', '--project', type=str, help='Project name.')
        parser.add_argument('-u', '--subproject', type=str, help='Subproject name.')
        parser.add_argument('-b', '--barcodes', action="store_true",
                            help='Use barcode column from CSV for multiplexing.')
        parser.add_argument("--overwrite", action="store_true", help="Overwrite existing records. Defaults to False.")

        # Quality control and update options
        parser.add_argument('-q', '--qc', action="store_true", help='Calculate QC statistics for input samples.')
        parser.add_argument('-x', '--update', action="store_true",
                            help='Update counts and abundances to the MMonitor DB.')

        # Abundance threshold
        parser.add_argument('-n', '--minabundance', type=float, default=0.01,
                            help='Minimal abundance threshold for 16s taxonomy. Default is 0.01 what means that all taxa'
                                 'below 1% abundance will not be uploaded to the MMonitor server.')

        # Verbose and logging level
        parser.add_argument('-v', '--verbose', action="store_true", help='Enable verbose output.')
        parser.add_argument('--loglevel', type=str, default='INFO', choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'],
                            help='Set the logging level.')

        return parser.parse_args()

    import argparse
    import os

    def valid_file(self, path):
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(f"{path} is not a valid file path")
        return path

    def valid_directory(self, path):
        if not os.path.isdir(path):
            raise argparse.ArgumentTypeError(f"{path} is not a valid directory path")
        return path

    def valid_date(self, date_str):
        try:
            return datetime.strptime(date_str, '%Y-%m-%d').date()
        except ValueError:
            raise argparse.ArgumentTypeError(f"{date_str} is not a valid date. Use YYYY-MM-DD format.")


    def concatenate_fastq_files(self, file_paths, output_file):
        """
        Concatenate all FASTQ files from the provided list into a single file.

        :param file_paths: List of paths to FASTQ files.
        :param output_file: Path to the output concatenated FASTQ file.
        """
        if not file_paths:
            print("No FASTQ files provided.")
            return

        # Use a set to track unique files to avoid duplication
        unique_files = set(file_paths)

        with gzip.open(output_file, 'wb') as outfile:
            for fastq_file in unique_files:
                try:
                    if fastq_file.endswith(".gz"):
                        with gzip.open(fastq_file, 'rb') as infile:
                            outfile.write(infile.read())
                    else:
                        with open(fastq_file, 'rb') as infile:
                            outfile.write(infile.read())
                except Exception as e:
                    print(f"Error processing file {fastq_file}: {e}")
                    continue

        print(f"Concatenated {len(unique_files)} files into {output_file}")

    def concatenate_files(self, files, sample_name):
        if not files:
            raise ValueError("The files list is empty.")

        file_extension = ".fastq.gz" if files[0].endswith(".gz") else ".fastq"

        # Ensure the first file has a valid directory path
        first_file_dir = os.path.dirname(files[0])
        if not first_file_dir:
            raise ValueError("The directory of the first file is invalid.")

        base_dir = first_file_dir
        concat_file_name = os.path.join(base_dir, f"{sample_name}_concatenated{file_extension}")

        # Debug: Print the directory and concatenated file path
        print(f"Base directory: {base_dir}")
        print(f"Concatenated file path: {concat_file_name}")

        # Ensure the directory is writable
        if not os.access(base_dir, os.W_OK):
            raise PermissionError(f"Directory {base_dir} is not writable.")

        # Check if the concatenated file already exists
        if not os.path.exists(concat_file_name):
            self.concatenate_fastq_files(files, concat_file_name)
        else:
            print(f"Concatenated file {concat_file_name} already exists. Skipping concatenation.")

        return concat_file_name
    def load_config(self):
        if os.path.exists(self.args.config):
            try:
                with open(self.args.config, "r") as f:
                    self.db_config = json.load(f)
                    print(f"DB config {self.args.config} loaded successfully.")
            except json.JSONDecodeError:
                print("Couldn't load DB config. Please make sure you provided the correct path to the file.")
        else:
            print(f"Config path doesn't exist")

    def add_statistics(self, fastq_file, sample_name, project_name, subproject_name, sample_date, multi=True):
        fastq_stats = FastqStatistics(fastq_file, multi=multi)
        # Calculate statistics
        fastq_stats.quality_statistics()
        fastq_stats.read_lengths_statistics()
        quality_vs_lengths_data = fastq_stats.qualities_vs_lengths()
        gc_contents = fastq_stats.gc_content_per_sequence()
        # qual_dist = fastq_stats.quality_score_distribution()
        # q20_q30 = fastq_stats.q20_q30_scores()



        data = {
            'sample_name': sample_name,
            'project_id': project_name,
            'subproject_id': subproject_name,
            'date': sample_date,
            'mean_gc_content': float(fastq_stats.gc_content()),  # Ensure float
            'mean_read_length': float(np.mean(fastq_stats.lengths)),  # Convert with float()
            'median_read_length': float(np.median(fastq_stats.lengths)),  # Convert with float()
            'mean_quality_score': float(np.mean([np.mean(q) for q in fastq_stats.qualities])),  # Ensure float
            'read_lengths': json.dumps(quality_vs_lengths_data['read_lengths'], cls=NumpyEncoder),
            # Use custom encoder if needed
            'avg_qualities': json.dumps(quality_vs_lengths_data['avg_qualities'], cls=NumpyEncoder),
            # Use custom encoder if needed
            'number_of_reads': int(fastq_stats.number_of_reads()),  # Ensure int
            'total_bases_sequenced': int(fastq_stats.total_bases_sequenced()),  # Ensure int
            'gc_contents_per_sequence': json.dumps(gc_contents, cls=NumpyEncoder)

        }

        self.django_db.send_sequencing_statistics(data)

    # method to check if a sample is already in the database, if user does not want to overwrite results
    # returns true if sample_name is in db and false if not
    def check_sample_in_db(self, sample_id):
        samples_in_db = self.django_db.get_unique_sample_ids()
        if samples_in_db is not None:
            return sample_id in samples_in_db
        else:
            return []

    def update_only_statistics(self):

        if not self.args.multicsv:
            sample_name = str(self.args.sample)
            project_name = str(self.args.project)
            subproject_name = str(self.args.subproject)
            sample_date = self.args.date.strftime('%Y-%m-%d')  # Convert date to string format
            files = self.args.input
            self.add_statistics(files, sample_name, project_name, subproject_name,
                                sample_date, multi=True)

        else:
            self.load_from_csv()
            print("Processing multiple samples")
            for index, file_path_list in enumerate(self.multi_sample_input["file_paths_lists"]):
                files = file_path_list
                sample_name = self.multi_sample_input["sample_names"][index]

                print(f"processing sample {sample_name}")
                project_name = self.multi_sample_input["project_names"][index]
                subproject_name = self.multi_sample_input["subproject_names"][index]
                sample_date = self.multi_sample_input["dates"][index]

                self.add_statistics(files, sample_name, project_name, subproject_name,
                                    sample_date, multi=True)

    def taxonomy_nanopore_16s(self):
        global sample_name, project_name, subproject_name, sample_date


        def add_sample_to_databases(sample_name, project_name, subproject_name, sample_date):
            emu_out_path = f"{ROOT}/src/resources/pipeline_out/{sample_name}/"
            self.django_db.update_django_with_emu_out(emu_out_path, "species", sample_name, project_name, sample_date,
                                                      subproject_name, self.args.overwrite)
            print(self.args.overwrite)


        if not os.path.exists(os.path.join(ROOT, "src", "resources", "emu_db", "taxonomy.tsv")):
            print("emu db not found")

        if not self.args.multicsv:
            sample_name = str(self.args.sample)
            print(f"Analyzing amplicon data for sample {sample_name}.")
            # when a sample is already in the database and user does not want to overwrite quit now
            if not self.args.overwrite:
                if self.check_sample_in_db(sample_name):
                    return
            project_name = str(self.args.project)
            subproject_name = str(self.args.subproject)
            sample_date = self.args.date.strftime('%Y-%m-%d')  # Convert date to string format
            files = self.args.input
            if self.args.update:
                print("Update parameter specified. Will only update results from file.")
                add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
                return

            self.emu_runner.run_emu(files, sample_name, self.args.minabundance)

            add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
            if self.args.qc:
                self.add_statistics(self.emu_runner.concat_file_name, sample_name, project_name, subproject_name,
                                    sample_date)
                print("adding statistics")

        else:
            self.load_from_csv()
            print("Processing multiple samples")
            for index, file_path_list in enumerate(self.multi_sample_input["file_paths_lists"]):
                files = file_path_list
                sample_name = self.multi_sample_input["sample_names"][index]
                # when a sample is already in the database and user does not want to overwrite quit now
                if not self.args.overwrite:
                    if self.check_sample_in_db(sample_name):
                        print(
                            f"Sample {sample_name} already in DB and overwrite not specified, continue with next sample...")
                        continue
                project_name = self.multi_sample_input["project_names"][index]
                subproject_name = self.multi_sample_input["subproject_names"][index]
                sample_date = self.multi_sample_input["dates"][index]
                if self.args.update:
                    print("Update parameter specified. Will only update results from file.")
                    add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
                    continue

                self.emu_runner.run_emu(files, sample_name, self.args.minabundance)
                add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
                print(f"Finished processing sample {index+1} of {len(file_path_list)}")
                if self.args.qc:
                    self.add_statistics(self.emu_runner.concat_file_name, sample_name, project_name, subproject_name,
                                        sample_date)
                    print("adding statistics")

        # calculate QC statistics if qc argument is given by user

        # emu_out_path = f"{ROOT}/src/resources/pipeline_out/subset/"

        # if self.db is not None:
        #     self.db.update_table_with_emu_out(emu_out_path, "species", sample_name, "project", self.sample_date)
        #
        # self.db_mysql.update_django_with_emu_out(emu_out_path, "species", sample_name, project_name, sample_date,
        #                                          subproject_name)

        print("Analysis complete. You can start monitoring now.")

    def load_from_csv(self):
        file_path = self.args.multicsv

        if not file_path:
            print("CSV file not provided")
            return

        # Check if CSV is empty
        if os.path.getsize(file_path) == 0:
            print("Selected CSV file is empty.")
            return

        error_messages = []  # list to accumulate error messages

        # Prepare dictionary to hold multiple sample data
        self.multi_sample_input = {
            "file_paths_lists": [],
            "sample_names": [],
            "dates": [],
            "project_names": [],
            "subproject_names": []
        }

        with open(file_path, 'r') as file:
            reader = csv.DictReader(file)

            for row in reader:
                # Check if provided path exists
                if not os.path.exists(row["sample folder"].strip()):
                    error_message = f"Invalid path from CSV: {row['sample folder'].strip()}"
                    print(error_message)
                    error_messages.append(error_message)
                    continue

                # Look for the fastq_pass folder in the provided path and its child directories
                folder_path = None
                for root, dirs, files in os.walk(row["sample folder"].strip()):
                    if "fastq_pass" in dirs:
                        folder_path = os.path.join(root, "fastq_pass")
                        break
                    else:
                        folder_path = os.path.join(root)

                if not folder_path:
                    error_message = f"'fastq_pass' directory not found for path: {row['sample folder'].strip()}"
                    print(error_message)
                    error_messages.append(error_message)
                    continue

                # If multiplexing is selected, navigate further to the barcode_x folder
                if self.args.barcodes:
                    barcode_id_string = str(row['Barcode ID'])

                    if len(barcode_id_string) == 1:
                        barcode_id_string = f"0{barcode_id_string}"

                    barcode_folder = f"barcode{barcode_id_string}"
                    folder_path = os.path.join(folder_path, barcode_folder)

                    if not os.path.exists(folder_path):
                        error_message = f"Barcode folder '{barcode_folder}' not found."
                        print(error_message)
                        error_messages.append(error_message)
                        continue

                files = self.get_files_from_folder(folder_path)

                # Extract attributes from CSV and check for errors
                required_columns = ["sample_name", "date", "project_name", "subproject_name"]
                for col in required_columns:
                    if not row[col]:
                        error_message = f"Missing {col} in CSV for path: {row['sample folder']}"
                        print(error_message)
                        error_messages.append(error_message)

                # Store extracted CSV information to multi_sample_input
                self.multi_sample_input["file_paths_lists"].append(files)
                self.multi_sample_input["sample_names"].append(row["sample_name"])
                self.multi_sample_input["dates"].append(row["date"])
                self.multi_sample_input["project_names"].append(row["project_name"])
                self.multi_sample_input["subproject_names"].append(row["subproject_name"])

            # Provide feedback on the number of samples to be processed
            num_samples_to_process = len(self.multi_sample_input["file_paths_lists"])
            if num_samples_to_process <= 5:
                print(f"Processing {num_samples_to_process} samples.")
            else:
                print(f"Processing {num_samples_to_process} samples. This may take a while.")

        # Handle error messages from loading the CSV

    def taxonomy_nanopore_wgs(self):
        cent_db_path = os.path.join(ROOT, 'src', 'resources', 'dec_22')
        today = datetime.now()
        # Format the date as "year month day"
        default_date = today.strftime("%Y-%m-%d")

        def add_sample_to_databases(sample_name, project_name, subproject_name, sample_date):
            kraken_out_path = f"{ROOT}/src/resources/pipeline_out/{sample_name}_kraken_out"
            self.django_db.send_nanopore_record_centrifuge(kraken_out_path, sample_name, project_name, subproject_name,
                                                           sample_date, self.args.overwrite)

        if not os.path.exists(os.path.join(ROOT, "src", "resources", "dec_22.1.cf")):
            print("centrifuge db not found")

        if not self.args.multicsv:
            sample_name = str(self.args.sample)
            # when a sample is already in the database and user does not want to overwrite quit now
            if not self.args.overwrite:
                if self.check_sample_in_db(sample_name):
                    print("Sample is already in DB use --overwrite to overwrite it...")
                    return
            project_name = str(self.args.project)
            subproject_name = str(self.args.subproject)
            # sample_date = self.args.date.strftime('%Y-%m-%d')  # Convert date to string format
            sample_date = self.args.date  # Convert date to string format
            if sample_date is None:
                sample_date = default_date

            files = self.get_files_from_folder(self.args.input)
            if self.args.update:
                print("Update parameter specified. Will only update results from file.")
                add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
                return
            # concat_file_name = f"{os.path.dirname(files[0])}/{sample_name}_concatenated.fastq.gz"
            if files[0].endswith(".gz"):
                concat_file_name = f"{os.path.dirname(files[0])}/{sample_name}_concatenated.fastq.gz"
                CentrifugeRunner.concatenate_gzipped_files(files,concat_file_name)
            else:
                concat_file_name = f"{os.path.dirname(files[0])}/{sample_name}_concatenated.fastq"
                CentrifugeRunner.concatenate_fastq_fast(files, concat_file_name, False)

            self.centrifuge_runner.run_centrifuge(concat_file_name, sample_name, cent_db_path)
            add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
            if self.args.qc:
                self.add_statistics(self.centrifuge_runner.concat_file_name, sample_name, project_name, subproject_name,
                                    sample_date)
                print("adding statistics")
            os.remove(concat_file_name)

        else:
            self.load_from_csv()
            print("Processing multiple samples")
            concat_files_list = []
            all_file_paths = []
            sample_names_to_process = []
            project_names = []
            subproject_names = []
            sample_dates = []
            for index, file_path_list in enumerate(self.multi_sample_input["file_paths_lists"]):
                files = file_path_list
                all_file_paths.append(files)
                sample_name = self.multi_sample_input["sample_names"][index]
                print(f"Analyzing amplicon data for sample {sample_name}.")
                # when a sample is already in the database and user does not want to overwrite quit now
                if not self.args.overwrite:
                    if self.check_sample_in_db(sample_name):
                        print(
                            f"Sample {sample_name} already in DB and overwrite not specified, continue with next sample...")
                        continue

                sample_names_to_process.append(sample_name)
                if files[index].endswith(".gz"):
                    concat_file_name = f"{os.path.dirname(files[index])}/{sample_name}_concatenated.fastq.gz"
                else:
                    concat_file_name = f"{os.path.dirname(files[index])}/{sample_name}_concatenated.fastq"
                concat_files_list.append(concat_file_name)


                project_name = self.multi_sample_input["project_names"][index]
                subproject_name = self.multi_sample_input["subproject_names"][index]
                sample_date = self.multi_sample_input["dates"][index]

                project_names.append(project_name)
                subproject_names.append(subproject_name)
                sample_dates.append(sample_date)

            for idx, files in enumerate(all_file_paths):
                sample_name = self.multi_sample_input["sample_names"][idx]
                print(f"Analyzing amplicon data for sample {sample_name}.")
                total_files = len(all_file_paths)
                print(f"Concatenating fastq files... ({idx + 1}/{total_files})")
                if files[idx].endswith(".gz"):
                    concat_file_name = f"{os.path.dirname(files[idx])}/{sample_name}_concatenated.fastq.gz"
                    CentrifugeRunner.concatenate_gzipped_files(files,concat_file_name)
                else:
                    concat_file_name = f"{os.path.dirname(files[idx])}/{sample_name}_concatenated.fastq"
                    CentrifugeRunner.concatenate_fastq_fast(files, concat_file_name, False)
                concat_files_list.append(concat_file_name)
            centrifuge_tsv_path = os.path.join(ROOT, "src", "resources", "centrifuge.tsv")
            print(f"Creating centrifuge tsv...")
            CentrifugeRunner.create_centrifuge_input_file(self.multi_sample_input["sample_names"],
                                                                concat_files_list,
                                                                centrifuge_tsv_path)
            print(f"Running centrifuge for multiple samples from tsv {centrifuge_tsv_path}...")
            CentrifugeRunner.run_centrifuge_multi_sample(centrifuge_tsv_path, cent_db_path)

            print(f"Make kraken report from centrifuge reports...")

            CentrifugeRunner.make_kraken_report_from_tsv(centrifuge_tsv_path, cent_db_path)


            print(f"Adding all samples to database...")
            for idx, sample in enumerate(sample_names_to_process):
                add_sample_to_databases(sample, project_names[idx], subproject_names[idx], sample_dates[idx])

                # calculate QC statistics if qc argument is given by user
                if self.args.qc:
                    print(f"Adding statistics for sample: {sample}...")
                    print(f"Loading files: {file_paths}")
                    print(f"concat files list {concat_files_list}")


                    self.add_statistics(concat_files_list[idx], sample_names_to_process[idx], project_names[idx],
                                        subproject_names[idx],
                                        sample_dates[idx])
            print(f"Removing concatenated files...")
            for concat_file in concat_files_list:
                os.remove(concat_file)

    def assembly_pipeline(self, s_name, p_name, sp_name, s_date, fils):
        concat_file_name = self.concatenate_files(fils, s_name)
        print(concat_file_name)
        contig_file_path = self.functional_runner.run_flye(concat_file_name, s_name, True, True)
        out_path = os.path.join(self.pipeline_out, s_name)
        self.functional_runner.run_medaka_consensus(contig_file_path, concat_file_name, out_path)

        print(f" contig_file_path: {contig_file_path}")
        print(f" concat_file_path: {concat_file_name}")
        print(f" out_path: {out_path}")
        os.startfile(out_path)
        self.functional_runner.run_metabat2_pipeline(contig_file_path, concat_file_name, out_path)
        bins_dir = os.path.join(out_path, "metabat_bins")
        bakta_dir = os.path.join(out_path, "bakta_results")

        # self.functional_runner.run_checkm2(bins_dir, out_path)
        # self.functional_runner.run_bakta(contig_file_path, bakta_dir)
        self.functional_runner.run_gtdb_tk(bins_dir, out_path)

class OutputLogger:
    def __init__(self, log_file_path):
        self.log_file_path = log_file_path
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.DEBUG)
        self.logger.addHandler(self._create_console_handler())
        self.logger.addHandler(self._create_file_handler())

    def _create_console_handler(self):
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.DEBUG)
        console_handler.setFormatter(logging.Formatter('%(message)s'))
        return console_handler

    def _create_file_handler(self):
        file_handler = logging.FileHandler(self.log_file_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter('%(message)s'))
        return file_handler

    def start_logging(self):
        """
        Redirects stdout and stderr to the logging system.
        """
        sys.stdout = self._StreamToLogger(self.logger, logging.INFO)
        sys.stderr = self._StreamToLogger(self.logger, logging.ERROR)

    def stop_logging(self):
        """
        Restores the original stdout and stderr.
        """
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

    class _StreamToLogger:
        def __init__(self, logger, log_level):
            self.logger = logger
            self.log_level = log_level
            self.linebuf = ''

        def write(self, buf):
            for line in buf.rstrip().splitlines():
                self.logger.log(self.log_level, line.rstrip())

        def flush(self):
            pass


# Example usage

if __name__ == "__main__":
    def run_user_choice(sample_name, project_name, subproject_name, sample_date, files):
        if command_runner.args.analysis == "taxonomy-16s":
            command_runner.taxonomy_nanopore_16s()
        if command_runner.args.analysis == "taxonomy-wgs":
            command_runner.taxonomy_nanopore_wgs()
        if command_runner.args.analysis == "stats":
            command_runner.update_only_statistics()
        if command_runner.args.analysis == "assembly":
            command_runner.assembly_pipeline(sample_name, project_name, subproject_name, sample_date, files)

    command_runner = MMonitorCMD()
    print(command_runner.args)
    command_runner.load_config()
    command_runner.args = command_runner.parse_arguments()
    command_runner.django_db = DjangoDBInterface(command_runner.args.config)

    if not command_runner.args.multicsv:
        sample_name = str(command_runner.args.sample)
        if not command_runner.args.overwrite and self.check_sample_in_db(sample_name):
            print("Sample is already in DB use --overwrite to overwrite it...")
            pass
        project_name = str(command_runner.args.project)
        subproject_name = str(command_runner.args.subproject)
        sample_date = command_runner.args.date if command_runner.args.date else datetime.now().strftime('%Y-%m-%d')

        files = command_runner.get_files_from_folder(command_runner.args.input, False)
        concat_file_name = command_runner.concatenate_files(files, sample_name)

        print(f" contig_file_path: {contig_file_path}")
        print(f" concat_file_path: {concat_file_name}")
        print(f" out_path: {out_path}")
        run_user_choice(sample_name, project_name, subproject_name, sample_date, files)
    else:
        command_runner.load_from_csv()
        print("Processing multiple samples")
        concat_files_list = []
        all_file_paths = []
        sample_names_to_process = []
        project_names = []
        subproject_names = []
        sample_dates = []
        for index, file_path_list in enumerate(command_runner.multi_sample_input["file_paths_lists"]):
            files = file_path_list
            all_file_paths.append(files)
            sample_name = command_runner.multi_sample_input["sample_names"][index]
            print(f"Analyzing amplicon data for sample {sample_name}.")
            # when a sample is already in the database and user does not want to overwrite quit now
            if not command_runner.args.overwrite:
                if command_runner.check_sample_in_db(sample_name):
                    print(
                        f"Sample {sample_name} already in DB and overwrite not specified, continue with next sample...")
                    continue

            sample_names_to_process.append(sample_name)
            if files[index].endswith(".gz"):
                concat_file_name = f"{os.path.dirname(files[index])}/{sample_name}_concatenated.fastq.gz"
            else:
                concat_file_name = f"{os.path.dirname(files[index])}/{sample_name}_concatenated.fastq"
            concat_files_list.append(concat_file_name)

            project_name = command_runner.multi_sample_input["project_names"][index]
            subproject_name = command_runner.multi_sample_input["subproject_names"][index]
            sample_date = command_runner.multi_sample_input["dates"][index]

            project_names.append(project_name)
            subproject_names.append(subproject_name)
            sample_dates.append(sample_date)

        for idx, files in enumerate(all_file_paths):
            sample_name = command_runner.multi_sample_input["sample_names"][idx]
            total_files = len(all_file_paths)
            print(f"Concatenating fastq files... ({idx + 1}/{total_files})")
            if files[idx].endswith(".gz"):
                concat_file_name = f"{os.path.dirname(files[idx])}/{sample_name}_concatenated.fastq.gz"
                CentrifugeRunner.concatenate_gzipped_files(files, concat_file_name)
            else:
                concat_file_name = f"{os.path.dirname(files[idx])}/{sample_name}_concatenated.fastq"
                CentrifugeRunner.concatenate_fastq_fast(files, concat_file_name, False)
            concat_files_list.append(concat_file_name)
        for idx, sample in enumerate(sample_names_to_process):
            run_user_choice(command_runner.multi_sample_input["sample_names"][idx],
                            command_runner.multi_sample_input["project_names"][idx],
                            command_runner.multi_sample_input["subproject_names"][idx],
                            command_runner.multi_sample_input["dates"][idx],
                            files
                            )





    logger = OutputLogger(os.path.join(ROOT, 'src', 'resources', 'mmonitor_log.txt'))
    # logger.start_logging()

    # Your script's operations here. All stdout will be written to 'output_log.txt'.


    # logger.stop_logging()

    # logger.stop_logging()
