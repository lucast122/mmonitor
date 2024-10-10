import argparse
import csv
import gzip
import json
import os
import sys
import numpy as np
import logging
import subprocess
from datetime import datetime
import tempfile
import shutil
import traceback

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

logger = logging.getLogger(__name__)

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
    def __init__(self):
        self.use_multiplexing = None
        self.multi_sample_input = None
        self.emu_runner = None  # Initialize to None
        self.centrifuge_runner = CentrifugeRunner()
        self.functional_runner = FunctionalRunner()
        self.db_config = {}
        self.pipeline_out = os.path.join(ROOT, "src", "resources", "pipeline_out")
        self.args = None
        self.django_db = None
        self.centrifuge_db = None
        self.db_path = os.path.join(ROOT, "src", "resources", "db_config.json")
        self.output_dir = None  # We'll set this in the initialize_from_args method
        self.config = {}  # Initialize config as an empty dictionary
        self.emu_db_path = os.path.join(ROOT, "src", "resources", "emu_db")

        # Add minimap2 from lib folder to PATH
        minimap2_path = os.path.join(ROOT, "lib", "minimap2")
        os.environ['PATH'] = f"{minimap2_path}:{os.environ['PATH']}"

        # Optionally, you can print the updated PATH for debugging
        print(f"Updated PATH: {os.environ['PATH']}")

    def initialize_from_args(self, args):
        self.args = args
        self.output_dir = os.path.join(self.pipeline_out, self.args.sample)
        os.makedirs(self.output_dir, exist_ok=True)
        if self.args.emu_db:
            self.emu_runner = EmuRunner(custom_db_path=self.args.emu_db)
        else:
            self.emu_runner = EmuRunner(custom_db_path=self.emu_db_path)
        try:
            self.django_db = DjangoDBInterface(self.args.config or self.db_path)
            # Ensure the django_db is logged in
            if not self.django_db.verify_credentials():
                print("Error: Not logged in or invalid credentials.")
                sys.exit(1)
        except (FileNotFoundError, ValueError) as e:
            print(f"Error initializing database interface: {e}")
            sys.exit(1)
        self.centrifuge_db = args.centrifuge_db or os.path.join(ROOT, "src", "resources", "centrifuge_db", "p_compressed")

    def valid_file(self, path):
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(f"{path} is not a valid file path")
        return path

    def valid_directory(self, path):
        if not os.path.isdir(path):
            raise argparse.ArgumentTypeError(f"{path} is not a valid directory path")
        return path

    @staticmethod
    def valid_date(date_str):
        try:
            return datetime.strptime(date_str, '%Y-%m-%d').date()
        except ValueError:
            raise argparse.ArgumentTypeError(f"Invalid date format: {date_str}. Use YYYY-MM-DD")



    @staticmethod
    def parse_arguments(args=None):
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
        group.add_argument('-m', '--multicsv', type=MMonitorCMD.valid_file,
                           help='Path to CSV containing information for multiple samples.')
        group.add_argument('-i', '--input', nargs='+', help='Input files for single sample processing')

        # Additional parameters
        parser.add_argument('-s', '--sample', type=str, help='Sample name.')
        parser.add_argument('-d', '--date', type=MMonitorCMD.valid_date, help='Sample date in YYYY-MM-DD format.')
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

        parser.add_argument('--emu-db', type=str, help='Path to custom Emu database')
        parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use for processing.')
        parser.add_argument('--centrifuge-db', type=str, help='Path to the Centrifuge database')

        parsed_args = parser.parse_args(args)

        # Check for required arguments
        if parsed_args.analysis in ['functional', 'assembly'] and not parsed_args.input:
            parser.error("the following arguments are required: -i/--input")

        return parsed_args

    import argparse
    import os

    
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
        
        sample_name = str(self.args.sample)
        project_name = str(self.args.project)
        subproject_name = str(self.args.subproject)
        sample_date = self.args.date.strftime('%Y-%m-%d')
        input_files = self.args.input  # This should be a list of input files

        print(f"Analyzing amplicon data for sample {sample_name}.")
        emu_success = self.emu_runner.run_emu(input_files, sample_name, self.args.minabundance)
        
        if not emu_success:
            print(f"Emu analysis failed for sample {sample_name}.")
            return

        print(f"Emu analysis completed for sample {sample_name}.")
        
        # Update the database after successful Emu analysis
        if self.update_database_16s(sample_name, project_name, subproject_name, sample_date):
            print(f"Database successfully updated with results for sample {sample_name}.")
        else:
            print(f"Failed to update database with results for sample {sample_name}.")

    def update_database_16s(self, sample_name, project_name, subproject_name, sample_date):
        emu_out_path = os.path.join(ROOT, "src", "resources", "pipeline_out", sample_name)
        if not os.path.exists(emu_out_path):
            print(f"No output directory found for sample {sample_name}. Skipping database update.")
            return False
        
        abundance_file = os.path.join(emu_out_path, f"{sample_name}_rel-abundance.tsv")
        if not os.path.exists(abundance_file):
            print(f"No abundance file found for sample {sample_name}. Skipping database update.")
            return False
        
        try:
            self.django_db.update_django_with_emu_out(emu_out_path, "species", sample_name, project_name, sample_date,
                                                      subproject_name, self.args.overwrite)
            print(f"Database successfully updated with Emu results for sample {sample_name}.")
            return True
        except Exception as e:
            print(f"Error updating database with Emu results for sample {sample_name}: {e}")
            return False

    def taxonomy_nanopore_wgs(self):
        cent_db_path = os.path.join(ROOT, 'src', 'resources', 'dec_22')
        today = datetime.now()
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

            # Use the get_files_from_folder method from CentrifugeRunner
            files = self.centrifuge_runner.get_files_from_folder(self.args.input[0])  # Assuming args.input is a list with at least one element
            
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

            self.centrifuge_runner.run_centrifuge(concat_file_name, sample_name, self.centrifuge_db)
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
            CentrifugeRunner.run_centrifuge_multi_sample(centrifuge_tsv_path, self.centrifuge_db)

            print(f"Make kraken report from centrifuge reports...")

            CentrifugeRunner.make_kraken_report_from_tsv(centrifuge_tsv_path, self.centrifuge_db)


            print(f"Adding all samples to database...")
            for idx, sample in enumerate(sample_names_to_process):
                add_sample_to_databases(sample, project_names[idx], subproject_names[idx], sample_dates[idx])

                # calculate QC statistics if qc argument is given by user
                if self.args.qc:
                    print(f"Adding statistics for sample: {sample}...")
                    # print(f"Loading files: {file_paths}")
                    print(f"concat files list {concat_files_list}")


                    self.add_statistics(concat_files_list[idx], sample_names_to_process[idx], project_names[idx],
                                        subproject_names[idx],
                                        sample_dates[idx])
            print(f"Removing concatenated files...")
            for concat_file in concat_files_list:
                os.remove(concat_file)

    def assembly_pipeline(self):
        output_dir = os.path.join(self.pipeline_out, self.args.sample)
        os.makedirs(output_dir, exist_ok=True)

        # Assembly with Flye
        assembly_out = self.functional_runner.run_flye(self.args.input, self.args.sample, output_dir, self.args.threads)

        # Correction with Medaka
        corrected_assembly = self.functional_runner.run_medaka(assembly_out, self.args.input[0], output_dir, self.args.threads)

        # Binning with MetaBAT2
        bins_dir = self.functional_runner.run_metabat2(corrected_assembly, self.args.input[0], output_dir, self.args.threads)

        print("Assembly pipeline completed successfully.")

    def functional_pipeline(self):
        self.assembly_pipeline()
        output_dir = os.path.join(self.pipeline_out, self.args.sample)
        bins_dir = os.path.join(output_dir, "metabat2_bins")

        # CheckM2
        checkm2_out = self.functional_runner.run_checkm2(bins_dir, output_dir, self.args.threads)

        # Annotation with Bakta
        bakta_out = self.functional_runner.run_bakta(bins_dir, output_dir, self.args.threads)

        # Taxonomy with GTDB-TK
        gtdbtk_out = self.functional_runner.run_gtdbtk(bins_dir, output_dir, self.args.threads)

        print("Functional analysis pipeline completed successfully.")

    def concatenate_files(self, input_files, output_file):
        try:
            with open(output_file, 'wb') as outfile:
                for filename in input_files:
                    if filename.endswith('.gz'):
                        with gzip.open(filename, 'rb') as infile:
                            shutil.copyfileobj(infile, outfile)
                    else:
                        with open(filename, 'rb') as infile:
                            shutil.copyfileobj(infile, outfile)
            logger.info(f"Concatenated files saved to {output_file}")
        except Exception as e:
            logger.error(f"Error concatenating files: {str(e)}")
            raise

    def run(self):
        logger.info(f"Starting run with analysis type: {self.args.analysis}")
        
        logger.info(f"Current PATH: {os.environ['PATH']}")

        if self.args.analysis == "taxonomy-16s":
            self.emu_runner = EmuRunner(self.emu_runner.custom_db_path)
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                self.emu_runner.concatenate_fastq_files(self.args.input, temp_file.name)
                success = self.emu_runner.run_emu(
                    input_file=temp_file.name,
                    output_dir=self.output_dir,
                    db_dir=self.emu_runner.custom_db_path,
                    threads=self.args.threads,
                    N=50,
                    K="500M",
                    minimap_type="map-ont"
                )
            os.unlink(temp_file.name)

            if not success:
                self.logger.error("Emu analysis failed")
                return

        elif self.args.analysis == "taxonomy-wgs":
            self.centrifuge_runner.run(self.args.input, self.output_dir, self.args.centrifuge_db, self.args.threads)
        elif self.args.analysis in ["assembly", "functional"]:
            self.functional_runner.run(self.args.input, self.output_dir, self.args.threads)
        else:
            logger.error(f"Unsupported analysis type: {self.args.analysis}")

        logger.info("Finished run")

    def run_single_sample(self):
        if self.args.analysis == "taxonomy-16s":
            self.taxonomy_nanopore_16s()
        elif self.args.analysis == "taxonomy-wgs":
            self.taxonomy_nanopore_wgs()
        elif self.args.analysis == "assembly":
            self.assembly_pipeline()
        elif self.args.analysis == "functional":
            self.functional_pipeline()
        elif self.args.analysis == "kegg":
            self.run_kegg_analysis()
        elif self.args.analysis == "mag-upload":
            self.upload_mag()

    def run_multi_sample(self):
        self.load_from_csv()
        for index, sample_name in enumerate(self.multi_sample_input["sample_names"]):
            self.args.sample = sample_name
            self.args.project = self.multi_sample_input["project_names"][index]
            self.args.subproject = self.multi_sample_input["subproject_names"][index]
            self.args.date = self.multi_sample_input["dates"][index]
            self.args.input = self.multi_sample_input["file_paths_lists"][index]
            self.run_single_sample()

    def run_emu(self, files, sample_name, min_abundance):
        print(f"Running Emu for sample {sample_name}")
        # ... (rest of the method)
        # Replace any sys.stdout.write() or sys.stderr.write() with print()
        print(f"Emu analysis complete for sample {sample_name}")

    # Apply similar changes to other methods that might be writing directly to stdout/stderr

    def run_pipeline(self, selected_steps, params):
        input_files = self.get_input_files()
        output_dir = self.get_output_directory()

        if "concatenate" in selected_steps:
            concatenated_file = self.concatenate_files(input_files, os.path.join(output_dir, "concatenated.fastq"))
        
        if "filter" in selected_steps:
            filtered_file = self.functional_runner.run_filtlong(concatenated_file, output_dir, params["min_length"], params["min_quality"])
        
        if "assembly" in selected_steps:
            assembly_file = self.functional_runner.run_metaflye(filtered_file, output_dir, params["threads"])
        
        if "binning" in selected_steps:
            bins_dir = self.functional_runner.run_metabat2(assembly_file, filtered_file, output_dir, params["threads"])
        
        if "taxonomy" in selected_steps:
            gtdb_out = self.functional_runner.run_gtdbtk(bins_dir, output_dir, params["threads"])
        
        if "quality" in selected_steps:
            checkm2_out = self.functional_runner.run_checkm2(bins_dir, output_dir, params["threads"])
        
        if "annotation" in selected_steps:
            bakta_out = self.functional_runner.run_bakta(bins_dir, output_dir, params["threads"])
        
        if "phylogeny" in selected_steps:
            tree_file = self.functional_runner.build_phylogenetic_tree(gtdb_out, output_dir, params["threads"])
        
        # TODO: Implement step 9 - Upload annotations to mmonitor web app

    def get_input_files(self):
        # Implement method to get input files from user
        pass

    def get_output_directory(self):
        # Implement method to get output directory from user
        pass

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
    cmd_runner = MMonitorCMD()
    args = cmd_runner.parse_arguments()
    cmd_runner.initialize_from_args(args)
    cmd_runner.run()