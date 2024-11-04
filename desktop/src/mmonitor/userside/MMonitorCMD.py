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
import keyring

from build_mmonitor_pyinstaller import ROOT
from mmonitor.userside.FastqStatistics import FastqStatistics
# from mmonitor.userside.CentrifugeRunner import CentrifugeRunner  # Update import
from src.mmonitor.database.django_db_interface import DjangoDBInterface
from src.mmonitor.userside.FunctionalRunner import FunctionalRunner
from src.mmonitor.userside.CentrifugerRunner import CentrifugerRunner
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
        self.emu_runner = None
        self.centrifuger_runner = CentrifugerRunner()  
        self.functional_runner = FunctionalRunner()
        self.db_config = {}
        self.pipeline_out = os.path.join(ROOT, "src", "resources", "pipeline_out")
        self.args = None
        self.django_db = None
        self.centrifuger_db = None
        self.db_path = os.path.join(ROOT, "src", "resources", "pipeline_config.json")
        self.output_dir = None
        self.config = {}
        self.emu_db_path = os.path.join(ROOT, "src", "resources", "emu_db")
        self.offline_mode = False
        self.logged_in = False
        self.current_user = None

        # Add minimap2 from lib folder to PATH
        minimap2_path = os.path.join(ROOT, "lib", "minimap2")
        os.environ['PATH'] = f"{minimap2_path}:{os.environ['PATH']}"
        print(f"Updated PATH: {os.environ['PATH']}")

        # Use separate config files
        self.db_config_path = os.path.join(ROOT, "src", "resources", "db_config.json")
        self.pipeline_config_path = os.path.join(ROOT, "src", "resources", "pipeline_config.json")
        
        # Load database config for authentication with error handling
        try:
            if os.path.exists(self.db_config_path):
                with open(self.db_config_path, 'r') as f:
                    content = f.read().strip()  # Remove any whitespace
                    try:
                        self.db_config = json.loads(content)
                    except json.JSONDecodeError as e:
                        print(f"Error parsing db_config.json: {e}")
                        print(f"Content causing error: {content}")
                        # Initialize with empty config
                        self.db_config = {}
        except Exception as e:
            print(f"Error loading db_config: {e}")
            self.db_config = {}
        
        # Initialize Django DB interface with proper config
        self.django_db = DjangoDBInterface(self.db_config_path)
        
        # Pass authentication info to Django DB interface
        if self.db_config:
            self.django_db.username = self.db_config.get('user')
            self.django_db.password = keyring.get_password("MMonitor", self.db_config.get('user'))

    def initialize_from_args(self, args):
        """Initialize with proper offline mode handling"""
        self.args = args
        self.output_dir = os.path.join(self.pipeline_out, self.args.sample)
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Load config file first
        try:
            with open(args.config, 'r') as f:
                self.config = json.load(f)
                print(f"Loaded config: {self.config}")
        except Exception as e:
            print(f"Error loading config file: {e}")
            self.config = {}
        
        # Initialize EMU runner with proper database path
        if self.args.emu_db:
            self.emu_runner = EmuRunner(custom_db_path=self.args.emu_db)
        else:
            self.emu_runner = EmuRunner(custom_db_path=self.emu_db_path)
        
        # Skip database initialization in offline mode
        if not self.offline_mode:
            try:
                self.django_db = DjangoDBInterface(self.args.config or self.db_path)
            except (FileNotFoundError, ValueError) as e:
                print(f"Error initializing database interface: {e}")
                sys.exit(1)
        
        # Set centrifuger_db path with proper fallback chain
        if hasattr(self.args, 'centrifuger_db') and self.args.centrifuger_db:
            self.centrifuger_db = os.path.abspath(self.args.centrifuger_db)
        elif 'centrifuger_db' in self.config:
            self.centrifuger_db = os.path.abspath(self.config['centrifuger_db'])
        else:
            self.centrifuger_db = os.path.abspath(os.path.join(ROOT, "src", "resources", "centrifuger_db"))
        
        print(f"Initializing MMonitorCMD with args:")
        print(f"Centrifuger DB: {self.centrifuger_db}")
        print(f"EMU DB: {self.args.emu_db}")

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
                            help='Type of analysis to perform.')

        # Configuration file
        parser.add_argument('-c', '--config', required=True, type=str,
                            help='Path to JSON config file.')

        # Input options: Multi CSV or single input folder
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-m', '--multicsv', type=MMonitorCMD.valid_file,
                           help='Path to CSV containing information for multiple samples.')
        group.add_argument('-i', '--input', nargs='+', help='Input files for single sample processing')

        # Sample information
        parser.add_argument('-s', '--sample', type=str, help='Sample name.')
        parser.add_argument('-d', '--date', type=MMonitorCMD.valid_date, help='Sample date in YYYY-MM-DD format.')
        parser.add_argument('-p', '--project', type=str, help='Project name.')
        parser.add_argument('-u', '--subproject', type=str, help='Subproject name.')
        
        # Quality control parameters
        parser.add_argument('--min-length', type=int, default=1000,
                           help='Minimum read length to include in analysis.')
        parser.add_argument('--min-quality', type=float, default=10.0,
                           help='Minimum read quality score to include in analysis.')
        parser.add_argument('--min-abundance', type=float, default=0.01,
                           help='Minimum abundance threshold for taxonomic classification.')
        
        # Processing options
        parser.add_argument('-b', '--barcodes', action="store_true",
                            help='Use barcode column from CSV for multiplexing.')
        parser.add_argument("--overwrite", action="store_true", 
                            help="Overwrite existing records.")
        parser.add_argument('-q', '--qc', action="store_true", 
                            help='Calculate QC statistics for input samples.')
        parser.add_argument('-x', '--update', action="store_true",
                            help='Update counts and abundances to the MMonitor DB.')
        
        # Database paths
        parser.add_argument('--emu-db', type=str, help='Path to custom Emu database')
        parser.add_argument('--centrifuger-db', type=str, help='Path to the Centrifuger database')  # Updated name
        
        # Performance options
        parser.add_argument('-t', '--threads', type=int, default=1, 
                            help='Number of threads to use for processing.')
        
        # Logging options
        parser.add_argument('-v', '--verbose', action="store_true", 
                            help='Enable verbose output.')
        parser.add_argument('--loglevel', type=str, default='INFO', 
                            choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'],
                            help='Set the logging level.')

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
        """Run WGS taxonomy analysis using Centrifuger"""
        if not self.args.multicsv:
            sample_name = str(self.args.sample)
            if not self.args.overwrite:
                if self.check_sample_in_db(sample_name):
                    print("Sample is already in DB use --overwrite to overwrite it...")
                    return
            project_name = str(self.args.project)
            subproject_name = str(self.args.subproject)
            sample_date = self.args.date
            
            files = self.centrifuger_runner.get_files_from_folder(self.args.input[0])
            
            if self.args.update:
                print("Update parameter specified. Will only update results from file.")
                self.add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
                return

            # Concatenate input files
            if files[0].endswith(".gz"):
                concat_file_name = f"{os.path.dirname(files[0])}/{sample_name}_concatenated.fastq.gz"
                self.centrifuger_runner.concatenate_gzipped_files(files, concat_file_name)
            else:
                concat_file_name = f"{os.path.dirname(files[0])}/{sample_name}_concatenated.fastq"
                self.centrifuger_runner.concatenate_fastq_fast(files, concat_file_name, False)

            # Run Centrifuger analysis
            self.centrifuger_runner.run(concat_file_name, sample_name, self.centrifuger_db)
            self.add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
            
            if self.args.qc:
                self.add_statistics(concat_file_name, sample_name, project_name, subproject_name, sample_date)
                print("adding statistics")
            os.remove(concat_file_name)

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
        """Run the analysis with proper authentication handling"""
        logger.info(f"Starting run with analysis type: {self.args.analysis}")
        
        try:
            if self.args.analysis == "taxonomy-16s":
                # Create temporary concatenated file for EMU
                with tempfile.NamedTemporaryFile(delete=False, suffix='.fastq.gz') as temp_file:
                    self.centrifuger_runner.concatenate_gzipped_files(self.args.input, temp_file.name)
                    
                    # Run EMU analysis with correct parameters
                    success = self.emu_runner.run_emu(
                        input_file=temp_file.name,
                        output_dir=os.path.join(self.pipeline_out, self.args.sample),
                        db_dir=self.args.emu_db,
                        threads=self.args.threads,
                        N=50,  # Default EMU parameters
                        K="500000000",  # Default EMU parameters
                        
                        minimap_type="map-ont"
                    )
                    os.unlink(temp_file.name)
                    
                    if not success:
                        logger.error("EMU analysis failed")
                        return False

            elif self.args.analysis == "taxonomy-wgs":
                # Concatenate input files first
                with tempfile.NamedTemporaryFile(delete=False, suffix='.fastq.gz') as temp_file:
                    self.centrifuger_runner.concatenate_gzipped_files(self.args.input, temp_file.name)
                    
                    # Run Centrifuger with concatenated file
                    success = self.centrifuger_runner.run_centrifuger(
                        input_file=temp_file.name,
                        sample_name=self.args.sample,
                        db_path=self.centrifuger_db
                    )
                    os.unlink(temp_file.name)
                    
                    if not success:
                        logger.error("Centrifuger analysis failed")
                        return False

            logger.info("Analysis completed successfully")
            return True

        except Exception as e:
            logger.error(f"Error during analysis: {str(e)}")
            traceback.print_exc()
            return False

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
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Initialize the EMU runner
        emu_runner = EmuRunner(
            input_files=files,
            output_dir=self.output_dir,
            sample_name=sample_name,
            min_abundance=min_abundance,
            threads=self.args.threads,
            db_path=self.emu_db
        )
        
        # Run EMU analysis
        success = emu_runner.run()
        
        if success:
            print(f"Emu analysis complete for sample {sample_name}")
            return True
        else:
            print(f"Emu analysis failed for sample {sample_name}")
            return False

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








