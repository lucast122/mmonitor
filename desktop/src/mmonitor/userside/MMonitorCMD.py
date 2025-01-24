import os
import sys
import json
import gzip
import shutil
import traceback
import argparse
import logging
import numpy as np
import keyring
import requests
from datetime import datetime, date
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import getpass
import time
import tempfile
import multiprocessing

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

# Add the src directory to Python path
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from mmonitor.database.django_db_interface import DjangoDBInterface
from mmonitor.userside.FastqStatistics import FastqStatistics
from mmonitor.userside.AssemblyPipeline import AssemblyPipeline
from mmonitor.userside.CentrifugerRunner import CentrifugerRunner
from mmonitor.userside.EmuRunner import EmuRunner
from mmonitor.userside.FunctionalRunner import FunctionalRunner

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
        self.pipeline_out = os.path.join(src_dir, "resources", "pipeline_out")
        self.args = None
        self.django_db = None
        self.centrifuger_db = None
        self.db_path = os.path.join(src_dir, "resources", "pipeline_config.json")
        self.output_dir = None
        self.config = {}
        self.emu_db_path = os.path.join(src_dir, "resources", "emu_db")
        self.offline_mode = False
        self.logged_in = False
        self.current_user = None

        # Add minimap2 from lib folder to PATH
        minimap2_path = os.path.join(src_dir, "lib", "minimap2")
        os.environ['PATH'] = f"{minimap2_path}:{os.environ['PATH']}"
        print(f"Updated PATH: {os.environ['PATH']}")

        # Use separate config files
        self.db_config_path = os.path.join(src_dir, "resources", "db_config.json")
        self.pipeline_config_path = os.path.join(src_dir, "resources", "pipeline_config.json")
        
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
        
        # Initialize Django DB interface with proper offline mode
        self.django_db = DjangoDBInterface(self.args.config or self.db_path)
        self.django_db.offline_mode = self.offline_mode
        if self.offline_mode:
            self.django_db.set_offline_mode(True)
        
        # Set centrifuger_db path with proper fallback chain
        if hasattr(self.args, 'centrifuger_db') and self.args.centrifuger_db:
            self.centrifuger_db = os.path.abspath(self.args.centrifuger_db)
        elif 'centrifuger_db' in self.config:
            self.centrifuger_db = os.path.abspath(self.config['centrifuger_db'])
        else:
            # Use default path from config or fallback
            default_db = os.path.join(src_dir, "resources", "custom_centrifuger_db", 
                                    "ncbi_build_20241030_090603", "cfr_ncbi")
            self.centrifuger_db = os.path.abspath(default_db)
        
        print(f"Initializing MMonitorCMD with args:")
        print(f"Centrifuger DB: {self.centrifuger_db}")
        print(f"EMU DB: {self.args.emu_db}")
        print(f"Offline mode: {self.offline_mode}")

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
        parser.add_argument('-s', '--sample', type=str, required='--multicsv' not in sys.argv,
                          help='Sample name (required when not using --multicsv).')
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

    def concatenate_fastq_files(self, file_paths, output_file):
        """
        Concatenate FASTQ files in parallel using multiprocessing.
        
        :param file_paths: List of paths to FASTQ files.
        :param output_file: Path to the output concatenated FASTQ file.
        """
        if not file_paths:
            print("No FASTQ files provided.")
            return

        # Use a set to track unique files to avoid duplication
        unique_files = list(set(file_paths))
        print(f"Number of unique files to process: {len(unique_files)}")
        
        # Determine chunk size for parallel processing
        num_cores = min(os.cpu_count(), len(unique_files))
        chunk_size = max(1, len(unique_files) // num_cores)
        print(f"Using {num_cores} CPU cores with chunk size of {chunk_size}")
        
        # Create temporary directory for chunks
        temp_dir = os.path.join(os.path.dirname(output_file), ".temp_chunks")
        os.makedirs(temp_dir, exist_ok=True)

        def process_chunk(chunk_files, chunk_id):
            """Process a chunk of files in a separate process"""
            chunk_output = os.path.join(temp_dir, f"chunk_{chunk_id}.gz")
            print(f"Process {os.getpid()} starting chunk {chunk_id} with {len(chunk_files)} files")
            
            try:
                with gzip.open(chunk_output, 'wb') as outfile:
                    for idx, fastq_file in enumerate(chunk_files, 1):
                        try:
                            if fastq_file.endswith(".gz"):
                                with gzip.open(fastq_file, 'rb') as infile:
                                    shutil.copyfileobj(infile, outfile, length=1024*1024)
                            else:
                                with open(fastq_file, 'rb') as infile:
                                    shutil.copyfileobj(infile, outfile, length=1024*1024)
                            print(f"Process {os.getpid()} completed file {idx}/{len(chunk_files)} in chunk {chunk_id}")
                        except Exception as e:
                            print(f"Error processing file {fastq_file}: {e}")
                            continue
            except Exception as e:
                print(f"Error creating chunk {chunk_id}: {e}")

        # Split files into chunks and create processes
        processes = []
        for i in range(0, len(unique_files), chunk_size):
            chunk = unique_files[i:i + chunk_size]
            p = multiprocessing.Process(target=process_chunk, args=(chunk, i // chunk_size))
            processes.append(p)

        print(f"Created {len(processes)} processes for parallel processing")

        # Start all processes
        for p in processes:
            p.start()
            print(f"Started process {p.pid}")

        # Wait for all processes to complete
        for p in processes:
            p.join()
            print(f"Process {p.pid} completed")

        # Get list of chunk files
        chunk_files = sorted([os.path.join(temp_dir, f) for f in os.listdir(temp_dir) if f.startswith("chunk_")])
        print(f"Found {len(chunk_files)} chunk files")

        try:
            # Concatenate chunks into final output
            print("Concatenating chunks into final output...")
            with gzip.open(output_file, 'wb') as outfile:
                for chunk_file in chunk_files:
                    if os.path.exists(chunk_file):  # Check if chunk file exists
                        with gzip.open(chunk_file, 'rb') as infile:
                            shutil.copyfileobj(infile, outfile, length=1024*1024)
                        os.remove(chunk_file)  # Clean up chunk file
                    else:
                        print(f"Warning: Chunk file {chunk_file} not found")
            print("Finished concatenating chunks")
        except Exception as e:
            print(f"Error in concatenating chunks: {e}")
            print("Falling back to single-process mode")
            # Fall back to single-process concatenation
            with gzip.open(output_file, 'wb') as outfile:
                for fastq_file in unique_files:
                    try:
                        if fastq_file.endswith(".gz"):
                            with gzip.open(fastq_file, 'rb') as infile:
                                shutil.copyfileobj(infile, outfile, length=1024*1024)
                        else:
                            with open(fastq_file, 'rb') as infile:
                                shutil.copyfileobj(infile, outfile, length=1024*1024)
                    except Exception as e:
                        print(f"Error processing file {fastq_file}: {e}")
                        continue
        finally:
            # Clean up temporary directory
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                print(f"Error cleaning up temporary directory: {e}")

        print(f"Successfully concatenated {len(unique_files)} files into {output_file}")

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
        emu_out_path = os.path.join(self.pipeline_out, sample_name)
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

    def add_sample_to_databases(self, sample_name, project_name, subproject_name, sample_date):
        # Add sample to Django DB
        self.django_db.add_sample(sample_name, project_name, subproject_name, sample_date)
        
        # Add sample to Centrifuger DB
        self.centrifuger_runner.add_sample_to_centrifuger_db(sample_name, self.centrifuger_db)

    def run(self):
        """Run the analysis with proper authentication handling"""
        logger.info(f"Starting run with analysis type: {self.args.analysis}")
        
        try:
            # Authenticate first
            if not self.authenticate_user():
                logger.error("Authentication failed")
                return False
            
            # Create output directory for this sample
            sample_output_dir = os.path.join(self.pipeline_out, self.args.sample)
            os.makedirs(sample_output_dir, exist_ok=True)
            
            success = False
            if self.args.analysis == "taxonomy-16s":
                # EMU analysis implementation...
                pass
                
            elif self.args.analysis == "taxonomy-wgs":
                # Centrifuge analysis implementation...
                pass
                
            elif self.args.analysis == "assembly":
                try:
                    # Initialize FunctionalRunner
                    runner = FunctionalRunner()
                    
                    # Get input files
                    input_files = self.args.input
                    if isinstance(input_files, str):
                        input_files = [input_files]
                    
                    print(f"Starting assembly pipeline for {self.args.sample}")
                    
                    # Create temporary concatenated file
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.fastq.gz') as temp_file:
                        print(f"Concatenating input files to {temp_file.name}")
                        with gzip.open(temp_file.name, 'wb') as outfile:
                            for filename in input_files:
                                with gzip.open(filename, 'rb') as infile:
                                    shutil.copyfileobj(infile, outfile)
                        
                        # Run Flye assembly with concatenated file
                        print("Running Flye assembly...")
                        assembly_file = runner.run_flye(
                            input_files=[temp_file.name],  # Pass as list with single concatenated file
                            sample_name=self.args.sample,
                            output_dir=sample_output_dir,
                            threads=self.args.threads
                        )
                        
                        # Run Medaka for assembly correction
                        print("Running Medaka correction...")
                        corrected_assembly = runner.run_medaka(
                            assembly_file,
                            reads=temp_file.name,  # Use concatenated file
                            output_dir=sample_output_dir,
                            threads=self.args.threads
                        )
                        
                        # Run MetaBAT2 for binning
                        print("Running MetaBAT2 binning...")
                        bins_dir = runner.run_metabat2(
                            assembly=corrected_assembly,
                            reads=temp_file.name,  # Use concatenated file
                            output_dir=sample_output_dir,
                            threads=self.args.threads
                        )
                        
                        # Clean up temp file
                        os.unlink(temp_file.name)
                        
                        success = True
                        print(f"Assembly pipeline completed successfully for {self.args.sample}")
                        
                except Exception as e:
                    print(f"Error in assembly pipeline: {str(e)}")
                    traceback.print_exc()
                    success = False
                    
            elif self.args.analysis == "functional":
                try:
                    # Initialize FunctionalRunner
                    runner = FunctionalRunner()
                    
                    # Create temporary concatenated file
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.fastq.gz') as temp_file:
                        print(f"Concatenating input files to {temp_file.name}")
                        with gzip.open(temp_file.name, 'wb') as outfile:
                            for filename in self.args.input:
                                with gzip.open(filename, 'rb') as infile:
                                    shutil.copyfileobj(infile, outfile)
                        
                        # Run assembly steps with concatenated file
                        assembly_file = runner.run_flye(
                            [temp_file.name],
                            self.args.sample,
                            sample_output_dir,
                            self.args.threads
                        )
                        
                        corrected_assembly = runner.run_medaka(
                            assembly_file,
                            temp_file.name,
                            sample_output_dir,
                            self.args.threads
                        )
                        
                        bins_dir = runner.run_metabat2(
                            corrected_assembly,
                            temp_file.name,
                            sample_output_dir,
                            self.args.threads
                        )
                        
                        # Clean up temp file
                        os.unlink(temp_file.name)
                        
                        # Run additional functional analysis steps
                        print("Running CheckM2 quality assessment...")
                        checkm2_out = runner.run_checkm2(
                            bins_dir=bins_dir,
                            output_dir=sample_output_dir,
                            threads=self.args.threads
                        )
                        
                        print("Running GTDB-TK taxonomy classification...")
                        gtdbtk_out = runner.run_gtdbtk(
                            bins_dir=bins_dir,
                            output_dir=sample_output_dir,
                            threads=self.args.threads
                        )
                        
                        print("Running Bakta annotation...")
                        bakta_out = runner.run_bakta(
                            bins_dir=bins_dir,
                            output_dir=sample_output_dir
                        )
                        
                        success = True
                        print(f"Functional analysis pipeline completed successfully for {self.args.sample}")
                        
                except Exception as e:
                    print(f"Error in functional pipeline: {str(e)}")
                    traceback.print_exc()
                    success = False
            
            return success

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

    # method to check if a sample is already in the database, if user does not want to overwrite results
    # returns true if sample_name is in db and false if not
    def check_sample_in_db(self, sample_id):
        samples_in_db = self.django_db.get_unique_sample_ids()
        if samples_in_db is not None:
            return sample_id in samples_in_db
        else:
            return []

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

    def authenticate_user(self, max_attempts=3):
        """Authenticate user before running analysis with multiple attempts"""
        if self.offline_mode:
            # Start local server for offline mode
            try:
                server_path = "/Users/timo/Downloads/home/minion-computer/mmonitor_production/MMonitor/server"
                if not os.path.exists(server_path):
                    raise FileNotFoundError(f"Server directory not found at {server_path}")
                
                print(f"Starting Django server at {server_path}")
                os.chdir(server_path)
                
                # Start server process
                cmd = [sys.executable, "manage.py", "runserver", "127.0.0.1:8000"]
                self.server_process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    bufsize=1
                )
                
                # Wait for server to start
                max_server_attempts = 30
                attempt = 0
                while attempt < max_server_attempts:
                    try:
                        response = requests.get("http://127.0.0.1:8000", timeout=1)
                        if response.status_code == 200:
                            print("Server detected as running!")
                            break
                    except requests.exceptions.RequestException:
                        attempt += 1
                        time.sleep(1)
                
                # Set up offline mode credentials
                self.django_db.set_offline_mode(True)
                return True
                
            except Exception as e:
                print(f"Error starting offline server: {e}")
                return False
        else:
            # Online mode authentication with multiple attempts
            attempts = 0
            while attempts < max_attempts:
                try:
                    attempts += 1
                    remaining_attempts = max_attempts - attempts
                    
                    if attempts > 1:
                        print(f"\nAuthentication attempt {attempts}/{max_attempts}")
                        if remaining_attempts > 0:
                            print(f"{remaining_attempts} attempts remaining")
                    
                    if not self.db_config.get('user'):
                        print("No user found in config. Please log in.")
                        username = input("Username: ")
                        password = getpass.getpass("Password: ")
                    else:
                        username = self.db_config['user']
                        password = getpass.getpass("Password: ")
                    
                    success = self.django_db.login(
                        username=username,
                        password=password,
                        host=self.db_config.get('host', 'mmonitor.org'),
                        port=self.db_config.get('port', '443'),
                        remember=False
                    )
                    
                    if success:
                        print("Authentication successful")
                        return True
                    else:
                        if remaining_attempts > 0:
                            retry = input("Authentication failed. Would you like to try again? (y/n): ")
                            if retry.lower() != 'y':
                                print("Authentication cancelled by user")
                                return False
                        else:
                            print("Authentication failed. Maximum attempts reached.")
                            return False
                        
                except KeyboardInterrupt:
                    print("\nAuthentication cancelled by user")
                    return False
                except Exception as e:
                    print(f"Error during authentication: {e}")
                    if remaining_attempts > 0:
                        retry = input("Would you like to try again? (y/n): ")
                        if retry.lower() != 'y':
                            return False
                    else:
                        return False
            
            return False

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
