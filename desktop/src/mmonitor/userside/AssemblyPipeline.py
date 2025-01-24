import os
import sys
import json
import gzip
import shutil
import traceback
import multiprocessing
import subprocess
import logging
from pathlib import Path
import customtkinter as ctk
from .PipelineStateTracker import PipelineStateTracker, PipelineStep
from .ToolInstaller import ToolInstaller
from .MetaBatRunner import MetaBatRunner
from ..paths import SRC_DIR, LIB_DIR, RESOURCES_DIR

logger = logging.getLogger(__name__)

class AssemblyPipeline:
    def __init__(self, sample_name, output_dir, config=None, progress_callback=None, offline_mode=False, is_isolate=False):
        """Initialize the assembly pipeline
        
        Args:
            sample_name (str): Name of the sample
            output_dir (str): Directory for pipeline outputs
            config (dict, optional): Configuration dictionary. Defaults to None.
            progress_callback (callable, optional): Callback for progress updates. Defaults to None.
            offline_mode (bool, optional): Whether to run in offline mode. Defaults to False.
            is_isolate (bool, optional): Whether sample is an isolate. Defaults to False.
        """
        self.sample_name = sample_name
        self.output_dir = output_dir
        self.config = config or {}
        self.progress_callback = progress_callback
        self.offline_mode = offline_mode
        self.is_isolate = is_isolate
        
        # Initialize state tracker
        self.state_tracker = PipelineStateTracker(output_dir, sample_name)
        
        # Set default threads to max available
        self.default_threads = multiprocessing.cpu_count()
        
        # Set up tool paths
        self.lib_path = LIB_DIR
        self.flye_path = os.path.join(LIB_DIR, "Flye-2.9.5", "bin", "flye")
        self.tool_installer = ToolInstaller()
        self.minimap2_path = os.path.join(LIB_DIR, "minimap2", "minimap2")
        self.metabat_path = self._check_and_install_tool('metabat2')
        self.checkm2_path = self._check_and_install_tool('checkm2')
        self.bakta_path = self._check_and_install_tool('bakta')
        self.gtdbtk_path = self._check_and_install_tool('gtdbtk')
        
        # Initialize runners
        self.metabat_runner = MetaBatRunner()

    def is_step_completed(self, step):
        """Check if a pipeline step is completed"""
        return self.state_tracker.step_completed(step)

    def complete_step(self, step, output=None):
        """Mark a pipeline step as completed"""
        self.state_tracker.mark_step_complete(step)
        if output:
            self.state_tracker.set_step_output(step, output)

    def _check_and_install_tool(self, tool_name, required=False):
        """Check if a tool is available and optionally install it"""
        tool_path = self.tool_installer.check_tool(tool_name)
        if tool_path is None and required:
            raise RuntimeError(f"{tool_name} is required but not installed and user declined installation")
        return tool_path

    def _check_and_install_checkm2_db(self):
        """Check for CheckM2 DB and offer to download if missing"""
        try:
            # Run checkm2 database check
            result = subprocess.run(['checkm2', 'database', 'status'], capture_output=True, text=True)
            if "not found" in result.stdout or "not found" in result.stderr:
                print("CheckM2 database not found.")
                user_input = input("Would you like to download the CheckM2 database? [y/N] ")
                if user_input.lower() == 'y':
                    print("Downloading CheckM2 database...")
                    subprocess.run(['checkm2', 'database', 'download'], check=True)
                    return True
                else:
                    print("Skipping CheckM2 database download")
                    return False
            return True
        except Exception as e:
            print(f"Error checking/installing CheckM2 database: {e}")
            return False

    def _check_and_install_gtdbtk_db(self):
        """Check for GTDB-TK DB and offer to download if missing"""
        try:
            # Check GTDB_DATA_PATH environment variable
            gtdb_path = os.environ.get('GTDBTK_DATA_PATH')
            if not gtdb_path:
                gtdb_path = os.path.join(self.output_dir, "databases", "gtdbtk")
                os.environ['GTDBTK_DATA_PATH'] = gtdb_path

            # Check if database exists
            if not os.path.exists(gtdb_path) or not os.listdir(gtdb_path):
                print("GTDB-TK database not found.")
                user_input = input("Would you like to download the GTDB-TK database? (This may take a while) [y/N] ")
                if user_input.lower() == 'y':
                    print("Downloading GTDB-TK database...")
                    os.makedirs(gtdb_path, exist_ok=True)
                    subprocess.run(['gtdbtk', 'download', '--out_dir', gtdb_path], check=True)
                    return True
                else:
                    print("Skipping GTDB-TK database download")
                    return False
            return True
        except Exception as e:
            print(f"Error checking/installing GTDB-TK database: {e}")
            return False

    def _check_and_install_bakta_db(self):
        """Check for Bakta DB and offer to download if missing"""
        try:
            db_dir = os.path.join(self.output_dir, "databases", "bakta")
            db_path = os.path.join(db_dir, "db")
            
            if not os.path.exists(db_path):
                print("Bakta database not found.")
                user_input = input("Would you like to download the Bakta database? [y/N] ")
                if user_input.lower() == 'y':
                    print("Downloading Bakta database...")
                    os.makedirs(db_dir, exist_ok=True)
                    subprocess.run(['bakta', 'download', '--output', db_dir, '--type', 'light'], check=True)
                    return db_path
                else:
                    print("Skipping Bakta database download")
                    return None
            return db_path
        except Exception as e:
            print(f"Error checking/installing Bakta database: {e}")
            return None

    def _get_auth_token(self):
        """Get authentication token from login session or offline mode"""
        try:
            auth_url = f"{self.server_url}/users/get_user_id/"
            headers = {'Content-Type': 'application/json'}
            
            if self.offline_mode:
                auth_data = {
                    'username': 'offlinemode',
                    'password': 'offline123'
                }
            else:
                # Use provided credentials from config
                auth_data = {
                    'username': self.config.get('username'),
                    'password': self.config.get('password')
                }
                
            if not auth_data['username'] or not auth_data['password']:
                print("No credentials provided")
                return None
                
            print(f"Attempting authentication with {auth_data['username']} at {auth_url}")
            response = requests.post(
                auth_url,
                json=auth_data,  # Use json parameter to set Content-Type automatically
                verify=False,
                timeout=10
            )
            print(f"Response status code: {response.status_code}")
            print(f"Response content: {response.content}")
            
            if response.status_code == 200:
                data = response.json()
                if 'token' in data:
                    print("Successfully got authentication token")
                    return data['token']
                print(f"No token in response: {data}")
            else:
                print(f"Failed to get session token: {response.status_code} {response.text}")
        except Exception as e:
            print(f"Error getting authentication token: {str(e)}")
            traceback.print_exc()
        return None

    def _upload_results(self, assembly_path=None, bin_dir=None):
        """Upload results to server"""
        print("Progress: Uploading results...")
        try:
            if not self.db_interface:
                raise Exception("Database interface not initialized. Make sure set_main_window was called.")
            
            # Prepare MAG data
            mag_data = {}
            if assembly_path:
                mag_data['assembly'] = assembly_path
            if bin_dir:
                mag_data['bins'] = bin_dir
            
            # Upload MAG data
            success = self.db_interface.upload_mag(
                name=f"{self.sample_name}_assembly",
                taxonomy="Unknown",
                sample_name=self.sample_name,
                gff_file_path=f"{assembly_path}.gff" if assembly_path else None,
                fasta_file_path=assembly_path
            )
            
            if not success:
                raise Exception("Failed to upload results")
            
            print("Progress: Results uploaded successfully")
            return True
            
        except Exception as e:
            print(f"Error uploading results: {str(e)}")
            traceback.print_exc()
            raise

    def set_main_window(self, main_window):
        """Set reference to main window for accessing shared resources"""
        self.main_window = main_window
        self.db_interface = main_window.django_db

    def run_assembly(self, input_files, threads=12, assembly_mode="nano-raw", progress_callback=None):
        """Run the complete assembly pipeline with checkpointing"""
        
        self.progress_callback = progress_callback
        results = {}
        
        try:
            # Ensure required tools are available
            self._check_required_tools()
            
            # Create sample directory
            self.sample_dir = os.path.join(self.output_dir, self.sample_name)
            os.makedirs(self.sample_dir, exist_ok=True)
            
            # 1. Flye Assembly
            if not self.is_step_completed(PipelineStep.FLYE_ASSEMBLY):
                print("Running Flye assembly...")
                flye_output = os.path.join(self.sample_dir, "flye_out")
                # Handle multiple input files
                if len(input_files) > 1:
                    # Concatenate input files
                    concat_file = os.path.join(self.sample_dir, f"{self.sample_name}_concatenated.fastq")
                    with open(concat_file, 'wb') as outfile:
                        for input_file in input_files:
                            if input_file.endswith('.gz'):
                                with gzip.open(input_file, 'rb') as infile:
                                    outfile.write(infile.read())
                            else:
                                with open(input_file, 'rb') as infile:
                                    outfile.write(infile.read())
                    input_file = concat_file
                else:
                    input_file = input_files[0]
                
                assembly_file = self.run_flye(input_file, self.sample_name, flye_output, threads)
                if not assembly_file:
                    raise RuntimeError("Flye assembly failed")
                results["assembly"] = assembly_file
                self.state_tracker.mark_step_complete(PipelineStep.FLYE_ASSEMBLY, assembly_file)
            else:
                print("Using existing Flye assembly")
                results["assembly"] = self.state_tracker.get_step_output(PipelineStep.FLYE_ASSEMBLY)
            
            # Run Medaka polishing if available
            medaka_path = self._get_medaka_binary()
            if medaka_path:
                self.state_tracker.set_current_step(PipelineStep.MEDAKA_POLISH)
                if progress_callback:
                    progress_callback("Running Medaka polishing...")
                medaka_results = self._run_medaka(results["assembly"], input_files[0], os.path.join(self.sample_dir, "medaka"), threads)
                if medaka_results:
                    results.update(medaka_results)
                    self.state_tracker.mark_step_complete(PipelineStep.MEDAKA_POLISH, medaka_results)
            else:
                print("Skipping Medaka polishing - tool not available")

            # Run MetaBAT2 binning if available and not in isolate mode
            if not self.is_isolate:
                metabat_path = self.tool_installer.check_tool('metabat2')
                if metabat_path:
                    self.state_tracker.set_current_step(PipelineStep.METABAT_BINNING)
                    if progress_callback:
                        progress_callback("Running MetaBAT2 binning...")
                    binning_results = self._run_metabat(results["assembly"], threads)
                    if binning_results:
                        results.update(binning_results)
                        self.state_tracker.mark_step_complete(PipelineStep.METABAT_BINNING, binning_results)
                else:
                    print("Skipping MetaBAT2 binning - tool not available")

            # Run CheckM2 quality assessment if available and bins exist
            if "bins" in results:
                checkm2_path = self.tool_installer.check_tool('checkm2')
                if checkm2_path:
                    self.state_tracker.set_current_step(PipelineStep.CHECKM2_QUALITY)
                    if progress_callback:
                        progress_callback("Running CheckM2 quality assessment...")
                    quality_results = self._run_checkm2(results["bins"], threads)
                    if quality_results:
                        results.update(quality_results)
                        self.state_tracker.mark_step_complete(PipelineStep.CHECKM2_QUALITY, quality_results)
                else:
                    print("Skipping CheckM2 quality assessment - tool not available")

            # Run Bakta annotation if available and bins exist
            if "bins" in results:
                bakta_path = self.tool_installer.check_tool('bakta')
                if bakta_path:
                    self.state_tracker.set_current_step(PipelineStep.BAKTA_ANNOTATION)
                    if progress_callback:
                        progress_callback("Running Bakta annotation...")
                    annotation_results = self._run_bakta(results["bins"], threads)
                    if annotation_results:
                        results.update(annotation_results)
                        self.state_tracker.mark_step_complete(PipelineStep.BAKTA_ANNOTATION, annotation_results)
                else:
                    print("Skipping Bakta annotation - tool not available")

            # Run GTDB-Tk taxonomy classification if available and bins exist
            if "bins" in results:
                gtdbtk_path = self.tool_installer.check_tool('gtdbtk')
                if gtdbtk_path:
                    self.state_tracker.set_current_step(PipelineStep.GTDBTK_TAXONOMY)
                    if progress_callback:
                        progress_callback("Running GTDB-Tk classification...")
                    taxonomy_results = self._run_gtdbtk(results["bins"], threads)
                    if taxonomy_results:
                        results.update(taxonomy_results)
                        self.state_tracker.mark_step_complete(PipelineStep.GTDBTK_TAXONOMY, taxonomy_results)
                else:
                    print("Skipping GTDB-Tk classification - tool not available")

            # Upload results
            self.state_tracker.set_current_step(PipelineStep.UPLOAD_TO_SERVER)
            if progress_callback:
                progress_callback("Uploading results...")
            try:
                if not self._upload_results(results["assembly"], results.get("bins")):
                    raise RuntimeError("Failed to upload results to server")
                self.state_tracker.mark_step_complete(PipelineStep.UPLOAD_TO_SERVER, {"uploaded_files": list(results.keys())})
            except Exception as e:
                self.state_tracker.mark_step_failed(PipelineStep.UPLOAD_TO_SERVER, str(e))
                raise

            return results

        except Exception as e:
            print(f"Error in assembly pipeline: {str(e)}")
            if progress_callback:
                progress_callback(f"Error: {str(e)}")
            # Mark pipeline as failed in state tracker
            if isinstance(self.state_tracker.current_step, PipelineStep):
                self.state_tracker.mark_step_failed(self.state_tracker.current_step, str(e))
            raise

    def stop(self):
        """Stop the pipeline execution"""
        self._stop_requested = True

    def _get_medaka_binary(self):
        """Get the path to Medaka binary, trying bundled version first then system"""
        # First try bundled version
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
        if sys.platform == "darwin":  # macOS
            bundled_binary = os.path.join(root_dir, "lib", "medaka", "mac", "medaka_consensus")
        elif sys.platform == "linux":  # Linux
            bundled_binary = os.path.join(root_dir, "lib", "medaka", "linux", "medaka_consensus")
        else:
            raise RuntimeError(f"Unsupported platform: {sys.platform}")
            
        # If bundled binary exists and is executable, use it
        if os.path.exists(bundled_binary) and os.access(bundled_binary, os.X_OK):
            print(f"Using bundled Medaka binary: {bundled_binary}")
            return bundled_binary
            
        # Otherwise try system binary
        system_binary = shutil.which("medaka_consensus")
        if system_binary:
            print(f"Using system Medaka binary: {system_binary}")
            return system_binary
            
        raise RuntimeError("Medaka binary not found in bundled location or system PATH")

    def _run_medaka(self, assembly_file, reads_file, output_dir, threads=None, model='r941_min_sup_g507'):
        """Run Medaka polishing on assembly
        
        Args:
            assembly_file (str): Path to assembly FASTA file
            reads_file (str): Path to reads FASTQ file
            output_dir (str): Output directory for Medaka
            threads (int, optional): Number of threads to use. Defaults to None.
            model (str, optional): Medaka model to use. Defaults to 'r941_min_sup_g507'.
            
        Returns:
            str: Path to polished assembly file
        """
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Use default threads if none specified
        if threads is None:
            threads = self.default_threads
            
        # Build Medaka command with spoa disabled
        cmd = [
            'medaka_consensus',
            '-i', reads_file,
            '-d', assembly_file,
            '-o', output_dir,
            '-t', str(threads),
            '-m', model,
            '-g'  # Disable gap filling which uses spoa
        ]
        
        print(f"Running Medaka command: {' '.join(cmd)}")
        
        # Run Medaka
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"Medaka stdout:\n{result.stdout}")
            if result.stderr:
                print(f"Medaka stderr:\n{result.stderr}")
        except subprocess.CalledProcessError as e:
            print(f"Medaka error output:\n{e.stderr}")
            raise RuntimeError(f"Medaka polishing failed: {e}")
            
        # Return path to polished assembly
        polished_file = os.path.join(output_dir, 'consensus.fasta')
        if not os.path.exists(polished_file):
            raise RuntimeError("Medaka completed but polished assembly file not found")
            
        return polished_file

    def _run_metabat(self, assembly_file, threads):
        """Run MetaBAT2 binning"""
        output_dir = os.path.join(self.output_dir, "metabat2")
        return self.metabat_runner.run_binning(
            assembly_file=assembly_file,
            output_dir=output_dir,
            threads=threads
        )

    def _check_existing_flye_assembly(self):
        """Check if a completed Flye assembly exists for this sample"""
        flye_dir = os.path.join(self.sample_dir, "flye_out")
        log_file = os.path.join(flye_dir, "flye.log")
        
        if not os.path.exists(log_file):
            return None
            
        # Check log file for successful completion
        try:
            with open(log_file, 'r') as f:
                log_content = f.read()
                
            # Split into separate assembly runs (if multiple exist)
            assembly_runs = log_content.split("/Users/timo/PycharmProjects/MMonitor/desktop/lib/Flye-")
            
            # Check the last run first, then go backwards
            for run in reversed(assembly_runs):
                lines = run.splitlines()
                for line in lines:
                    if "INFO: Final assembly:" in line:
                        # Extract assembly path from log
                        assembly_path = line.split("INFO: Final assembly:")[-1].strip()
                        if os.path.exists(assembly_path):
                            assembly_info = os.path.join(flye_dir, "assembly_info.txt")
                            print(f"\nFound completed assembly in log:")
                            print(f"Assembly path: {assembly_path}")
                            return {
                                "assembly": assembly_path,
                                "assembly_info": assembly_info if os.path.exists(assembly_info) else None,
                                "output_dir": flye_dir
                            }
        except Exception as e:
            print(f"Error reading Flye log: {e}")
            
        return None

    def _run_checkm2(self, bins_dir, threads):
        """Run CheckM2 quality assessment"""
        try:
            if not self._check_and_install_checkm2_db():
                print("CheckM2 database not available, skipping quality assessment")
                return None

            output_dir = os.path.join(self.output_dir, "checkm2")
            os.makedirs(output_dir, exist_ok=True)

            # Get list of proper bin files (exclude tooShort and unbinned)
            bin_files = [f for f in os.listdir(bins_dir) if f.endswith('.fa') 
                        and not any(x in f for x in ['tooShort', 'unbinned'])]
            
            if not bin_files:
                print("No valid bin files found for CheckM2 analysis")
                return None

            # Create a temporary directory for filtered bins
            temp_bins_dir = os.path.join(output_dir, "filtered_bins")
            os.makedirs(temp_bins_dir, exist_ok=True)

            # Copy valid bin files to temp directory
            for bin_file in bin_files:
                src = os.path.join(bins_dir, bin_file)
                dst = os.path.join(temp_bins_dir, bin_file)
                shutil.copy2(src, dst)

            cmd = [
                'checkm2',
                'predict',
                '--threads', str(threads),
                '--input', temp_bins_dir,
                '--output-directory', output_dir
            ]

            print(f"Running CheckM2 command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"CheckM2 error: {result.stderr}")
                return None

            # Clean up temporary directory
            shutil.rmtree(temp_bins_dir)

            return {
                'quality_dir': output_dir,
                'predictions': os.path.join(output_dir, 'predictions.tsv')
            }

        except Exception as e:
            print(f"Error running CheckM2: {str(e)}")
            return None

    def _run_gtdbtk(self, bins_dir, threads):
        """Run GTDB-Tk classification"""
        try:
            if not self._check_and_install_gtdbtk_db():
                print("GTDB-TK database not available, skipping classification")
                return None

            output_dir = os.path.join(self.sample_dir, "gtdbtk")
            os.makedirs(output_dir, exist_ok=True)

            cmd = [
                'gtdbtk',
                'classify_wf',
                '--genome_dir', bins_dir,
                '--out_dir', output_dir,
                '--cpus', str(threads),
                '--extension', 'fa'
            ]

            print(f"Running GTDB-TK command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"GTDB-TK error: {result.stderr}")
                return None

            return {
                'taxonomy_dir': output_dir,
                'summary': os.path.join(output_dir, 'gtdbtk.summary.tsv')
            }

        except Exception as e:
            print(f"Error running GTDB-TK: {str(e)}")
            return None

    def _run_bakta(self, bins_dir, threads):
        """Run Bakta annotation"""
        try:
            db_path = self._check_bakta_db()
            if not db_path:
                print("Bakta database not available, skipping annotation")
                return None

            output_dir = os.path.join(self.sample_dir, "bakta")
            os.makedirs(output_dir, exist_ok=True)

            # Run Bakta on each bin
            results = {}
            for bin_file in os.listdir(bins_dir):
                if bin_file.endswith('.fa') or bin_file.endswith('.fasta'):
                    bin_path = os.path.join(bins_dir, bin_file)
                    bin_name = os.path.splitext(bin_file)[0]
                    bin_output = os.path.join(output_dir, bin_name)
                    
                    cmd = [
                        'bakta',
                        '--db', db_path,
                        '--threads', str(threads),
                        '--output', bin_output,
                        '--prefix', bin_name,
                        bin_path
                    ]

                    print(f"Running Bakta on {bin_file}: {' '.join(cmd)}")
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    
                    if result.returncode != 0:
                        print(f"Bakta error for {bin_file}: {result.stderr}")
                        continue

                    results[bin_name] = {
                        'gff3': f"{bin_output}.gff3",
                        'faa': f"{bin_output}.faa",
                        'ffn': f"{bin_output}.ffn",
                        'tsv': f"{bin_output}.tsv",
                        'gbff': f"{bin_output}.gbff",
                        'embl': f"{bin_output}.embl",
                    }

            return {
                'annotation_dir': output_dir,
                'bin_results': results
            } if results else None

        except Exception as e:
            print(f"Error running Bakta: {str(e)}")
            return None

    def _check_bakta_db(self):
        """Check if Bakta database exists and download if needed
        
        Returns:
            str: Path to the Bakta database
        """
        # Get the database path from config or use default
        db_path = os.path.join(os.getcwd(), "db-light")
        
        # Check if database exists
        if not os.path.exists(db_path):
            if self.progress_callback:
                response = self.progress_callback({
                    'type': 'user_input',
                    'message': 'Bakta database not found. Would you like to download it? This may take a while (~1.5GB).',
                    'options': ['yes', 'no']
                })
                if response != 'yes':
                    raise RuntimeError("Bakta database is required but user declined download")
            
            print("Downloading Bakta database...")
            
            # Use bakta_db download command
            cmd = [
                'bakta_db',
                'download',
                '--type', 'light'
            ]
            
            try:
                # Run the command and ignore AMRFinderPlus errors
                process = subprocess.run(cmd, capture_output=True, text=True)
                if process.returncode != 0 and "AMRFinderPlus failed" not in process.stderr:
                    raise subprocess.CalledProcessError(process.returncode, cmd, process.stderr)
                print("Bakta database downloaded successfully")
            except subprocess.CalledProcessError as e:
                error_msg = f"Error downloading Bakta database: {str(e)}"
                logger.error(error_msg)
                raise RuntimeError(error_msg)
            except Exception as e:
                error_msg = f"Unexpected error downloading Bakta database: {str(e)}"
                logger.error(error_msg)
                raise RuntimeError(error_msg)
        
        return db_path

    def run_bakta_annotation(self, assembly_file, output_dir, threads=None):
        """Run Bakta annotation on the assembly

        Args:
            assembly_file (str): Path to assembly file
            output_dir (str): Output directory for Bakta results
            threads (int, optional): Number of threads to use. Defaults to None.

        Returns:
            str: Path to the Bakta output directory
        """
        if not assembly_file or not os.path.exists(assembly_file):
            error_msg = f"Assembly file not found: {assembly_file}"
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
            
        os.makedirs(output_dir, exist_ok=True)
        
        # Use default threads if not specified
        threads = threads or self.default_threads
        
        # Check and download Bakta database if needed
        db_path = self._check_bakta_db()
        
        # Build the command
        cmd = [
            str(self.bakta_path),  # Convert to string in case it's a Path object
            '--db', str(db_path),
            '--threads', str(threads),
            '--output', str(output_dir),
            '--prefix', str(self.sample_name),
            str(assembly_file)  # Convert to string in case it's a Path object
        ]
        
        # Run Bakta
        try:
            subprocess.run(cmd, check=True)
            return output_dir
        except subprocess.CalledProcessError as e:
            error_msg = f"Error running Bakta: {str(e)}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

    def run_flye(self, input_files, sample_name, output_dir, threads):
        """Run Flye assembly"""
        flye_out = os.path.join(output_dir, "flye_out")
        os.makedirs(flye_out, exist_ok=True)

        # Get the Flye bin directory
        flye_bin_dir = os.path.dirname(self.flye_path)
        flye_lib_dir = os.path.dirname(os.path.dirname(flye_bin_dir))

        # Create a symlink to minimap2 with the name flye-minimap2 in Flye's bin directory
        try:
            # First check if minimap2 exists in Flye's lib directory
            flye_minimap2 = os.path.join(flye_lib_dir, "lib", "minimap2", "minimap2")
            if os.path.exists(flye_minimap2):
                source_minimap2 = flye_minimap2
                print(f"Using minimap2 from Flye's lib directory: {flye_minimap2}")
            else:
                source_minimap2 = "/usr/bin/minimap2"
                print(f"Using system minimap2: {source_minimap2}")

            # Create the symlink with the correct name
            flye_minimap2_target = os.path.join(flye_bin_dir, "flye-minimap2")
            if not os.path.exists(flye_minimap2_target):
                os.symlink(source_minimap2, flye_minimap2_target)
                print(f"Created symlink: {source_minimap2} -> {flye_minimap2_target}")
                os.chmod(flye_minimap2_target, 0o755)
        except Exception as e:
            print(f"Warning: Could not create minimap2 symlink: {e}")

        # Add Flye's bin directory to PATH
        os.environ['PATH'] = f"{flye_bin_dir}:{os.environ['PATH']}"
        print(f"Updated PATH: {os.environ['PATH']}")

        # Run Flye
        cmd = [self.flye_path, "--nano-raw"] + input_files + ["--out-dir", flye_out, "--threads", str(threads)]
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True, env=os.environ)
            print(f"Flye stdout: {result.stdout}")
        except subprocess.CalledProcessError as e:
            print(f"Error running Flye: {e}")
            print(f"Flye stderr: {e.stderr}")
            raise

        return os.path.join(flye_out, "assembly.fasta")

    def _check_required_tools(self):
        """Check for required tools and raise error if essential ones are missing"""
        tools = ['minimap2', 'flye', 'medaka', 'bakta', 'gtdbtk', 'checkm2', 'metabat2']
        for tool in tools:
            tool_path = self.tool_installer.check_tool(tool)
            if not tool_path:
                print(f"Tool {tool} not available, some pipeline steps may be skipped")
                if tool in ['minimap2', 'flye']:  # These are essential
                    raise RuntimeError(f"Essential tool {tool} is not available")