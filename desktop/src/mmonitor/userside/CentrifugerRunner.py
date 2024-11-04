import os
import subprocess
import gzip
import shutil
import multiprocessing
import logging
from build_mmonitor_pyinstaller import ROOT
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import concurrent
from threading import Lock
import gzip
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager
import platform
import subprocess

class CentrifugerRunner:
    def __init__(self):
        self.centrifuger_path = "centrifuger"
        self.logger = logging.getLogger('timestamp')
        self.pipeline_out = os.path.join(ROOT, "src", "resources", "pipeline_out")
        self.temp_dir = os.path.join(ROOT, "src", "resources", "temp")
        os.makedirs(self.temp_dir, exist_ok=True)

    def get_files_from_folder(self, folder):
        """Get all FASTQ files from a folder"""
        print(f"\nScanning folder for FASTQ files: {folder}")
        files = []
        for root, _, filenames in os.walk(folder):
            for filename in filenames:
                if filename.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                    file_path = os.path.join(root, filename)
                    files.append(file_path)
                    print(f"Found file: {file_path}")
        
        print(f"Total files found: {len(files)}")
        return sorted(files)

    def concatenate_gzipped_files(self, files, output_file):
        """Concatenate gzipped FASTQ files with detailed logging"""
        print(f"\nConcatenating {len(files)} files:")
        for file in files:
            print(f"  - {file}")
        
        print(f"\nCreating temporary concatenated file: {output_file}")
        
        try:
            total_size = 0
            with gzip.open(output_file, 'wb') as outfile:
                for file in files:
                    try:
                        if file.endswith('.gz'):
                            with gzip.open(file, 'rb') as infile:
                                # Read and write in chunks to handle large files
                                chunk_size = 1024 * 1024  # 1MB chunks
                                while True:
                                    chunk = infile.read(chunk_size)
                                    if not chunk:
                                        break
                                    outfile.write(chunk)
                                    total_size += len(chunk)
                        else:
                            with open(file, 'rb') as infile:
                                while True:
                                    chunk = infile.read(chunk_size)
                                    if not chunk:
                                        break
                                    outfile.write(gzip.compress(chunk))
                                    total_size += len(chunk)
                        print(f"Added {file}")
                    except Exception as e:
                        print(f"Error processing file {file}: {e}")
                        continue
            
            print(f"Concatenation complete. Total size: {total_size:,} bytes")
            return True
        except Exception as e:
            print(f"Error during concatenation: {e}")
            return False

    def run_centrifuger(self, input_file, sample_name, db_path):
        """Run Centrifuger analysis on input file"""
        os.makedirs(self.pipeline_out, exist_ok=True)
        
        # Create sample-specific output directory
        sample_out_dir = os.path.join(self.pipeline_out, sample_name)
        os.makedirs(sample_out_dir, exist_ok=True)
        
        
        # Output files
        output_file = os.path.join(sample_out_dir, f"{sample_name}_centrifuger_classifications.txt")
        report_file = os.path.join(sample_out_dir, f"{sample_name}_centrifuger_report.tsv")
        log_file = os.path.join(sample_out_dir, f"{sample_name}_centrifuger.log")
        
        print(f"\nStarting Centrifuger analysis for sample: {sample_name}")
        print(f"Database path: {db_path}")
        print(f"Input file: {input_file}")
        print(f"Output directory: {sample_out_dir}")
        
        # Construct command with correct parameters for Centrifuger
        cmd = [
            "centrifuger",
            "-x", db_path,           # Index prefix
            "-u", input_file,        # Single-end read file
            "-t", str(multiprocessing.cpu_count()),  # Number of threads
            "-k", "1"                # Report up to k distinct assignments
        ]
        
        print(f"\nRunning Centrifuger command: {' '.join(cmd)}")
        
        try:
            # Run Centrifuger and capture all output
            with open(log_file, 'w') as log:
                process = subprocess.run(
                    cmd,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                # Save stdout to output file
                with open(output_file, 'w') as f:
                    f.write(process.stdout)
                
                # Log both stdout and stderr
                log.write("=== STDOUT ===\n")
                log.write(process.stdout)
                log.write("\n=== STDERR ===\n")
                log.write(process.stderr)
                
                # Also print to console
                print("\nCentrifuger STDOUT:")
                print(process.stdout)
                print("\nCentrifuger STDERR:")
                print(process.stderr)
            
            # Generate report using centrifuger-kreport and capture its output
            print("\nGenerating Kraken-style report...")
            kreport_cmd = [
                "centrifuger-kreport",
                "-x", db_path,
                output_file
            ]
            
            try:
                kreport_result = subprocess.run(
                    kreport_cmd,
                    check=True,
                    capture_output=True,
                    text=True
                )
                
                # Write kreport output to file
                with open(report_file, 'w') as f:
                    f.write(kreport_result.stdout)
                
                # Append kreport logs
                with open(log_file, 'a') as log:
                    log.write("\n=== KREPORT STDOUT ===\n")
                    log.write(kreport_result.stdout)
                    log.write("\n=== KREPORT STDERR ===\n")
                    log.write(kreport_result.stderr)
                
                print("Centrifuger analysis completed successfully")
                print(f"Results saved to: {sample_out_dir}")
                return True
                
            except subprocess.CalledProcessError as e:
                print(f"Error generating Kraken-style report: {e}")
                print(f"kreport stderr: {e.stderr}")
                return False
                
        except subprocess.CalledProcessError as e:
            print(f"Error running Centrifuger: {e}")
            print(f"Centrifuger stderr: {e.stderr}")
            if "command not found" in str(e):
                print("Centrifuger is not installed or not in the system PATH")
            return False
            
        except Exception as e:
            print(f"Unexpected error running Centrifuger: {e}")
            return False

    def run_multi_sample(self, sample_sheet, db_path):
        """Run Centrifuger on multiple samples using a sample sheet"""
        cmd = [
            "centrifuger",
            "-x", db_path,
            "--sample-sheet", sample_sheet,
            "-p", str(multiprocessing.cpu_count())
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info("Multi-sample Centrifuger analysis completed successfully")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error running multi-sample Centrifuger: {e.stderr}")
            raise

    def create_sample_sheet(self, samples, output_file):
        """Create a sample sheet for multi-sample analysis"""
        with open(output_file, 'w') as f:
            for sample_name, file_path in samples.items():
                output_path = os.path.join(self.pipeline_out, f"{sample_name}_centrifuger_out.txt")
                report_path = os.path.join(self.pipeline_out, f"{sample_name}_centrifuger_report.tsv")
                f.write(f"1\t{file_path}\t\t{output_path}\t{report_path}\n")
        return output_file

    def make_kraken_report(self, centrifuger_out, db_path, output_file):
        """Convert Centrifuger output to Kraken-style report"""
        cmd = [
            "centrifuger-kreport",
            "-x", db_path,
            centrifuger_out
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            with open(output_file, 'w') as f:
                f.write(result.stdout)
            self.logger.info("Kraken report generated successfully")
            return output_file
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error generating Kraken report: {e.stderr}")
            raise
