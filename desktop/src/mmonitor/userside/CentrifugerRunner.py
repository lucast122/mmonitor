import os
import subprocess
import gzip
import shutil
import multiprocessing
import logging
import threading
import time
from ..paths import SRC_DIR, RESOURCES_DIR, ROOT, PROJECT_ROOT, CENTRIFUGER_DIR
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import concurrent
from threading import Lock
import gzip
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager
import platform
import subprocess
import datetime
import sys

class CentrifugerRunner:
    def __init__(self):
        # Find the custom centrifuger binary
        self.centrifuger_path = self._find_custom_centrifuger()
        self.logger = logging.getLogger('timestamp')
        self.pipeline_out = os.path.join(SRC_DIR, "resources", "pipeline_out")
        self.temp_dir = os.path.join(SRC_DIR, "resources", "temp")
        os.makedirs(self.temp_dir, exist_ok=True)
        self.centrifuger_db_path = None
        self.stop_log_monitoring = False
        self.log_callback = None

    def _find_custom_centrifuger(self):
        """Find the custom centrifuger binary in possible locations"""
        from ..paths import PROJECT_ROOT, CENTRIFUGER_DIR
        
        # First, check if we have a specifically defined CENTRIFUGER_DIR
        if CENTRIFUGER_DIR and os.path.exists(CENTRIFUGER_DIR):
            print(f"Using CENTRIFUGER_DIR from paths.py: {CENTRIFUGER_DIR}")
            return CENTRIFUGER_DIR
        
        # Check for bundled app with _internal structure first (most specific path)
        internal_path = None
        if getattr(sys, 'frozen', False):
            # Determine app location
            if sys.platform == 'darwin':
                # For macOS .app bundle
                app_dir = os.path.dirname(sys.executable)
                if 'Contents/MacOS' in app_dir:
                    # We're in the .app bundle, check for _internal
                    app_root = os.path.abspath(os.path.join(app_dir, '..', '..', '..'))
                    dist_dir = os.path.dirname(app_root)
                    internal_path = os.path.join(dist_dir, "MMonitor", "_internal", "lib", "centrifuger-realtime", "centrifuger", "centrifuger")
                    if os.path.exists(internal_path):
                        print(f"Found centrifuger in _internal: {internal_path}")
                        return internal_path
                else:
                    # We might be inside the _internal directory itself
                    app_dir_parts = app_dir.split(os.sep)
                    if '_internal' in app_dir_parts:
                        # Determine the _internal directory and look from there
                        idx = app_dir_parts.index('_internal')
                        internal_dir = os.sep.join(app_dir_parts[:idx+1])
                        internal_path = os.path.join(internal_dir, "lib", "centrifuger-realtime", "centrifuger", "centrifuger")
                        if os.path.exists(internal_path):
                            print(f"Found centrifuger in _internal (from path): {internal_path}")
                            return internal_path
        
        # Check project lib directory (absolute path)
        custom_binary = os.path.join(PROJECT_ROOT, "lib", "centrifuger-realtime", "centrifuger", "centrifuger")
        if os.path.exists(custom_binary):
            print(f"Found custom centrifuger at: {custom_binary}")
            return custom_binary
            
        # Check relative lib directory 
        relative_path = os.path.join("lib", "centrifuger-realtime", "centrifuger", "centrifuger")
        if os.path.exists(relative_path):
            abs_path = os.path.abspath(relative_path)
            print(f"Found custom centrifuger at: {abs_path}")
            return abs_path
            
        # For development environment
        desktop_path = os.path.join(PROJECT_ROOT, "desktop", "lib", "centrifuger-realtime", "centrifuger", "centrifuger")
        if os.path.exists(desktop_path):
            print(f"Found custom centrifuger at: {desktop_path}")
            return desktop_path
            
        # For bundled application (Contents/Resources/bin)
        if sys.platform == 'darwin' and getattr(sys, 'frozen', False):
            app_dir = os.path.dirname(sys.executable)
            if 'Contents/MacOS' in app_dir:
                resources_dir = os.path.abspath(os.path.join(app_dir, '..', 'Resources'))
                resources_bin = os.path.join(resources_dir, "lib", "centrifuger-realtime", "centrifuger")
                if os.path.exists(resources_bin):
                    print(f"Found bundled centrifuger at: {resources_bin}")
                    return resources_bin
            
        # Check in src/bin folder
        src_bin = os.path.join(SRC_DIR, "bin", "centrifuger")
        if os.path.exists(src_bin):
            print(f"Found src/bin centrifuger at: {src_bin}")
            return src_bin
        
        # Final attempt - check specific paths for bundled apps
        if getattr(sys, 'frozen', False):
            # Path pattern found in user's message
            specific_path = "/Users/timo/PycharmProjects/MMonitor/dist/MMonitor/_internal/lib/centrifuger-realtime/centrifuger/centrifuger"
            if os.path.exists(specific_path):
                print(f"Found specific bundled centrifuger at: {specific_path}")
                return specific_path
            
            # Try to generalize the path
            if sys.platform == 'darwin':
                # Check in relative paths from executable
                app_dir = os.path.dirname(sys.executable)
                if 'Contents/MacOS' in app_dir:
                    # We're in a .app bundle
                    project_dir = os.path.abspath(os.path.join(app_dir, '..', '..', '..', '..'))
                    dist_dir = os.path.join(project_dir, "dist")
                    if os.path.exists(dist_dir):
                        for subdir in os.listdir(dist_dir):
                            candidate = os.path.join(dist_dir, subdir, "_internal", "lib", "centrifuger-realtime", "centrifuger", "centrifuger")
                            if os.path.exists(candidate):
                                print(f"Found candidate centrifuger at: {candidate}")
                                return candidate
            
        # Fall back to PATH-based binary as last resort
        print("WARNING: Could not find custom centrifuger binary, falling back to PATH-based centrifuger")
        print("Searched paths:")
        print(f"  PROJECT_ROOT: {PROJECT_ROOT}")
        print(f"  CENTRIFUGER_DIR: {CENTRIFUGER_DIR}")
        print(f"  Custom binary: {custom_binary}")
        print(f"  Desktop path: {desktop_path}")
        print(f"  Internal path: {internal_path}")
        
        # Try to use 'which' to locate centrifuger in PATH
        try:
            import subprocess
            which_result = subprocess.run(["which", "centrifuger"], capture_output=True, text=True)
            if which_result.returncode == 0:
                path = which_result.stdout.strip()
                print(f"Found centrifuger in PATH: {path}")
                return path
        except:
            pass
        
        return "centrifuger"

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

    def get_full_db_path(self, db_path):
        """Get the full path to the Centrifuger database"""
        if not db_path:
            raise ValueError("Database path not provided")
            
        # If it's already a full path, check if it's a directory or a file prefix
        if os.path.isabs(db_path):
            # If it's a directory, look for .cfr files and return the prefix
            if os.path.isdir(db_path):
                print(f"Looking for Centrifuger database files in directory: {db_path}")
                # Find all .cfr files
                cfr_files = []
                for file in os.listdir(db_path):
                    if file.endswith('.cfr'):
                        cfr_files.append(file)
                
                if cfr_files:
                    # Sort to get predictable ordering
                    cfr_files.sort()
                    print(f"Found {len(cfr_files)} .cfr files: {cfr_files}")
                    
                    # Extract prefix by removing the .N.cfr part
                    # For example, if we have cfr_ncbi.1.cfr, get cfr_ncbi
                    first_file = cfr_files[0]
                    prefix_parts = first_file.split('.')
                    if len(prefix_parts) >= 3:  # Should be like prefix.number.cfr
                        prefix = '.'.join(prefix_parts[:-2])
                        full_prefix_path = os.path.join(db_path, prefix)
                        print(f"Using database prefix: {full_prefix_path}")
                        return full_prefix_path
            
            # If it's already a file prefix or we couldn't find .cfr files, return it as is
            return db_path
            
        # Otherwise, look for it in the library directory
        lib_db_path = os.path.join(ROOT, "lib", "centrifuger_db", db_path)
        if os.path.exists(lib_db_path):
            return lib_db_path
            
        # If not found, return the original path
        return db_path

    def monitor_log_file(self, log_file, callback=None):
        """Monitor log file for updates and call the callback function with new lines"""
        self.stop_log_monitoring = False
        self.log_callback = callback
        
        def _monitor():
            last_position = 0
            last_check_time = time.time()
            last_update_time = time.time()
            last_db_check_time = time.time()
            last_output_size = 0
            last_species_check_time = time.time()
            output_dir = os.path.dirname(log_file)
            classifications_file = None
            report_file = None
            last_processed_species = set()
            species_counts = {}
            
            # Find the classifications file and report file in the same directory
            for f in os.listdir(output_dir):
                if f.endswith("_centrifuger_classifications.txt"):
                    classifications_file = os.path.join(output_dir, f)
                if f.endswith("_centrifuger_report.tsv"):
                    report_file = os.path.join(output_dir, f)
            
            print(f"Starting log file monitoring for: {log_file}")
            if not os.path.exists(log_file):
                print(f"WARNING: Log file does not exist yet: {log_file}")
            
            while not self.stop_log_monitoring:
                try:
                    current_time = time.time()
                    
                    # Print a status message every 15 seconds even if no updates
                    if current_time - last_check_time > 15:
                        if os.path.exists(log_file):
                            file_size = os.path.getsize(log_file)
                            print(f"Log monitor checking file: {log_file} (size: {file_size} bytes, last pos: {last_position})")
                        else:
                            print(f"Log file still does not exist: {log_file}")
                        last_check_time = current_time
                    
                    # Check for database updates after processing completes
                    if classifications_file and os.path.exists(classifications_file):
                        current_output_size = os.path.getsize(classifications_file)
                        
                        # If file has grown since last check, data has been processed
                        if current_output_size > last_output_size:
                            db_update_msg = f"New classification data detected: {current_output_size - last_output_size} bytes added"
                            print(db_update_msg)
                            last_output_size = current_output_size
                            
                            # Count classified reads and extract species information
                            try:
                                # Parse classification file to get read count and species
                                new_species_data = self._parse_classifications(classifications_file, last_processed_species)
                                read_count = new_species_data["total_reads"]
                                new_species = new_species_data["new_species"]
                                updated_species_counts = new_species_data["species_counts"]
                                
                                # Update our tracking of species counts
                                species_counts.update(updated_species_counts)
                                
                                # Create detailed message about species
                                species_msg = ""
                                if new_species:
                                    species_msg = f"\nNEW SPECIES DETECTED:"
                                    for species, count in sorted([(s, updated_species_counts[s]) for s in new_species], 
                                                               key=lambda x: x[1], reverse=True):
                                        species_msg += f"\n  - {species}: {count} reads"
                                    
                                # Add summary of top species
                                top_species = sorted(species_counts.items(), key=lambda x: x[1], reverse=True)[:10]
                                species_msg += f"\n\nTOP 10 SPECIES (out of {len(species_counts)}):"
                                for species, count in top_species:
                                    percent = (count / read_count * 100) if read_count > 0 else 0
                                    species_msg += f"\n  - {species}: {count} reads ({percent:.2f}%)"
                                
                                # Create full message
                                db_update_msg = f"Database update: {read_count} total classified reads available for upload{species_msg}"
                                print(db_update_msg)
                                if self.log_callback:
                                    self.log_callback(db_update_msg + "\n")
                                
                                # Update our tracking of what we've seen
                                last_processed_species.update(new_species)
                                
                            except Exception as e:
                                print(f"Error analyzing classifications: {e}")
                        
                        # Check report file for more detailed taxonomic information
                        elif current_time - last_species_check_time > 30 and os.path.exists(report_file):
                            try:
                                # Only process report if it has content
                                if os.path.getsize(report_file) > 0:
                                    abundance_data = self._parse_report_file(report_file)
                                    if abundance_data and len(abundance_data) > 0:
                                        # Show phylum-level distribution
                                        phylum_msg = "\nTAXONOMIC DISTRIBUTION (PHYLUM LEVEL):"
                                        for phylum, percent in abundance_data["phylum"].items():
                                            phylum_msg += f"\n  - {phylum}: {percent:.2f}%"
                                        
                                        # Show genus-level distribution (top 10)
                                        genus_items = sorted(abundance_data["genus"].items(), 
                                                           key=lambda x: x[1], reverse=True)[:10]
                                        genus_msg = "\n\nTOP 10 GENERA:"
                                        for genus, percent in genus_items:
                                            genus_msg += f"\n  - {genus}: {percent:.2f}%"
                                        
                                        # Show species-level distribution (top 10)
                                        species_items = sorted(abundance_data["species"].items(), 
                                                            key=lambda x: x[1], reverse=True)[:10]
                                        species_msg = "\n\nTOP 10 SPECIES:"
                                        for species, percent in species_items:
                                            species_msg += f"\n  - {species}: {percent:.2f}%"
                                        
                                        full_msg = f"Taxonomic abundance summary:{phylum_msg}{genus_msg}{species_msg}"
                                        print(full_msg)
                                        if self.log_callback:
                                            self.log_callback(full_msg + "\n")
                            except Exception as e:
                                print(f"Error analyzing taxonomic report: {e}")
                            
                            last_species_check_time = current_time
                        
                        # Every 30 seconds, check if we've processed new data and should update database
                        elif current_time - last_db_check_time > 30 and current_output_size > 0:
                            status_msg = f"Database status check: {current_output_size} bytes of classification data available"
                            print(status_msg)
                            if self.log_callback:
                                self.log_callback(status_msg + "\n")
                            last_db_check_time = current_time
                    
                    if os.path.exists(log_file):
                        file_size = os.path.getsize(log_file)
                        if file_size > last_position:  # Only read if file has grown
                            with open(log_file, 'r') as f:
                                f.seek(last_position)
                                new_content = f.read()
                                if new_content:
                                    print(f"Log monitor: Read {len(new_content)} new bytes from log file")
                                    if self.log_callback:
                                        self.log_callback(new_content)
                                    else:
                                        print(new_content, end='')
                                    last_position = f.tell()
                                    last_update_time = current_time
                        elif current_time - last_update_time > 60:  # Check every minute
                            # This would be a good time to update the database if needed
                            status_msg = f"Waiting for new files... (Current log size: {file_size} bytes)"
                            print(status_msg)
                            if self.log_callback:
                                self.log_callback(status_msg + "\n")
                            last_update_time = current_time
                    time.sleep(1)  # Check every second
                except Exception as e:
                    error_msg = f"Error monitoring log file: {e}"
                    print(error_msg)
                    if self.log_callback:
                        self.log_callback(f"ERROR: {error_msg}\n")
                    time.sleep(5)  # Longer delay if there's an error
            
            print(f"Log monitoring stopped for: {log_file}")
        
        # Start the monitoring thread
        monitor_thread = threading.Thread(target=_monitor, name="LogMonitorThread")
        monitor_thread.daemon = True
        monitor_thread.start()
        
        print(f"Log monitor thread started for: {log_file}")
        return monitor_thread

    def _parse_classifications(self, classifications_file, last_processed_species):
        """Parse the classifications file to extract species information"""
        result = {
            "total_reads": 0,
            "new_species": set(),
            "species_counts": {}
        }
        
        try:
            with open(classifications_file, 'r') as f:
                # Skip header line
                header = f.readline()
                
                # Process all other lines
                for line in f:
                    result["total_reads"] += 1
                    parts = line.strip().split('\t')
                    
                    # Format depends on version, but typically:
                    # readID, seq-id, taxID, score, 2nd-best-score, hit-length, query-length, num-matches
                    if len(parts) >= 3:
                        tax_id = parts[2]
                        
                        # Get species name (use taxID if name lookup fails)
                        species_name = self._get_species_name(tax_id)
                        
                        # Count occurrences
                        if species_name in result["species_counts"]:
                            result["species_counts"][species_name] += 1
                        else:
                            result["species_counts"][species_name] = 1
                            
                        # Check if this is a new species
                        if species_name not in last_processed_species:
                            result["new_species"].add(species_name)
            
            return result
        except Exception as e:
            print(f"Error parsing classifications file: {e}")
            return result

    def _get_species_name(self, tax_id):
        """Get species name from tax ID - simple implementation"""
        # This is a placeholder - ideally we would use a real taxonomy database
        # For now, just return the tax ID as a string
        try:
            # If tax_id is 0 or not provided, it's unclassified
            if not tax_id or tax_id == '0':
                return "Unclassified"
                
            return f"taxID:{tax_id}"
        except:
            return f"Unknown"

    def _parse_report_file(self, report_file):
        """Parse the Kraken-style report to extract taxonomic abundances"""
        result = {
            "phylum": {},
            "genus": {},
            "species": {}
        }
        
        try:
            with open(report_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 6:
                        continue
                        
                    # Kraken-style report format:
                    # percent, clade_reads, direct_reads, rank_code, taxID, name
                    percent = float(parts[0])
                    rank = parts[3]
                    name = parts[5].strip()
                    
                    # Remove leading spaces and dashes from name
                    clean_name = name.lstrip(' -')
                    
                    # Store information by rank
                    if rank == 'P':  # Phylum
                        result["phylum"][clean_name] = percent
                    elif rank == 'G':  # Genus
                        result["genus"][clean_name] = percent
                    elif rank == 'S':  # Species
                        result["species"][clean_name] = percent
                        
            return result
        except Exception as e:
            print(f"Error parsing report file: {e}")
            return result

    def stop_monitoring(self):
        """Stop the log monitoring thread"""
        self.stop_log_monitoring = True

    def run_centrifuger(self, input_file, sample_name, db_path, realtime_mode=False, watch_dir=None, min_files=5, log_callback=None):
        """Run Centrifuger analysis on input file"""
        os.makedirs(self.pipeline_out, exist_ok=True)
        
        # Create sample-specific output directory
        sample_out_dir = os.path.join(self.pipeline_out, sample_name)
        os.makedirs(sample_out_dir, exist_ok=True)
        
        try:
            # Get full database path with correct prefix
            full_db_path = self.get_full_db_path(db_path)
            
            # Output files
            output_file = os.path.join(sample_out_dir, f"{sample_name}_centrifuger_classifications.txt")
            report_file = os.path.join(sample_out_dir, f"{sample_name}_centrifuger_report.tsv")
            log_file = os.path.join(sample_out_dir, f"{sample_name}_centrifuger.log")
            
            print(f"\n{'='*50}")
            print(f"CENTRIFUGER ANALYSIS: {sample_name}")
            print(f"{'='*50}")
            print(f"Starting Centrifuger analysis for sample: {sample_name}")
            print(f"Using custom centrifuger binary: {self.centrifuger_path}")
            print(f"Database path: {full_db_path}")
            
            # In realtime mode, set up the watch directory properly
            if realtime_mode:
                # If watch_dir is not provided or doesn't exist, create it
                if not watch_dir or not os.path.exists(watch_dir):
                    watch_dir = os.path.join(sample_out_dir, "watch_dir")
                    os.makedirs(watch_dir, exist_ok=True)
                    print(f"Created watch directory: {watch_dir}")
                
                print(f"REALTIME MODE ENABLED")
                print(f"Watch directory: {watch_dir}")
                print(f"Min files threshold: {min_files}")
                print(f"Output directory: {sample_out_dir}")
                print(f"Output classification file: {output_file}")
                print(f"Output report file: {report_file}")
                print(f"Log file: {log_file}")
                
                # Check for existing FASTQ files in the watch directory
                files_in_dir = [f for f in os.listdir(watch_dir) if os.path.isfile(os.path.join(watch_dir, f))]
                fastq_files = [f for f in files_in_dir if f.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz'))]
                print(f"Current files in watch directory: {len(fastq_files)} FASTQ files out of {len(files_in_dir)} total files")
            else:
                print(f"STANDARD MODE")
                print(f"Input file: {input_file}")
                print(f"Output directory: {sample_out_dir}")
                print(f"Output classification file: {output_file}")
                print(f"Output report file: {report_file}")
                print(f"Log file: {log_file}")
            
            # Verify binary exists
            if not os.path.exists(self.centrifuger_path) and not shutil.which(self.centrifuger_path):
                raise FileNotFoundError(f"Centrifuger binary not found: {self.centrifuger_path}")
            
            # Base command with common parameters
            cmd = [
                self.centrifuger_path,
                "-x", full_db_path,           # Index prefix
                "-t", str(multiprocessing.cpu_count()),  # Number of threads
                "-k", "1",                    # Report up to k distinct assignments
                "-o", output_file             # Output file
            ]
            
            if realtime_mode:
                if not watch_dir:
                    raise ValueError("watch_dir is required in realtime mode")
                cmd.extend([
                    "--realtime",
                    "--watch-dir", watch_dir,
                    "--min-files", str(min_files)
                ])
            else:
                # The C++ version uses -u for single-end reads
                cmd.extend(["-u", input_file])
            
            print(f"\nRunning Centrifuger command:")
            print(f"  {' '.join(cmd)}")
            print(f"\nStarting analysis. This may take some time...")
            print(f"Check log file for detailed progress: {log_file}")
            
            # Start monitoring the log file in a separate thread if in realtime mode
            if realtime_mode:
                with open(log_file, 'w') as log:
                    log.write(f"=== CENTRIFUGER ANALYSIS FOR {sample_name} ===\n")
                    log.write(f"Command: {' '.join(cmd)}\n")
                    log.write(f"Started at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    log.write(f"Watch directory: {watch_dir}\n")
                    log.write(f"Min files threshold: {min_files}\n\n")
                
                # Setup result monitoring thread to update progress
                last_file_count = len(fastq_files)
                last_output_size = 0
                if os.path.exists(output_file):
                    last_output_size = os.path.getsize(output_file)
                
                def check_progress():
                    nonlocal last_file_count, last_output_size
                    stop_checking = False
                    
                    while not stop_checking and not self.stop_log_monitoring:
                        try:
                            # Check watch directory for new files
                            current_files = [f for f in os.listdir(watch_dir) if f.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz'))]
                            current_file_count = len(current_files)
                            
                            # Check output file for new classifications
                            current_output_size = 0
                            classified_reads = 0
                            if os.path.exists(output_file):
                                current_output_size = os.path.getsize(output_file)
                                if current_output_size > 0:
                                    try:
                                        with open(output_file, 'r') as f:
                                            classified_reads = sum(1 for _ in f) - 1  # Subtract header
                                    except:
                                        pass
                            
                            progress_info = []
                            
                            # Report if new files have been detected
                            if current_file_count != last_file_count:
                                progress_info.append(f"Files in watch directory: {current_file_count} (changed from {last_file_count})")
                                last_file_count = current_file_count
                            
                            # Report if output has grown
                            if current_output_size != last_output_size:
                                progress_info.append(f"Output file size: {current_output_size} bytes (changed from {last_output_size})")
                                progress_info.append(f"Classified reads: {classified_reads}")
                                last_output_size = current_output_size
                            
                            # Update progress
                            if progress_info and log_callback:
                                log_callback("\n".join(progress_info) + "\n")
                                
                                # Also append to the log file since the centrifuger process isn't writing much
                                with open(log_file, 'a') as log:
                                    log.write(f"\n--- PROGRESS UPDATE: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ---\n")
                                    for info in progress_info:
                                        log.write(f"{info}\n")
                            
                            # Check if we should stop (process might have finished)
                            if not os.path.exists(os.path.join("/proc", str(centrifuger_process.pid))) and platform.system() == "Linux":
                                stop_checking = True
                            
                            time.sleep(10)  # Check every 10 seconds
                        except Exception as e:
                            error_msg = f"Error in progress monitoring: {e}"
                            print(error_msg)
                            if log_callback:
                                log_callback(f"ERROR: {error_msg}\n")
                            time.sleep(15)  # Wait a bit longer after an error
                
                # Start progress checking thread
                progress_thread = threading.Thread(target=check_progress, name="ProgressMonitorThread")
                progress_thread.daemon = True
                progress_thread.start()
            
            # Run Centrifuger as a subprocess
            with open(log_file, 'a' if realtime_mode else 'w') as log:
                if not realtime_mode:
                    log.write(f"=== CENTRIFUGER ANALYSIS FOR {sample_name} ===\n")
                    log.write(f"Command: {' '.join(cmd)}\n")
                    log.write(f"Started at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                try:
                    # Use Popen to be able to capture output in real-time
                    if realtime_mode:
                        centrifuger_process = subprocess.Popen(
                            cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                            bufsize=1  # Line buffered
                        )
                        
                        # Monitor output in real-time
                        def monitor_output(stream, is_error=False):
                            nonlocal log_callback
                            prefix = "STDERR: " if is_error else "STDOUT: "
                            for line in iter(stream.readline, ''):
                                print(prefix + line.strip())
                                log.write(prefix + line)
                                log.flush()  # Force write to disk
                                if log_callback:
                                    try:
                                        log_callback(prefix + line)
                                    except Exception as e:
                                        # Handle case where the progress window has been closed
                                        if "invalid command name" in str(e) or "TclError" in str(e.__class__):
                                            print("Progress window closed, stopping callbacks")
                                            # Replace log_callback with a no-op function
                                            log_callback = lambda x: None
                                        else:
                                            # Re-raise unexpected errors
                                            print(f"Unexpected error in log callback: {e}")
                                            raise
                        
                        # Start threads to monitor stdout and stderr
                        stdout_thread = threading.Thread(target=monitor_output, args=(centrifuger_process.stdout,))
                        stderr_thread = threading.Thread(target=monitor_output, args=(centrifuger_process.stderr, True))
                        stdout_thread.daemon = True
                        stderr_thread.daemon = True
                        stdout_thread.start()
                        stderr_thread.start()
                        
                        # Wait for process to complete
                        return_code = centrifuger_process.wait()
                        stdout_thread.join()
                        stderr_thread.join()
                        
                        if return_code != 0:
                            raise subprocess.CalledProcessError(return_code, cmd)
                    else:
                        # For standard mode, use run directly
                        process = subprocess.run(
                            cmd,
                            check=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True
                        )
                        
                        # Log both stdout and stderr
                        log.write("=== STDOUT ===\n")
                        log.write(process.stdout)
                        log.write("\n=== STDERR ===\n")
                        log.write(process.stderr)
                        
                        # Also print to console
                        print("\nCentrifuger finished processing. Output summary:")
                        print("\nCentrifuger STDOUT:")
                        print(process.stdout)
                        print("\nCentrifuger STDERR:")
                        print(process.stderr)
                    
                    # After process completes
                    if realtime_mode:
                        print("Realtime mode has completed.")
                        self.stop_monitoring()  # Stop the monitoring thread
                    
                    # Show result file information
                    if os.path.exists(output_file):
                        file_size = os.path.getsize(output_file)
                        print(f"\nOutput file created: {output_file} ({file_size:,} bytes)")
                        
                        # Count classified reads
                        try:
                            with open(output_file) as f:
                                line_count = sum(1 for _ in f) - 1  # Subtract header line
                            print(f"Total classified reads: {line_count:,}")
                            
                            # If we reached this point, the analysis completed successfully
                            status_msg = f"Analysis completed successfully with {line_count:,} classified reads"
                            if log_callback:
                                log_callback(status_msg + "\n")
                            log.write(f"\n=== COMPLETED: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===\n")
                            log.write(f"{status_msg}\n")
                        except Exception as e:
                            print(f"Could not count reads in output file: {e}")
                    else:
                        print(f"Warning: Output file was not created: {output_file}")
                        
                except FileNotFoundError as e:
                    error_msg = f"Binary not found: {cmd[0]}. Make sure it's installed and in your PATH."
                    print(error_msg)
                    log.write(f"ERROR: {error_msg}\n")
                    raise e
                except subprocess.CalledProcessError as e:
                    error_msg = f"Centrifuger command failed with exit code {e.returncode}"
                    print(error_msg)
                    if hasattr(e, 'stdout') and e.stdout:
                        print(f"STDOUT: {e.stdout}")
                        log.write(f"STDOUT: {e.stdout}\n")
                    if hasattr(e, 'stderr') and e.stderr:
                        print(f"STDERR: {e.stderr}")
                        log.write(f"STDERR: {e.stderr}\n")
                    log.write(f"ERROR: {error_msg}\n")
                    raise e
            
            # Generate report using centrifuger-kreport and capture its output
            kreport_path = os.path.join(os.path.dirname(self.centrifuger_path), "centrifuger-kreport")
            if not os.path.exists(kreport_path):
                kreport_path = "centrifuger-kreport"  # Fall back to PATH
                
            print("\nGenerating Kraken-style report...")
            kreport_cmd = [
                kreport_path,
                "-x", full_db_path,
                output_file
            ]
            
            print(f"Kreport command: {' '.join(kreport_cmd)}")
            
            try:
                print("Running centrifuger-kreport...")
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
                
                # Print summary information about the report
                print(f"Kraken-style report generated: {report_file}")
                if os.path.exists(report_file):
                    file_size = os.path.getsize(report_file)
                    print(f"Report file size: {file_size:,} bytes")
                    
                    # Count lines in report
                    try:
                        with open(report_file) as f:
                            line_count = sum(1 for _ in f)
                        print(f"Report contains {line_count:,} taxonomic entries")
                    except Exception as e:
                        print(f"Could not analyze report file: {e}")
                        
                print("\nAnalysis Summary:")
                print(f"- Sample: {sample_name}")
                print(f"- Classification file: {output_file}")
                print(f"- Report file: {report_file}")
                print(f"- Log file: {log_file}")
                if realtime_mode:
                    print(f"- Watch directory: {watch_dir}")
                    print(f"- Min files threshold: {min_files}")
                    print("- Mode: Realtime monitoring (completed)")
                else:
                    print("- Mode: Standard (one-time analysis)")
                
                print("\nCentrifuger analysis completed successfully")
                return True
            except subprocess.CalledProcessError as e:
                print(f"Error running centrifuger-kreport: {e}")
                print(f"STDOUT: {e.stdout}")
                print(f"STDERR: {e.stderr}")
                return False
            
        except subprocess.CalledProcessError as e:
            print(f"Error running Centrifuger: {e}")
            print(f"STDOUT: {e.stdout if hasattr(e, 'stdout') else 'None'}")
            print(f"STDERR: {e.stderr if hasattr(e, 'stderr') else 'None'}")
            return False
        except FileNotFoundError as e:
            print(f"File not found error: {e}")
            return False
        except Exception as e:
            print(f"Unexpected error: {e}")
            return False
        finally:
            # Make sure to stop the monitoring thread if it was started
            if realtime_mode:
                self.stop_monitoring()

    def run_multi_sample(self, sample_sheet, db_path):
        """Run Centrifuger on multiple samples using a sample sheet"""
        cmd = [
            self.centrifuger_path,
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
        """Generate a Kraken-style report from Centrifuger output"""
        full_db_path = self.get_full_db_path(db_path)
        
        # Use the kreport tool from the same directory as the centrifuger binary
        kreport_path = os.path.join(os.path.dirname(self.centrifuger_path), "centrifuger-kreport")
        if not os.path.exists(kreport_path):
            kreport_path = "centrifuger-kreport"  # Fall back to PATH
            
        cmd = [
            kreport_path,
            "-x", full_db_path,
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
