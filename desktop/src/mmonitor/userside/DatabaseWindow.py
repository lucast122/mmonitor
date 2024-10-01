import os
import json
import subprocess
import threading
import shutil
import customtkinter as ctk
from tkinter import filedialog, messagebox
from build_mmonitor_pyinstaller import ROOT
import traceback
import datetime
import multiprocessing
import queue
import tempfile
import time
import gzip
import urllib.request
import tarfile
import ssl

def process_chunk(chunk):
    result = []
    skipped = 0
    for line in chunk:
        if line.startswith(">"):
            parts = line.split()
            seq_id = parts[0][1:]  # Remove the '>' character
            species_name = " ".join(parts[1:3])  # Use only genus and species
            taxid = get_taxid_from_species_name(species_name)
            if taxid:
                result.append(f"{seq_id}\t{taxid}\n")
            else:
                skipped += 1
                print(f"Skipped line (no taxid found): {line.strip()}")
        else:
            # This is a sequence line, not a header line
            continue
    print(f"Processed chunk: {len(chunk)} lines, {len(result)} results, {skipped} skipped")
    return result, skipped

def download_and_process_taxdump():
    taxdump_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    local_file = "taxdump.tar.gz"
    
    print("Downloading taxdump.tar.gz...")
    # Create an SSL context that doesn't verify certificates
    context = ssl._create_unverified_context()
    
    try:
        with urllib.request.urlopen(taxdump_url, context=context) as response, open(local_file, 'wb') as out_file:
            out_file.write(response.read())
    except Exception as e:
        print(f"Error downloading file: {e}")
        return {}

    print("Extracting taxdump.tar.gz...")
    with tarfile.open(local_file, "r:gz") as tar:
        tar.extractall()
    
    print("Processing names.dmp...")
    species_to_taxid = {}
    try:
        with open('names.dmp', 'r') as f:
            for line in f:
                fields = line.split('|')
                if len(fields) > 3 and fields[3].strip() == 'scientific name':
                    taxid = fields[0].strip()
                    name = fields[1].strip()
                    species_to_taxid[name] = taxid
    except Exception as e:
        print(f"Error processing names.dmp: {e}")
        return {}
    
    print("Saving processed mapping...")
    with open('species_taxid_mapping.json', 'w') as f:
        json.dump(species_to_taxid, f)
    
    print("Cleaning up...")
    os.remove(local_file)
    os.remove('names.dmp')
    if os.path.exists('nodes.dmp'):
        os.remove('nodes.dmp')
    if os.path.exists('delnodes.dmp'):
        os.remove('delnodes.dmp')
    if os.path.exists('merged.dmp'):
        os.remove('merged.dmp')
    
    return species_to_taxid

def get_taxid_from_species_name(species_name):
    mapping_file = 'species_taxid_mapping.json'
    
    if not os.path.exists(mapping_file):
        print("Mapping file not found. Downloading and processing...")
        species_taxid_map = download_and_process_taxdump()
        if not species_taxid_map:
            print("Failed to create species to taxid mapping.")
            return None
    else:
        with open(mapping_file, 'r') as f:
            species_taxid_map = json.load(f)
    
    taxid = species_taxid_map.get(species_name)
    
    if taxid:
        return taxid
    else:
        print(f"No taxid found for species: {species_name}")
        return None

def create_species_taxid_mapping(output_file):
    # This function should be called once to create the mapping file
    # You can implement the logic to create the mapping here
    # For example, you could use a local taxonomy database or a faster API
    # to create a dictionary of species names to taxids
    species_taxid_map = {}  # Populate this dictionary
    
    with open(output_file, 'w') as f:
        json.dump(species_taxid_map, f)
    
    print(f"Created species to taxid mapping file: {output_file}")

class DatabaseWindow(ctk.CTkToplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.title("Database Management")
        self.geometry("500x400")

        # Use a directory in the user's home folder or a temporary directory
        self.base_dir = os.path.join(os.path.expanduser("~"), ".mmonitor")
        if not os.path.exists(self.base_dir):
            os.makedirs(self.base_dir)

        self.emu_db_path = os.path.join(self.base_dir, "emu_db")
        self.centrifuge_db_path = os.path.join(self.base_dir, "centrifuge_db")
        self.config_file = os.path.join(self.base_dir, "db_config.json")
        self.log_file = os.path.join(self.base_dir, "database_build.log")

        self.load_config()
        self.create_widgets()

        self.message_queue = queue.Queue()
        self.after(100, self.process_message_queue)

    def create_widgets(self):
        frame = ctk.CTkFrame(self)
        frame.pack(padx=20, pady=20, fill="both", expand=True)

        ctk.CTkLabel(frame, text="Emu Database", font=("Helvetica", 16, "bold")).pack(pady=5)
        ctk.CTkButton(frame, text="Select Emu Database", command=self.select_emu_db).pack(pady=5)
        ctk.CTkButton(frame, text="Build Emu Database", command=self.confirm_build_emu_db).pack(pady=5)

        ctk.CTkLabel(frame, text="Centrifuge Database", font=("Helvetica", 16, "bold")).pack(pady=5)
        ctk.CTkButton(frame, text="Select Centrifuge Database", command=self.select_centrifuge_db).pack(pady=5)
        ctk.CTkButton(frame, text="Build Centrifuge Database", command=self.confirm_build_centrifuge_db).pack(pady=5)

    def load_config(self):
        if os.path.exists(self.config_file):
            with open(self.config_file, 'r') as f:
                config = json.load(f)
                self.emu_db_path = config.get('emu_db', self.emu_db_path)
                self.centrifuge_db_path = config.get('centrifuge_db', self.centrifuge_db_path)

    def save_config(self):
        config = {
            'emu_db': self.emu_db_path,
            'centrifuge_db': self.centrifuge_db_path
        }
        with open(self.config_file, 'w') as f:
            json.dump(config, f)

    def select_emu_db(self):
        db_path = filedialog.askdirectory(title="Select Emu Database Directory")
        if db_path:
            self.emu_db_path = db_path
            self.save_config()

    def select_centrifuge_db(self):
        db_path = filedialog.askdirectory(title="Select Centrifuge Database Directory")
        if db_path:
            self.centrifuge_db_path = db_path
            self.save_config()

    def confirm_build_emu_db(self):
        response = messagebox.askyesno("Confirm Build", "Building the Emu database requires an internet connection and may take a while depending on your internet speed and computer. Do you want to proceed?")
        if response:
            self.build_emu_db()

    def confirm_build_centrifuge_db(self):
        response = messagebox.askyesno("Confirm Build", "Building the Centrifuge database requires an internet connection and may take a while depending on your internet speed and computer. Do you want to proceed?")
        if response:
            self.build_centrifuge_db()

    def build_emu_db(self):
        self.show_progress_window("Building Emu Database")
        threading.Thread(target=self._build_emu_db).start()

    def build_centrifuge_db(self):
        self.show_progress_window("Building Centrifuge Database")
        threading.Thread(target=self._build_centrifuge_db).start()

    def show_progress_window(self, title):
        self.progress_window = ctk.CTkToplevel(self)
        self.progress_window.title(title)
        self.progress_window.geometry("300x100")
        self.progress_window.transient(self)
        self.progress_window.grab_set()

        ctk.CTkLabel(self.progress_window, text="Building database...").pack(pady=10)
        self.progress_bar = ctk.CTkProgressBar(self.progress_window, mode="indeterminate")
        self.progress_bar.pack(pady=10)
        self.progress_bar.start()

    def close_progress_window(self):
        if hasattr(self, 'progress_window'):
            self.progress_window.destroy()

    def _build_emu_db(self):
        try:
            # Create a new directory for the custom database
            custom_db_path = os.path.join(self.base_dir, "custom_emu_db")
            os.makedirs(custom_db_path, exist_ok=True)
            self.log_progress("Created custom database directory: " + custom_db_path)

            # Use a temporary directory for intermediate files
            with tempfile.TemporaryDirectory() as temp_dir:
                # Download NCBI taxonomy
                self.log_progress("Downloading NCBI taxonomy...")
                taxdump_path = os.path.join(temp_dir, "taxdump.tar.gz")
                result = subprocess.run(["wget", "-O", taxdump_path, "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"], 
                                        capture_output=True, text=True, check=True)
                self.log_progress(result.stdout)
                self.log_progress(result.stderr)
                
                # Extract taxonomy files
                self.log_progress("Extracting taxonomy files...")
                result = subprocess.run(["tar", "-xzvf", taxdump_path, "-C", custom_db_path], 
                                        capture_output=True, text=True, check=True)
                self.log_progress(result.stdout)
                self.log_progress(result.stderr)

                # Download 16S rRNA sequences
                self.log_progress("Downloading 16S rRNA sequences...")
                rna_path = os.path.join(temp_dir, "16S.fna.gz")
                result = subprocess.run(["wget", "-O", rna_path, "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"], 
                                        capture_output=True, text=True, check=True)
                self.log_progress(result.stdout)
                self.log_progress(result.stderr)
                
                result = subprocess.run(["gunzip", "-f", rna_path], 
                                        capture_output=True, text=True, check=True)
                self.log_progress(result.stdout)
                self.log_progress(result.stderr)

                # Create seq2taxid map
                self.log_progress("Creating seq2taxid map...")
                seq2taxid_path = os.path.join(custom_db_path, "seq2taxid.map")
                unzipped_rna_path = os.path.join(temp_dir, "16S.fna")
                skipped_entries, total_entries = self.create_seq2taxid_map(unzipped_rna_path, seq2taxid_path)

                self.log_progress(f"Seq2taxid map creation complete. Processed {total_entries} entries, skipped {skipped_entries} entries due to missing taxid.")

                if os.path.getsize(seq2taxid_path) == 0:
                    with open("16S.fna", "r") as f:
                        first_lines = f.readlines(1000)
                    self.log_progress(f"First 1000 lines of 16S.fna:\n{''.join(first_lines)}")
                    raise ValueError(f"seq2taxid map is empty. Cannot proceed with database building. Input file size: {os.path.getsize('16S.fna')} bytes, Output file size: {os.path.getsize(seq2taxid_path)} bytes")

                # Build Emu database
                self.log_progress("Building Emu database...")
                emu_output_path = os.path.join(custom_db_path, "emu_custom_db")
                if os.path.exists(emu_output_path):
                    shutil.rmtree(emu_output_path)
                
                # Get the number of available CPU cores
                num_cores = multiprocessing.cpu_count()
                self.log_progress(f"Using {num_cores} CPU cores for database building")
                emu_path = os.path.join(ROOT, "lib", "emu.py")
                emu_build_command = [
                    sys.executable, emu_path, "build-database", "emu_custom_db",
                    "--sequences", "16S.fna",
                    "--seq2tax", seq2taxid_path,
                    "--ncbi-taxonomy", custom_db_path,
                    "--output", emu_output_path,
                    "--threads", str(num_cores)  # Use all available cores
                ]
                self.log_progress(f"Running command: {' '.join(emu_build_command)}")
                
                process = subprocess.Popen(emu_build_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                
                while True:
                    output = process.stdout.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        self.log_progress(output.strip())
                
                rc = process.poll()
                if rc != 0:
                    error_output = process.stderr.read()
                    raise subprocess.CalledProcessError(rc, emu_build_command, error_output)

                self.log_progress("Emu database build completed.")

                # Clean up
                self.log_progress("Cleaning up temporary files...")
                os.remove("taxdump.tar.gz")
                os.remove("16S.fna")

                self.close_progress_window()
                messagebox.showinfo("Success", f"Emu database built successfully at {emu_output_path}")
                
                # Update the emu_db_path to the new custom database
                self.emu_db_path = emu_output_path
                self.save_config()
                
        except subprocess.CalledProcessError as e:
            self.close_progress_window()
            error_message = f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.\n"
            if e.stdout:
                error_message += f"Stdout: {e.stdout}\n"
            if e.stderr:
                error_message += f"Stderr: {e.stderr}\n"
            self.log_progress(f"Error: {error_message}")
            self.show_error_message(error_message)
        except Exception as e:
            self.close_progress_window()
            error_message = f"An unexpected error occurred: {str(e)}\n\nTraceback:\n{traceback.format_exc()}"
            self.log_progress(f"Error: {error_message}")
            self.show_error_message(error_message)
        finally:
            # Ensure cleanup happens even if an error occurs
            if os.path.exists("taxdump.tar.gz"):
                os.remove("taxdump.tar.gz")
            if os.path.exists("16S.fna"):
                os.remove("16S.fna")

    def create_seq2taxid_map(self, input_file, output_file):
        try:
            if not os.path.exists(input_file):
                raise FileNotFoundError(f"Input file not found: {input_file}")

            file_size = os.path.getsize(input_file)
            self.log_progress(f"Input file size: {file_size} bytes")

            with open(input_file, 'r') as infile:
                # Read and log the first few lines of the file
                first_lines = infile.readlines(20)
                self.log_progress(f"First 20 lines of input file:\n{''.join(first_lines)}")
                infile.seek(0)  # Reset file pointer to the beginning

                chunks = []
                current_chunk = []
                line_count = 0
                for line in infile:
                    line_count += 1
                    if line.startswith(">"):
                        if current_chunk:
                            chunks.append(current_chunk)
                            current_chunk = []
                    current_chunk.append(line)
                    
                    if line_count % 100000 == 0:
                        self.log_progress(f"Processed {line_count} lines...")

                if current_chunk:
                    chunks.append(current_chunk)

            self.log_progress(f"Total lines read: {line_count}")
            self.log_progress(f"Total chunks created: {len(chunks)}")

            if not chunks:
                raise ValueError("No chunks created from input file. File might be empty or have an unexpected format.")

            self.log_progress("Starting multiprocessing pool...")
            with multiprocessing.Pool() as pool:
                results = pool.map(process_chunk, chunks)
            self.log_progress("Multiprocessing pool completed.")

            total_processed = 0
            total_skipped = 0

            self.log_progress("Writing results to output file...")
            with open(output_file, 'w') as outfile:
                for result, skipped in results:
                    outfile.writelines(result)
                    total_processed += len(result)
                    total_skipped += skipped

            self.log_progress(f"Total processed: {total_processed} entries")
            self.log_progress(f"Total skipped: {total_skipped} entries")
            self.log_progress(f"Output file size: {os.path.getsize(output_file)} bytes")

            if total_processed == 0:
                raise ValueError(f"No entries processed in input file. File might be empty or have an unexpected format.")

            return total_skipped, total_processed
        except Exception as e:
            self.log_progress(f"Error in create_seq2taxid_map: {str(e)}")
            if os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    self.log_progress(f"First 1000 characters of output file:\n{f.read(1000)}")
            raise

    def show_error_message(self, message):
        def show_message():
            messagebox.showerror("Error", message)
        self.after(0, show_message)

    def log_progress(self, message):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_message = f"[{timestamp}] {message}"
        print(log_message)  # Print to console
        self.update_progress_window(log_message)
        
        # Write to log file
        with open(self.log_file, "a") as f:
            f.write(log_message + "\n")

    def update_progress_window(self, message):
        try:
            self.message_queue.put(message)
        except tk.TclError:
            print(f"GUI update failed: {message}")

    def process_message_queue(self):
        try:
            while True:
                message = self.message_queue.get_nowait()
                if hasattr(self, 'progress_window') and hasattr(self, 'progress_text'):
                    self.progress_text.insert("end", message + "\n")
                    self.progress_text.see("end")
        except queue.Empty:
            pass
        finally:
            self.after(100, self.process_message_queue)

    def _build_centrifuge_db(self):
        try:
            index_name = "centrifuge_custom_index"
            taxonomy_dir = os.path.join(self.centrifuge_db_path, "taxonomy")
            library_dir = os.path.join(self.centrifuge_db_path, "library")

            # Remove existing directories if they exist
            if os.path.exists(taxonomy_dir):
                shutil.rmtree(taxonomy_dir)
            if os.path.exists(library_dir):
                shutil.rmtree(library_dir)

            os.makedirs(taxonomy_dir)
            os.makedirs(library_dir)

            # Download NCBI taxonomy
            subprocess.run(["centrifuge-download", "-o", taxonomy_dir, "taxonomy"], check=True)

            # Download genomes
            with open("seqid2taxid.map", "w") as f:
                subprocess.run(["centrifuge-download", "-o", library_dir, "-m", "-d", "archaea,bacteria,fungi", "refseq"], stdout=f, check=True)

            # Concatenate sequences
            with open("input-sequences.fna", "w") as outfile:
                for root, dirs, files in os.walk(library_dir):
                    for file in files:
                        if file.endswith(".fna"):
                            with open(os.path.join(root, file), "r") as infile:
                                outfile.write(infile.read())

            # Build Centrifuge index
            subprocess.run([
                "centrifuge-build",
                "-p", str(os.cpu_count()),
                "--conversion-table", "seqid2taxid.map",
                "--taxonomy-tree", os.path.join(taxonomy_dir, "nodes.dmp"),
                "--name-table", os.path.join(taxonomy_dir, "names.dmp"),
                "input-sequences.fna",
                index_name
            ], check=True)

            # Clean up
            os.remove("input-sequences.fna")
            os.remove("seqid2taxid.map")

            self.close_progress_window()
            messagebox.showinfo("Success", "Centrifuge database built successfully")
        except subprocess.CalledProcessError as e:
            self.close_progress_window()
            messagebox.showerror("Error", f"Failed to build Centrifuge database: {str(e)}")
        except Exception as e:
            self.close_progress_window()
            messagebox.showerror("Error", f"An unexpected error occurred: {str(e)}")