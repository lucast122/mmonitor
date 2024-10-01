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
import sys
import requests

def get_taxid_from_species_name(species_name):
    # First, try to get the taxid from our existing mapping
    taxid = species_taxid_map.get(species_name)
    if taxid:
        return taxid
    
    # If not found, try to fetch from NCBI
    # try:
    #     print(f"Fetching taxid for {species_name} from NCBI...")
    #     url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={species_name}&retmode=json"
    #     response = requests.get(url)
    #     data = response.json()
    #     id_list = data['esearchresult']['idlist']
    #     if id_list:
    #         return id_list[0]
    # except Exception as e:
    #     print(f"Error fetching taxid for {species_name}: {e}")
    else:
        return None

def download_and_process_taxdump():
    taxdump_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    local_file = "taxdump.tar.gz"
    
    print("Downloading taxdump.tar.gz...")
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

def process_chunk(chunk, species_taxid_map):
    result = []
    skipped = 0
    missing_keys = set()
    for line in chunk:
        if line.startswith(">"):
            parts = line.split()
            seq_id = parts[0][1:]  # Remove the '>' character
            species_name = " ".join(parts[1:3])  # Use only genus and species
            taxid = species_taxid_map.get(species_name)
            if taxid:
                result.append(f"{seq_id}\t{taxid}\n")
            else:
                skipped += 1
                missing_keys.add(species_name)
    print(f"Processed chunk: {len(chunk)} lines, {len(result)} results, {skipped} skipped")
    if missing_keys:
        print(f"Missing keys: {', '.join(missing_keys)}")
    return result, skipped, missing_keys

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
        selected_domains = self.create_domain_selection_dialog()
        if selected_domains:
            response = messagebox.askyesno("Confirm Build", f"Building the Emu database for {', '.join(selected_domains)} requires an internet connection and may take a while. MMonitor will become unresponsive during this process. Do you want to proceed?")
            if response:
                self.build_emu_db(selected_domains)

    def confirm_build_centrifuge_db(self):
        selected_domains = self.create_domain_selection_dialog()
        if selected_domains:
            response = messagebox.askyesno("Confirm Build", f"Building the Centrifuge database for {', '.join(selected_domains)} requires an internet connection and may take a while. MMonitor will become unresponsive during this process. Do you want to proceed?")
            if response:
                self.build_centrifuge_db(selected_domains)

    def build_emu_db(self, selected_domains):
        self.show_progress_window("Building Emu Database")
        threading.Thread(target=self._build_emu_db, args=(selected_domains,)).start()

    def build_centrifuge_db(self, selected_domains):
        self.show_progress_window("Building Centrifuge Database")
        threading.Thread(target=self._build_centrifuge_db, args=(selected_domains,)).start()

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

    def _build_emu_db(self, selected_domains):
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

                # Download 16S rRNA sequences for selected domains
                self.log_progress("Downloading 16S rRNA sequences...")
                rna_files = []
                for domain in selected_domains:
                    if domain == "bacteria":
                        url = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
                    elif domain == "archaea":
                        url = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz"
                    elif domain == "fungi":
                        url = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.18SrRNA.fna.gz"
                    
                    rna_path = os.path.join(temp_dir, f"{domain}_rRNA.fna.gz")
                    result = subprocess.run(["wget", "-O", rna_path, url], 
                                            capture_output=True, text=True, check=True)
                    self.log_progress(result.stdout)
                    self.log_progress(result.stderr)
                    
                    result = subprocess.run(["gunzip", "-f", rna_path], 
                                            capture_output=True, text=True, check=True)
                    self.log_progress(result.stdout)
                    self.log_progress(result.stderr)
                    
                    rna_files.append(rna_path[:-3])  # Remove .gz extension

                # Concatenate all RNA files
                combined_rna_path = os.path.join(temp_dir, "combined_rRNA.fna")
                with open(combined_rna_path, 'w') as outfile:
                    for rna_file in rna_files:
                        with open(rna_file) as infile:
                            outfile.write(infile.read())

                # Create seq2taxid map
                self.log_progress("Creating seq2taxid map...")
                seq2taxid_path = os.path.join(custom_db_path, "seq2taxid.map")
                skipped_entries, total_entries = self.create_seq2taxid_map(combined_rna_path, seq2taxid_path)

                self.log_progress(f"Seq2taxid map creation complete. Processed {total_entries} entries, skipped {skipped_entries} entries due to missing taxid.")

                if os.path.getsize(seq2taxid_path) == 0:
                    raise ValueError(f"seq2taxid map is empty. Cannot proceed with database building.")

                # Build Emu database
                self.log_progress("Building Emu database...")
                emu_output_path = os.path.join(ROOT, "src", "resources", "emu_custom_db")
                if os.path.exists(emu_output_path):
                    shutil.rmtree(emu_output_path)
                
                num_cores = multiprocessing.cpu_count()
                self.log_progress(f"Preprocessing finished. Starting database building...")

                emu_script_path = os.path.join(ROOT, "lib", "emu.py")
                emu_build_command = [
                    sys.executable,
                    emu_script_path,
                    "build-database",
                    "--sequences", combined_rna_path,
                    "--seq2tax", seq2taxid_path,
                    "--ncbi-taxonomy", custom_db_path,
                    emu_output_path
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

            self.close_progress_window()
            messagebox.showinfo("Success", f"Emu database built successfully at {emu_output_path}")
            
            # Update the emu_db_path to the new custom database
            self.emu_db_path = emu_output_path
            self.save_config()

            # Ask user if they want to open the folder
            self.ask_to_open_folder(custom_db_path)
            
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

    def create_seq2taxid_map(self, input_file, output_file):
        try:
            if not os.path.exists(input_file):
                raise FileNotFoundError(f"Input file not found: {input_file}")

            file_size = os.path.getsize(input_file)
            self.log_progress(f"Input file size: {file_size} bytes")

            self.log_progress("Downloading and processing taxonomy data...")
            species_taxid_map = download_and_process_taxdump()

            self.log_progress("Processing input file and creating seq2taxid map...")
            total_processed = 0
            total_skipped = 0
            all_missing_keys = set()

            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(input_file + '.filtered', 'w') as filtered_file:
                current_seq = []
                include_current_seq = False
                for line in infile:
                    if line.startswith('>'):
                        if current_seq:
                            if include_current_seq:
                                filtered_file.writelines(current_seq)
                            current_seq = []
                        
                        parts = line.split()
                        seq_id = parts[0][1:]
                        species_name = " ".join(parts[1:3])
                        taxid = species_taxid_map.get(species_name)
                        
                        if taxid:
                            outfile.write(f"{seq_id}\t{taxid}\n")
                            total_processed += 1
                            include_current_seq = True
                        else:
                            total_skipped += 1
                            all_missing_keys.add(species_name)
                            include_current_seq = False
                    
                    current_seq.append(line)
                
                if current_seq and include_current_seq:
                    filtered_file.writelines(current_seq)

            # Replace the original input file with the filtered one
            os.replace(input_file + '.filtered', input_file)

            self.log_progress(f"Total processed: {total_processed} entries")
            self.log_progress(f"Total skipped: {total_skipped} entries")
            self.log_progress(f"Output file size: {os.path.getsize(output_file)} bytes")
            self.log_progress(f"Total missing keys: {len(all_missing_keys)}")
            if all_missing_keys:
                self.log_progress(f"Sample of missing keys: {', '.join(list(all_missing_keys)[:10])}")

            # Write error details to a separate file
            error_file = os.path.join(os.path.dirname(output_file), "seq2taxid_errors.log")
            with open(error_file, 'w') as f:
                for species in all_missing_keys:
                    f.write(f"Missing taxid for species: {species}\n")
            self.log_progress(f"Detailed error log written to: {error_file}")

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

    def _build_centrifuge_db(self, selected_domains):
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

            # Download genomes for selected domains
            with open("seqid2taxid.map", "w") as f:
                subprocess.run(["centrifuge-download", "-o", library_dir, "-m", "-d", ",".join(selected_domains), "refseq"], stdout=f, check=True)

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
            messagebox.showinfo("Success", f"Centrifuge database built successfully at {self.centrifuge_db_path}")

            # Ask user if they want to open the folder
            self.ask_to_open_folder(self.centrifuge_db_path)

        except subprocess.CalledProcessError as e:
            self.close_progress_window()
            messagebox.showerror("Error", f"Failed to build Centrifuge database: {str(e)}")
        except Exception as e:
            self.close_progress_window()
            messagebox.showerror("Error", f"An unexpected error occurred: {str(e)}")

    def create_domain_selection_dialog(self):
        dialog = ctk.CTkToplevel(self)
        dialog.title("Select Domains")
        dialog.geometry("300x200")

        var_bacteria = ctk.BooleanVar(value=True)
        var_archaea = ctk.BooleanVar(value=True)
        var_fungi = ctk.BooleanVar(value=True)

        ctk.CTkCheckBox(dialog, text="Bacteria", variable=var_bacteria).pack(pady=5)
        ctk.CTkCheckBox(dialog, text="Archaea", variable=var_archaea).pack(pady=5)
        ctk.CTkCheckBox(dialog, text="Fungi", variable=var_fungi).pack(pady=5)

        selected_domains = []

        def on_confirm():
            if var_bacteria.get():
                selected_domains.append("bacteria")
            if var_archaea.get():
                selected_domains.append("archaea")
            if var_fungi.get():
                selected_domains.append("fungi")
            dialog.destroy()

        ctk.CTkButton(dialog, text="Confirm", command=on_confirm).pack(pady=20)

        self.wait_window(dialog)
        return selected_domains

    def ask_to_open_folder(self, folder_path):
        response = messagebox.askyesno("Open Folder", f"Do you want to open the folder containing the built database?\n\nPath: {folder_path}")
        if response:
            try:
                if sys.platform == "win32":
                    os.startfile(folder_path)
                elif sys.platform == "darwin":
                    subprocess.Popen(["open", folder_path])
                else:  # assume linux or unix
                    subprocess.Popen(["xdg-open", folder_path])
            except Exception as e:
                messagebox.showerror("Error", f"Failed to open folder: {str(e)}")