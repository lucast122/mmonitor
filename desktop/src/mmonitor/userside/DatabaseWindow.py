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
import wget
import ftplib
from io import BytesIO

# Add the parent directory of 'mmonitor' to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Remove or comment out this line:
# from mmonitor.userside.DatabaseWindow import DatabaseWindow

# Check if process_chunk is defined here
def process_chunk(chunk):
    # Implement the function
    pass

class DatabaseWindow(ctk.CTkFrame):
    def __init__(self, parent):
        print("Initializing DatabaseWindow")
        super().__init__(parent)
        self.parent = parent

        # Use a directory in the user's home folder or a temporary directory
        self.base_dir = os.path.join(os.path.expanduser("~"), ".mmonitor")
        if not os.path.exists(self.base_dir):
            os.makedirs(self.base_dir)

        self.emu_db_path = os.path.join(self.base_dir, "emu_db")
        self.centrifuge_db_path = os.path.join(self.base_dir, "centrifuge_db")
        self.gtdb_db_path = os.path.join(self.base_dir, "gtdb_db")
        self.config_file = os.path.join(self.base_dir, "db_config.json")
        self.log_file = os.path.join(self.base_dir, "database_build.log")

        self.load_config()
        self.create_widgets()
        print("DatabaseWindow initialized")

        self.message_queue = queue.Queue()
        self.after(100, self.process_message_queue)
        self.species_taxid_map = {}  # Initialize the map

    def create_widgets(self):
        print("Creating DatabaseWindow widgets")
        
        # Main frame
        main_frame = ctk.CTkFrame(self)
        main_frame.pack(fill="both", expand=True, padx=20, pady=20)

        # Title and explanation
        title_label = ctk.CTkLabel(main_frame, text="Database Management", font=("Helvetica", 24, "bold"))
        title_label.pack(pady=(0, 10))
        
        explanation = ("Here you can manage the databases used for taxonomic classification.\n"
                       "To select an existing custom database chose either Emu or Centrifuge tab and click 'Browse'."
                       "To build a new database, select the reference domains and click 'Build Database'."
                       "This will download the required files and create a custom database in the selected directory."
                       "After the database is built, you can select it in the 'Database Path' field and click 'Confirm'."
                       "The database will be used for taxonomic classification in the analysis pipeline."
                       "")
        ctk.CTkLabel(main_frame, text=explanation, wraplength=500).pack(pady=(0, 20))

        # Create a notebook (tabbed interface)
        self.notebook = ctk.CTkTabview(main_frame)
        self.notebook.pack(fill="both", expand=True)

        # Emu Database Tab
        emu_tab = self.notebook.add("Emu Database (16S)")
        self.create_emu_tab(emu_tab)

        # Centrifuge Database Tab
        centrifuge_tab = self.notebook.add("Centrifuge Database (WGS)")
        self.create_centrifuge_tab(centrifuge_tab)

        print("DatabaseWindow widgets created")

    def create_emu_tab(self, parent):
        frame = ctk.CTkFrame(parent)
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        ctk.CTkLabel(frame, text="Emu Database Path:", anchor="w").pack(fill="x", pady=(0, 5))
        path_frame = ctk.CTkFrame(frame)
        path_frame.pack(fill="x")

        self.emu_path_entry = ctk.CTkEntry(path_frame)
        self.emu_path_entry.pack(side="left", fill="x", expand=True)
        self.emu_path_entry.insert(0, self.emu_db_path)

        ctk.CTkButton(path_frame, text="Browse", command=self.select_emu_db).pack(side="right", padx=(5, 0))

        domains_frame = ctk.CTkFrame(frame)
        domains_frame.pack(fill="x", pady=(20, 0))

        ctk.CTkLabel(domains_frame, text="Select Domains:", anchor="w").pack(fill="x", pady=(0, 5))

        self.emu_domains = {
            "Bacteria": ctk.BooleanVar(value=True),
            "Archaea": ctk.BooleanVar(value=True),
            "Fungi": ctk.BooleanVar(value=False)
        }

        for domain, var in self.emu_domains.items():
            ctk.CTkCheckBox(domains_frame, text=domain, variable=var).pack(anchor="w", pady=2)

        ctk.CTkButton(frame, text="Build Emu Database", command=self.confirm_build_emu_db).pack(fill="x", pady=(20, 0))

    def create_centrifuge_tab(self, parent):
        frame = ctk.CTkFrame(parent)
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        ctk.CTkLabel(frame, text="Centrifuge Database Path:", anchor="w").pack(fill="x", pady=(0, 5))
        path_frame = ctk.CTkFrame(frame)
        path_frame.pack(fill="x")

        self.centrifuge_path_entry = ctk.CTkEntry(path_frame)
        self.centrifuge_path_entry.pack(side="left", fill="x", expand=True)
        self.centrifuge_path_entry.insert(0, self.centrifuge_db_path)

        ctk.CTkButton(path_frame, text="Browse", command=self.select_centrifuge_db).pack(side="right", padx=(5, 0))

        # Database Type Selection
        db_type_frame = ctk.CTkFrame(frame)
        db_type_frame.pack(fill="x", pady=(20, 0))

        ctk.CTkLabel(db_type_frame, text="Select Database Type:", anchor="w").pack(fill="x", pady=(0, 5))

        self.centrifuge_db_type = ctk.StringVar(value="NCBI")
        ctk.CTkRadioButton(db_type_frame, text="NCBI", variable=self.centrifuge_db_type, value="NCBI").pack(anchor="w", pady=2)
        ctk.CTkRadioButton(db_type_frame, text="GTDB", variable=self.centrifuge_db_type, value="GTDB").pack(anchor="w", pady=2)

        # Index File/Directory Selection
        index_frame = ctk.CTkFrame(frame)
        index_frame.pack(fill="x", pady=(20, 0))

        ctk.CTkLabel(index_frame, text="Index File/Directory:", anchor="w").pack(fill="x", pady=(0, 5))
        self.index_entry = ctk.CTkEntry(index_frame)
        self.index_entry.pack(side="left", fill="x", expand=True)
        ctk.CTkButton(index_frame, text="Browse", command=self.select_index).pack(side="right", padx=(5, 0))

        domains_frame = ctk.CTkFrame(frame)
        domains_frame.pack(fill="x", pady=(20, 0))

        ctk.CTkLabel(domains_frame, text="Select Domains:", anchor="w").pack(fill="x", pady=(0, 5))

        self.centrifuge_domains = {
            "Bacteria": ctk.BooleanVar(value=True),
            "Archaea": ctk.BooleanVar(value=True),
            "Viruses": ctk.BooleanVar(value=False),
            "Human": ctk.BooleanVar(value=False)
        }

        for domain, var in self.centrifuge_domains.items():
            ctk.CTkCheckBox(domains_frame, text=domain, variable=var).pack(anchor="w", pady=2)

        ctk.CTkButton(frame, text="Build Centrifuge Database", command=self.confirm_build_centrifuge_db).pack(fill="x", pady=(20, 0))

    def select_emu_db(self):
        db_path = filedialog.askdirectory(title="Select Emu Database Directory")
        if db_path:
            self.emu_db_path = db_path
            self.emu_path_entry.delete(0, ctk.END)
            self.emu_path_entry.insert(0, db_path)
            self.save_config()

    def select_centrifuge_db(self):
        db_path = filedialog.askdirectory(title="Select Centrifuge Database Directory")
        if db_path:
            self.centrifuge_db_path = db_path
            self.centrifuge_path_entry.delete(0, ctk.END)
            self.centrifuge_path_entry.insert(0, db_path)
            self.save_config()

    def select_index(self):
        index_path = filedialog.askopenfilename(title="Select Index File")
        if index_path:
            self.index_entry.delete(0, ctk.END)
            self.index_entry.insert(0, index_path)

    def confirm_build_emu_db(self):
        selected_domains = [domain for domain, var in self.emu_domains.items() if var.get()]
        if selected_domains:
            response = messagebox.askyesno("Confirm Build", f"Building the Emu database for {', '.join(selected_domains)} requires an internet connection and may take a while. Do you want to proceed?")
            if response:
                self.build_emu_db(selected_domains)

    def confirm_build_centrifuge_db(self):
        selected_domains = [domain for domain, var in self.centrifuge_domains.items() if var.get()]
        if selected_domains:
            db_type = self.centrifuge_db_type.get()
            index_path = self.index_entry.get()
            response = messagebox.askyesno("Confirm Build", f"Building the Centrifuge database ({db_type}) for {', '.join(selected_domains)} requires an internet connection and may take a while. Do you want to proceed?")
            if response:
                self.build_centrifuge_db(selected_domains, db_type, index_path)

    def build_emu_db(self, selected_domains):
        self.show_progress_window("Building Emu Database")
        threading.Thread(target=self._build_emu_db, args=(selected_domains,)).start()

    def build_centrifuge_db(self, selected_domains, db_type, index_path):
        self.show_progress_window("Building Centrifuge Database")
        threading.Thread(target=self._build_centrifuge_db, args=(selected_domains, db_type, index_path)).start()

    def show_progress_window(self, title):
        self.progress_window = ctk.CTkToplevel(self)
        self.progress_window.title(title)
        self.progress_window.geometry("400x150")
        self.progress_window.transient(self)
        self.progress_window.grab_set()

        self.progress_label = ctk.CTkLabel(self.progress_window, text="Initializing...")
        self.progress_label.pack(pady=10)

        self.progress_bar = ctk.CTkProgressBar(self.progress_window, mode="indeterminate")
        self.progress_bar.pack(pady=10, padx=20, fill="x")
        self.progress_bar.start()

    def close_progress_window(self):
        if hasattr(self, 'progress_window'):
            self.progress_window.destroy()

    def log_progress(self, message):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_message = f"[{timestamp}] {message}"
        print(log_message)  # Print to console
        self.update_progress_window(log_message)
        
        # Write to log file
        with open(self.log_file, "a") as f:
            f.write(log_message + "\n")
        
        # Update log display in the GUI
        self.log_text.insert(ctk.END, log_message + "\n")
        self.log_text.see(ctk.END)

    def update_progress_window(self, message):
        if hasattr(self, 'progress_label'):
            self.progress_label.configure(text=message)

    def show_error_message(self, message):
        def show_message():
            messagebox.showerror("Error", message)
        self.after(0, show_message)

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
                seq2taxid_path = os.path.join(temp_dir, "seq2taxid.map")
                
                # Ensure the combined RNA file is created
                with open(combined_rna_path, 'w') as f:
                    f.write("Dummy content for testing")
                
                skipped_entries, total_entries = self.create_seq2taxid_map(combined_rna_path, seq2taxid_path)
                
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
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        try:
            raise Exception("Test exception")
            
            
            return skipped, total
        except Exception as e:
            self.log_progress(f"Error in create_seq2taxid_map: {str(e)}")
            raise

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

    def _build_centrifuge_db(self, selected_domains, db_type, index_path):
        try:
            self.log_progress(f"Creating Centrifuge database directory...")
            os.makedirs(self.centrifuge_db_path, exist_ok=True)

            with tempfile.TemporaryDirectory() as temp_dir:
                if db_type == "NCBI":
                    # Build Centrifuge database using NCBI data
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

                    # Set paths for centrifuge executables
                    centrifuge_download = os.path.join(ROOT,  "lib", "centrifuge_mac", "centrifuge-download")
                    centrifuge_build = os.path.join(ROOT, "lib", "centrifuge_mac", "centrifuge-build")

                    # Download NCBI taxonomy
                    self.update_progress("Downloading NCBI taxonomy...")
                    subprocess.run([centrifuge_download, "-o", taxonomy_dir, "taxonomy"], check=True)

                    # Download genomes for selected domains
                    domains = ",".join(selected_domains)
                    self.update_progress(f"Downloading {domains} genomes... This may take a while.")
                    subprocess.run([centrifuge_download, "-o", library_dir, "-m", "-d", domains, "refseq"], 
                                check=True, stdout=subprocess.PIPE)

                    # Concatenate all downloaded sequences
                    self.update_progress("Concatenating all downloaded sequences...")
                    input_sequences = os.path.join(self.centrifuge_db_path, "input-sequences.fna")
                    seqid2taxid_map = os.path.join(self.centrifuge_db_path, "seqid2taxid.map")

                    try:
                        with open(input_sequences, 'wb') as outfile:
                            for root, dirs, files in os.walk(library_dir):
                                for file in files:
                                    if file.endswith('.fna'):
                                        file_path = os.path.join(root, file)
                                        with open(file_path, 'rb') as infile:
                                            shutil.copyfileobj(infile, outfile)
                        
                        self.update_progress(f"Concatenation complete. File size: {os.path.getsize(input_sequences)} bytes")

                        if os.path.getsize(input_sequences) == 0:
                            raise ValueError(f"Concatenated file is empty: {input_sequences}")

                    except Exception as e:
                        raise RuntimeError(f"Error during sequence concatenation: {str(e)}")

                    # Check if files were created and have content
                    if not os.path.exists(input_sequences) or os.path.getsize(input_sequences) == 0:
                        raise ValueError(f"Input sequences file is empty or does not exist: {input_sequences}")
                    if not os.path.exists(seqid2taxid_map) or os.path.getsize(seqid2taxid_map) == 0:
                        raise ValueError(f"Seqid2taxid map file is empty or does not exist: {seqid2taxid_map}")

                    self.update_progress("Files created successfully. Starting Centrifuge index build...")

                    # Build Centrifuge index
                    self.update_progress("Building Centrifuge index... This may take a while.")
                    subprocess.run([
                        centrifuge_build,
                        "-p", str(os.cpu_count()),
                        "--conversion-table", seqid2taxid_map,
                        "--taxonomy-tree", os.path.join(taxonomy_dir, "nodes.dmp"),
                        "--name-table", os.path.join(taxonomy_dir, "names.dmp"),
                        input_sequences,
                        index_name
                    ], check=True)

                elif db_type == "GTDB":
                    # Build Centrifuge database using GTDB data
                    self.log_progress("Creating GTDB database directory...")
                    os.makedirs(self.gtdb_db_path, exist_ok=True)

                    # Step 1: Download the latest GTDB release
                    self.log_progress("Downloading GTDB metadata...")
                    gtdb_metadata_url = "https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_representative/"
                    gtdb_metadata_file = os.path.join(temp_dir, "gtdb_metadata.tar.gz")
                    wget.download(gtdb_metadata_url, gtdb_metadata_file)

                    # Step 2: Extract GTDB metadata
                    self.log_progress("Extracting GTDB metadata...")
                    with tarfile.open(gtdb_metadata_file, "r:gz") as tar:
                        tar.extractall(path=self.gtdb_db_path)

                    # Step 3: Prepare GTDB sequences for Centrifuger
                    self.log_progress("Preparing GTDB sequences for Centrifuger...")
                    gtdb_sequences_file = os.path.join(self.gtdb_db_path, "gtdb_sequences.fna")
                    self._prepare_gtdb_sequences(gtdb_sequences_file)

                    # Step 4: Build GTDB index using Centrifuger
                    self.log_progress("Building GTDB index using Centrifuger...")
                    centrifuger_build_command = [
                        "centrifuger-build",
                        "-p", str(multiprocessing.cpu_count()),
                        "--conversion-table", os.path.join(self.gtdb_db_path, "seqid2taxid.map"),
                        gtdb_sequences_file,
                        os.path.join(self.gtdb_db_path, "gtdb_index")
                    ]
                    subprocess.run(centrifuger_build_command, check=True)

                    self.log_progress("GTDB database build completed successfully.")

            self.close_progress_window()
            messagebox.showinfo("Success", f"Centrifuge database built successfully at {self.centrifuge_db_path}")

            # Ask user if they want to open the folder
            self.ask_to_open_folder(self.centrifuge_db_path)

        except subprocess.CalledProcessError as e:
            self.close_progress_window()
            error_message = f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.\nStderr: {e.stderr}"
            self.log_progress(f"Error: {error_message}")
            self.show_error_message(error_message)
        except Exception as e:
            self.close_progress_window()
            error_message = f"An unexpected error occurred: {str(e)}\n\nTraceback:\n{traceback.format_exc()}"
            self.log_progress(f"Error: {error_message}")
            self.show_error_message(error_message)

    def _prepare_gtdb_sequences(self, output_file):
        self.log_progress("Generating GTDB sequences from metadata...")
        # Placeholder logic for preparing GTDB sequences
        with open(output_file, 'w') as f:
            f.write(">Sample_Genome_1\nATCGATCGATCG\n")  # Dummy data for demonstration

    def process_message_queue(self):
        try:
            while True:
                message = self.message_queue.get_nowait()
                self.log_text.insert(ctk.END, message + "\n")
                self.log_text.see(ctk.END)
        except queue.Empty:
            pass
        finally:
            self.after(100, self.process_message_queue)