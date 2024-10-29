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
from mmonitor.userside.utils import create_tooltip
import psutil
import logging

# Remove this import
# from mmonitor.userside.PipelineConfig import PipelineConfig

class DatabaseWindow(ctk.CTkFrame):
    def __init__(self, parent):
        print("Initializing DatabaseWindow")
        super().__init__(parent)
        self.parent = parent

        # Initialize paths in resources directory
        self.resources_dir = os.path.join(ROOT, "src", "resources")
        self.emu_db_path = os.path.join(self.resources_dir, "custom_emu_db")
        self.centrifuger_db_path = os.path.join(self.resources_dir, "custom_centrifuger_db")
        self.gtdb_db_path = os.path.join(self.resources_dir, "custom_gtdb_db")
        self.config_file = os.path.join(self.resources_dir, "pipeline_config.json")  # Changed from db_config.json
        self.log_file = os.path.join(self.resources_dir, "database_build.log")

        # Create necessary directories
        for path in [self.resources_dir, self.emu_db_path, self.centrifuger_db_path, self.gtdb_db_path]:
            os.makedirs(path, exist_ok=True)

        # Configure logging
        self.logger = logging.getLogger('DatabaseWindow')
        self.logger.setLevel(logging.INFO)
        
        # File handler
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        self.logger.addHandler(file_handler)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter('%(message)s'))
        self.logger.addHandler(console_handler)

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
                       "You can build and update both Emu (16S) and Centrifuge (WGS) databases.")
        ctk.CTkLabel(main_frame, text=explanation, wraplength=500).pack(pady=(0, 20))

        # Create a notebook (tabbed interface)
        self.notebook = ctk.CTkTabview(main_frame)
        self.notebook.pack(fill="both", expand=True)

        # Emu Database Tab
        emu_tab = self.notebook.add("Emu Database (16S)")
        self.create_emu_tab(emu_tab)

        # Centrifuger Database Tab (renamed)
        centrifuger_tab = self.notebook.add("Centrifuger Database (WGS)")  # Changed
        self.create_centrifuger_tab(centrifuger_tab)  # Changed

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
        create_tooltip(self.emu_path_entry, "Path to the Emu database directory")

        browse_button = ctk.CTkButton(path_frame, text="Browse", command=self.select_emu_db)
        browse_button.pack(side="right", padx=(5, 0))
        create_tooltip(browse_button, "Select the Emu database directory")

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

    def create_centrifuger_tab(self, parent):
        frame = ctk.CTkFrame(parent)
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Database path section
        path_section = ctk.CTkFrame(frame)
        path_section.pack(fill="x", pady=(0, 20))
        
        ctk.CTkLabel(path_section, text="Centrifuger Database Path:", 
                     font=("Helvetica", 12, "bold")).pack(fill="x", pady=(0, 5))
        
        path_frame = ctk.CTkFrame(path_section)
        path_frame.pack(fill="x")
        
        self.centrifuger_path_entry = ctk.CTkEntry(path_frame)
        self.centrifuger_path_entry.pack(side="left", fill="x", expand=True)
        self.centrifuger_path_entry.insert(0, self.centrifuger_db_path)
        
        browse_button = ctk.CTkButton(path_frame, text="Browse", 
                                     command=self.select_centrifuger_db)
        browse_button.pack(side="right", padx=(5, 0))

        # Database Type Selection
        db_type_frame = ctk.CTkFrame(frame)
        db_type_frame.pack(fill="x", pady=(0, 20))
        
        ctk.CTkLabel(db_type_frame, text="Database Type:", 
                     font=("Helvetica", 12, "bold")).pack(fill="x", pady=(0, 5))
        
        self.centrifuger_db_type = ctk.StringVar(value="NCBI")
        
        ctk.CTkRadioButton(db_type_frame, text="NCBI RefSeq", 
                          variable=self.centrifuger_db_type, 
                          value="NCBI",
                          command=self.toggle_domain_selection).pack(anchor="w", pady=2)
        
        ctk.CTkRadioButton(db_type_frame, text="GTDB", 
                          variable=self.centrifuger_db_type, 
                          value="GTDB",
                          command=self.toggle_domain_selection).pack(anchor="w", pady=2)

        # Assembly Level Selection (NCBI only)
        self.assembly_frame = ctk.CTkFrame(frame)
        self.assembly_frame.pack(fill="x", pady=(0, 20))

        # Title and explanation in a visually appealing box
        assembly_header = ctk.CTkFrame(self.assembly_frame)
        assembly_header.pack(fill="x", pady=(0, 10))
        
        ctk.CTkLabel(assembly_header, text="Assembly Level (NCBI only):", 
                     font=("Helvetica", 12, "bold")).pack(anchor="w", pady=(0, 5))
        
        # Create an info box with a border
        info_box = ctk.CTkFrame(assembly_header)
        info_box.pack(fill="x", pady=(0, 5))
        
        # Complete Genome explanation
        complete_frame = ctk.CTkFrame(info_box)
        complete_frame.pack(fill="x", padx=10, pady=5)
        ctk.CTkLabel(complete_frame, text="Complete Genome:", 
                     font=("Helvetica", 11, "bold")).pack(anchor="w")
        ctk.CTkLabel(complete_frame, text="Higher quality reference genomes only",
                     font=("Helvetica", 10)).pack(anchor="w", padx=20)
        
        # Any explanation
        any_frame = ctk.CTkFrame(info_box)
        any_frame.pack(fill="x", padx=10, pady=5)
        ctk.CTkLabel(any_frame, text="Any:", 
                     font=("Helvetica", 11, "bold")).pack(anchor="w")
        ctk.CTkLabel(any_frame, text="Includes both complete and draft genomes",
                     font=("Helvetica", 10)).pack(anchor="w", padx=20)
        
        # Radio buttons for selection
        self.assembly_level = ctk.StringVar(value="Complete Genome")
        
        radio_frame = ctk.CTkFrame(self.assembly_frame)
        radio_frame.pack(fill="x", pady=(5, 0))
        
        ctk.CTkRadioButton(radio_frame, text="Complete Genome", 
                          variable=self.assembly_level, 
                          value="Complete Genome").pack(anchor="w", pady=2)
        
        ctk.CTkRadioButton(radio_frame, text="Any", 
                          variable=self.assembly_level, 
                          value="Any").pack(anchor="w", pady=2)

        # Domain Selection (NCBI only)
        self.domains_frame = ctk.CTkFrame(frame)
        self.domains_frame.pack(fill="x", pady=(20, 0))

        ctk.CTkLabel(self.domains_frame, text="Select Domains (NCBI only):", anchor="w").pack(fill="x", pady=(0, 5))

        self.centrifuger_domains = {
            "bacteria": ctk.BooleanVar(value=False),
            "archaea": ctk.BooleanVar(value=False),
            "viral": ctk.BooleanVar(value=False),
            "fungi": ctk.BooleanVar(value=True),
        }

        for domain, var in self.centrifuger_domains.items():
            ctk.CTkCheckBox(self.domains_frame, text=domain.capitalize(), variable=var).pack(anchor="w", pady=2)

        # Build button
        ctk.CTkButton(frame, text="Build Centrifuger Database", 
                     command=self.build_centrifuger_db).pack(fill="x", pady=(20, 0))

    def toggle_domain_selection(self):
        """Enable/disable domain selection based on database type"""
        if self.centrifuger_db_type.get() == "GTDB":
            for widget in self.domains_frame.winfo_children():
                if isinstance(widget, ctk.CTkCheckBox):
                    widget.configure(state="disabled")
            for widget in self.assembly_frame.winfo_children():
                widget.configure(state="disabled")
        else:
            for widget in self.domains_frame.winfo_children():
                if isinstance(widget, ctk.CTkCheckBox):
                    widget.configure(state="normal")
            for widget in self.assembly_frame.winfo_children():
                widget.configure(state="normal")

    def select_emu_db(self):
        db_path = filedialog.askdirectory(title="Select Emu Database Directory")
        if db_path:
            self.emu_db_path = db_path
            self.emu_path_entry.delete(0, ctk.END)
            self.emu_path_entry.insert(0, db_path)
            self.save_config()

    def select_centrifuger_db(self):  # Changed
        """Select Centrifuger database directory"""  # Changed
        db_path = filedialog.askdirectory(title="Select Centrifuger Database Directory")  # Changed
        if db_path:
            self.centrifuger_path_entry.delete(0, tk.END)  # Changed
            self.centrifuger_path_entry.insert(0, db_path)  # Changed

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

    def confirm_build_centrifuger_db(self):  # Changed
        selected_domains = [domain for domain, var in self.centrifuger_domains.items() if var.get()]  # Changed
        if selected_domains:
            db_type = self.centrifuger_db_type.get()  # Changed
            index_path = self.index_entry.get()
            response = messagebox.askyesno("Confirm Build", 
                f"Building the Centrifuger database ({db_type}) for {', '.join(selected_domains)} requires an internet connection and may take a while. Do you want to proceed?")  # Changed
            if response:
                self.build_centrifuger_db(selected_domains, db_type, index_path)  # Changed

    def build_emu_db(self, selected_domains):
        self.show_progress_window("Building Emu Database")
        threading.Thread(target=self._build_emu_db, args=(selected_domains,)).start()

    def build_centrifuger_db(self, selected_domains, db_type, index_path):  # Changed
        try:
            self.log_progress(f"Creating Centrifuger database directory...")  # Changed
            os.makedirs(self.centrifuger_db_path, exist_ok=True)  # Changed

            with tempfile.TemporaryDirectory() as temp_dir:
                if db_type == "NCBI":
                    # Build Centrifuger database using NCBI data  # Changed
                    index_name = "centrifuger_custom_index"  # Changed
                    taxonomy_dir = os.path.join(self.centrifuger_db_path, "taxonomy")  # Changed
                    library_dir = os.path.join(self.centrifuger_db_path, "library")  # Changed

                    # Set paths for centrifuger executables  # Changed
                    centrifuger_download = os.path.join("centrifuger-download")  # Changed
                    centrifuger_build = os.path.join("centrifuger-build")  # Changed

                    # Download NCBI taxonomy
                    self.update_progress("Downloading NCBI taxonomy...")
                    subprocess.run([centrifuger_download, "-o", taxonomy_dir, "taxonomy"], check=True)  # Changed

                    # Download genomes for selected domains
                    domains = ",".join(selected_domains)
                    self.update_progress(f"Downloading {domains} genomes... This may take a while.")
                    subprocess.run([centrifuger_download, "-o", library_dir, "-m", "-d", domains, "refseq"],  # Changed
                                check=True, stdout=subprocess.PIPE)

                    # Build database
                    num_threads = multiprocessing.cpu_count()
                    self.update_progress("Building Centrifuger database...")  # Changed
                    build_cmd = [
                        centrifuger_build,  # Changed
                        "-t", str(num_threads),
                        "--conversion-table", "seqid2taxid.map",
                        "--taxonomy-tree", "taxonomy/nodes.dmp",
                        "--name-table", "taxonomy/names.dmp",
                        "-f", "file.list",
                        "-o", "refseq_abv",
                        "--build-mem", str(build_mem)
                    ]
                    subprocess.run(build_cmd, check=True, capture_output=True, text=True)

                elif db_type == "GTDB":
                    # Build Centrifuger database using GTDB data
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
            messagebox.showinfo("Success", f"Centrifuger database built successfully at {self.centrifuger_db_path}")

            # Ask user if they want to open the folder
            self.ask_to_open_folder(self.centrifuger_db_path)

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

    def show_progress_window(self, title):
        """Create and show progress window safely"""
        try:
            self.progress_window = ctk.CTkToplevel(self)
            self.progress_window.title(title)
            self.progress_window.geometry("400x150")
            self.progress_window.transient(self)
            self.progress_window.grab_set()
            
            # Create label first
            self.progress_label = ctk.CTkLabel(self.progress_window, text="Initializing...")
            self.progress_label.pack(pady=20)
            
            # Create progress bar
            self.progress_bar = ctk.CTkProgressBar(self.progress_window, mode="indeterminate")
            self.progress_bar.pack(pady=20, padx=20, fill="x")
            self.progress_bar.start()
            
            # Update immediately to ensure widgets are created
            self.progress_window.update()
        except Exception as e:
            print(f"Error creating progress window: {e}")

    def update_progress(self, message):
        """Update progress window message safely and log it"""
        def safe_update():
            try:
                if hasattr(self, 'progress_label') and self.progress_label and self.progress_label.winfo_exists():
                    self.progress_label.configure(text=message)
                self.log_progress(message)
            except Exception as e:
                print(f"Error updating progress: {e}")
        
        if self.winfo_exists():
            self.after(0, safe_update)

    def log_progress(self, message):
        """Log message to file and console"""
        self.logger.info(message)

    def close_progress_window(self):
        """Safely close the progress window"""
        try:
            if hasattr(self, 'progress_window') and self.progress_window is not None:
                self.progress_window.destroy()
                self.progress_window = None
        except Exception as e:
            print(f"Error closing progress window: {e}")

    def show_error_message(self, message):
        def show_message():
            messagebox.showerror("Error", message)
        self.after(0, show_message)

    def load_config(self):
        if os.path.exists(self.config_file):
            with open(self.config_file, 'r') as f:
                config = json.load(f)
                self.emu_db_path = config.get('emu_db', self.emu_db_path)
                self.centrifuger_db_path = config.get('centrifuger_db', self.centrifuger_db_path)

    def save_config(self):
        """Save current configuration to pipeline_config.json"""
        try:
            # First read existing config to preserve other settings
            if os.path.exists(self.config_file):
                with open(self.config_file, 'r') as f:
                    config = json.load(f)
            else:
                config = {}

            # Update database paths
            config.update({
                'emu_db': self.emu_db_path,
                'centrifuge_db': self.centrifuger_db_path,
                'gtdb_db': self.gtdb_db_path
            })

            # Save updated config
            with open(self.config_file, 'w') as f:
                json.dump(config, f, indent=4)

            # Notify PipelineConfig to refresh if it exists
            if hasattr(self.parent, 'pipeline_config'):
                self.parent.pipeline_config.load_config()
                
        except Exception as e:
            messagebox.showerror("Error", f"Could not save configuration: {str(e)}")

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

    def build_centrifuger_db(self):
        """Build Centrifuger database using selected domains"""
        try:
            selected_domains = [domain for domain, var in self.centrifuger_domains.items() 
                              if var.get()]
            
            if not selected_domains:
                messagebox.showerror("Error", "Please select at least one domain.")
                return

            # Create progress window
            self.show_progress_window("Building Centrifuger Database")
            
            if self.centrifuger_db_type.get() == "NCBI":
                threading.Thread(target=self._build_ncbi_db, args=(selected_domains,)).start()
            else:  # GTDB
                threading.Thread(target=self._build_gtdb_db).start()
                
        except Exception as e:
            self.close_progress_window()
            error_msg = f"An error occurred:\n{str(e)}"
            messagebox.showerror("Error", error_msg)

    def _build_ncbi_db(self, selected_domains):
        """Build Centrifuger database using NCBI data"""
        build_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        build_dir = os.path.join(self.centrifuger_db_path, f"ncbi_build_{build_time}")
        
        try:
            # Create build directory
            os.makedirs(build_dir, exist_ok=True)
            self.update_progress(f"Created build directory: {build_dir}")
            
            # Change to build directory
            original_dir = os.getcwd()
            os.chdir(build_dir)
            
            try:
                # Download taxonomy files
                self.update_progress("Downloading NCBI taxonomy files...")
                taxonomy_result = subprocess.run(
                    ["centrifuger-download", "-o", "taxonomy", "taxonomy"], 
                    check=True, capture_output=True, text=True
                )
                if taxonomy_result.stdout:
                    self.log_progress(taxonomy_result.stdout)
                
                # Create library directory
                os.makedirs("library", exist_ok=True)
                
                # Download sequences for each domain
                domain_string = ",".join(selected_domains)
                self.update_progress(f"Downloading {domain_string} sequences...")
                
                with open("seqid2taxid.map", "w") as seq2taxid_file:
                    # Run download command and capture output
                    download_cmd = [
                        "centrifuger-download",
                        "-d", domain_string,
                        "-o", "library",
                        "-P", str(multiprocessing.cpu_count()),
                        "-a", self.assembly_level.get(),
                        "refseq",  
                    ]
                    
                    process = subprocess.Popen(
                        download_cmd,
                        stdout=seq2taxid_file,
                        stderr=subprocess.PIPE,
                        text=True,
                        universal_newlines=True,
                        bufsize=1
                    )
                    
                    # Read stderr for progress updates
                    while True:
                        error = process.stderr.readline()
                        if error:
                            if not error.startswith("Progress"):
                                self.log_progress(error.strip())
                        if process.poll() is not None:
                            break
                    
                    # Wait for process to complete
                    process.wait()
                    
                    # Check return code
                    if process.returncode != 0:
                        raise subprocess.CalledProcessError(
                            process.returncode,
                            download_cmd
                        )
                
                # Create file list using glob after download completes
                self.update_progress("Creating sequence file list...")
                import glob
                sequence_files = glob.glob(os.path.join("library", "*", "*.fna.gz"))
                
                if not sequence_files:
                    raise Exception("No sequence files found in library directory")
                
                with open("file.list", "w") as f:
                    for file_path in sequence_files:
                        f.write(f"{file_path}\n")
                
                file_count = len(sequence_files)
                self.log_progress(f"Found {file_count} sequence files")
                
                # Continue with database build...
                self.update_progress("Building Centrifuger database (this may take a while)...")
                build_mem = int(psutil.virtual_memory().total * 0.9 / (1024 * 1024 * 1024))
                
                build_cmd = [
                    "centrifuger-build",
                    "-l", "file.list",
                    "--taxonomy-tree", "taxonomy/nodes.dmp",
                    "--name-table", "taxonomy/names.dmp",
                    "-t", str(multiprocessing.cpu_count()),
                    "--build-mem", f"{build_mem}G",
                    "--conversion-table", "seqid2taxid.map",
                    "-o", "cfr_ncbi"
                ]
                
                build_result = subprocess.run(build_cmd, check=True, capture_output=True, text=True)
                if build_result.stdout:
                    self.log_progress(build_result.stdout)
                if build_result.stderr and not build_result.stderr.startswith("Progress"):
                    self.log_progress(build_result.stderr)
                
                # Move necessary files to final location
                self.update_progress("Moving index files to final location...")
                index_files = ["cfr_ncbi.1.cf", "cfr_ncbi.2.cf", "cfr_ncbi.3.cf", "cfr_ncbi.4.cf"]
                for file in index_files:
                    if os.path.exists(file):
                        shutil.move(file, os.path.join(build_dir, file))
                
                success_message = (
                    "Database built successfully!\n\n"
                    f"Location: {build_dir}\n\n"
                    f"Memory used: {build_mem}GB\n"
                    f"Sequences processed: {file_count}\n\n"
                    "The database path has been saved in your configuration."
                )
                
                def show_success():
                    if self.winfo_exists():
                        # Close progress window first
                        self.close_progress_window()
                        # Then show success message
                        response = messagebox.showinfo(
                            "Success", 
                            success_message,
                            detail="Would you like to open the folder?"
                        )
                        if response:
                            self.open_folder(build_dir)
                self.after(0, show_success)
                
                # Update configuration with the new database path
                self.centrifuger_db_path = build_dir
                self.save_config()
                
            finally:
                os.chdir(original_dir)
                
        except Exception as error:
            error_msg = f"Error building database: {str(error)}"
            self.log_progress(error_msg)
            
            def show_error():
                if self.winfo_exists():
                    # Close progress window first
                    self.close_progress_window()
                    # Then show error message
                    messagebox.showerror(
                        "Error", 
                        "Database build failed.\n\n"
                        f"Error details:\n{error_msg}"
                    )
            self.after(0, show_error)

    def _build_gtdb_db(self):
        """Build Centrifuger database using GTDB data"""
        try:
            build_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            build_dir = os.path.join(self.centrifuger_db_path, f"gtdb_build_{build_time}")
            os.makedirs(build_dir, exist_ok=True)
            
            # Change to build directory
            original_dir = os.getcwd()
            os.chdir(build_dir)
            
            try:
                self.update_progress("Starting GTDB database build...")
                
                # Base URL for GTDB data
                base_url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"
                
                # Download and extract genome representatives
                self.update_progress("Downloading GTDB genome representatives...")
                genome_file = "gtdb_genomes_reps.tar.gz"
                self._download_file(f"{base_url}genomic_files_reps/{genome_file}", genome_file)
                
                self.update_progress("Extracting genome files...")
                with tarfile.open(genome_file, "r:gz") as tar:
                    tar.extractall()
                
                # Get GTDB version
                self.update_progress("Getting GTDB version...")
                version_text = requests.get(f"{base_url}VERSION.txt").text.strip()
                gtdb_version = version_text[1:]  # Remove 'R' prefix
                
                # Download and combine metadata
                self.update_progress("Downloading and processing metadata...")
                meta_file = "gtdb_meta.tsv"
                
                # Download and decompress bacterial metadata
                bac_response = requests.get(f"{base_url}bac120_metadata.tsv.gz")
                bac_data = gzip.decompress(bac_response.content).decode()
                
                # Download and decompress archaeal metadata
                arch_response = requests.get(f"{base_url}ar53_metadata.tsv.gz")
                arch_data = gzip.decompress(arch_response.content).decode()
                
                # Combine metadata
                with open(meta_file, 'w') as f:
                    f.write(bac_data)
                    # Add archaeal data without header
                    arch_lines = arch_data.split('\n')
                    f.write('\n'.join(arch_lines[1:]))
                
                # Create taxonomy files
                self.update_progress("Creating taxonomy files...")
                self._create_gtdb_taxonomy(meta_file, gtdb_version)
                
                # Build Centrifuger index
                self.update_progress("Building Centrifuger index...")
                num_threads = multiprocessing.cpu_count()
                build_mem = int(psutil.virtual_memory().total * 0.9 / (1024 * 1024 * 1024))
                
                build_cmd = [
                    "centrifuger-build",
                    "-l", "gtdb_fname_to_taxid.map",
                    "--taxonomy-tree", "gtdb_nodes.dmp",
                    "--name-table", "gtdb_names.dmp",
                    "-t", str(num_threads),
                    "-o", "cfr_gtdb",
                    "--build-mem", f"{build_mem}G",
                    "--offrate", "5",
                    "--rbbwt-b", "1"
                ]
                
                subprocess.run(build_cmd, check=True)
                
                # Cleanup temporary files
                self.update_progress("Cleaning up temporary files...")
                for temp_file in [genome_file, meta_file]:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
                
                success_message = (
                    "GTDB database built successfully!\n\n"
                    f"Location: {build_dir}\n"
                    f"GTDB Version: {gtdb_version}\n"
                    f"Memory used: {build_mem}GB\n"
                    "The database path has been saved in your configuration."
                )
                
                def show_success():
                    if self.winfo_exists():
                        response = messagebox.showinfo(
                            "Success", 
                            success_message,
                            detail="Would you like to open the folder?"
                        )
                        if response:
                            self.open_folder(build_dir)
                self.after(0, show_success)
                
                # Update configuration
                self.gtdb_db_path = build_dir
                self.save_config()
                
            finally:
                os.chdir(original_dir)
                
        except Exception as error:
            error_msg = f"Error building GTDB database: {str(error)}\n{traceback.format_exc()}"
            self.log_progress(error_msg)
            
            # Cleanup on failure
            self.update_progress("Cleaning up after error...")
            try:
                if os.path.exists(build_dir):
                    shutil.rmtree(build_dir)
                self.log_progress("Cleanup completed")
            except Exception as cleanup_error:
                self.log_progress(f"Error during cleanup: {cleanup_error}")
            
            def show_error():
                if self.winfo_exists():
                    messagebox.showerror(
                        "Error", 
                        "GTDB database build failed. All temporary files have been cleaned up.\n\n"
                        f"Error details:\n{error_msg}"
                    )
            self.after(0, show_error)

    def _download_file(self, url, filename):
        """Download file with progress updates"""
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        block_size = 8192
        
        with open(filename, 'wb') as f:
            downloaded = 0
            for data in response.iter_content(block_size):
                downloaded += len(data)
                f.write(data)
                progress = (downloaded / total_size) * 100 if total_size > 0 else 0
                self.update_progress(f"Downloading {filename}: {progress:.1f}%")

    def _create_gtdb_taxonomy(self, meta_file, version):
        """Create GTDB taxonomy files from metadata"""
        taxonomy_data = {}
        name_to_taxid = {}
        next_taxid = 1
        
        with open(meta_file) as f:
            header = f.readline().strip().split('\t')
            gtdb_taxonomy_idx = header.index('gtdb_taxonomy')
            
            for line in f:
                fields = line.strip().split('\t')
                taxonomy = fields[gtdb_taxonomy_idx].split(';')
                
                current_taxid = None
                parent_taxid = 0
                
                for rank, name in enumerate(taxonomy):
                    if name not in name_to_taxid:
                        name_to_taxid[name] = next_taxid
                        next_taxid += 1
                    
                    current_taxid = name_to_taxid[name]
                    if current_taxid not in taxonomy_data:
                        taxonomy_data[current_taxid] = {
                            'name': name,
                            'parent': parent_taxid,
                            'rank': self._get_rank_name(rank)
                        }
                    parent_taxid = current_taxid
        
        # Write nodes.dmp
        with open('gtdb_nodes.dmp', 'w') as f:
            for taxid, data in taxonomy_data.items():
                f.write(f"{taxid}\t|\t{data['parent']}\t|\t{data['rank']}\t|\n")
        
        # Write names.dmp
        with open('gtdb_names.dmp', 'w') as f:
            for taxid, data in taxonomy_data.items():
                f.write(f"{taxid}\t|\t{data['name']}\t|\t\t|\tscientific name\t|\n")

    def _get_rank_name(self, rank):
        """Get taxonomic rank name based on position"""
        ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        return ranks[rank] if rank < len(ranks) else 'no rank'

    def open_folder(self, path):
        """Open the folder in the system's file explorer"""
        try:
            if sys.platform == "win32":
                os.startfile(path)
            elif sys.platform == "darwin":  # macOS
                subprocess.run(["open", path])
            else:  # linux
                subprocess.run(["xdg-open", path])
        except Exception as e:
            messagebox.showerror("Error", f"Could not open folder: {str(e)}")

    def save_config(self):
        """Save current configuration to file"""
        config = {
            'emu_db': self.emu_db_path,
            'centrifuger_db': self.centrifuger_db_path,
            'gtdb_db': self.gtdb_db_path
        }
        try:
            with open(self.config_file, 'w') as f:
                json.dump(config, f, indent=4)
        except Exception as e:
            messagebox.showerror("Error", f"Could not save configuration: {str(e)}")


















