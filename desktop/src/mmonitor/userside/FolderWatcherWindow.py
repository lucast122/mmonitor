import os
import datetime
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import customtkinter as ctk
from .FolderMonitor import FolderMonitor
import subprocess
import sys
import gzip
import shutil
import queue
import threading
from build_mmonitor_pyinstaller import ROOT
from tkcalendar import Calendar
import re
import json
from mmonitor.userside.utils import create_tooltip
from mmonitor.userside.MMonitorCMD import MMonitorCMD
from mmonitor.userside.PipelineConfig import PipelineConfig
from mmonitor.database.django_db_interface import DjangoDBInterface
from Bio import SeqIO
import numpy as np
import time
import multiprocessing
import csv
import tempfile


class FolderWatcherWindow(ctk.CTkFrame):
    def __init__(self, parent, gui_ref):
        print("Initializing FolderWatcherWindow")
        super().__init__(parent)
        self.parent = parent
        self.gui = gui_ref
        
        # Get authentication state from main window
        self.offline_mode = self.gui.offline_mode
        self.logged_in = self.gui.logged_in
        self.current_user = self.gui.current_user
        
        # Initialize variables
        self.samples = {}
        self.watching = False
        self.folder_monitor = None
        self.auto_analyze_var = tk.BooleanVar(value=False)
        
        self.config = {}
        self.watch_start_time = None
        self.all_samples = {}  # All samples including existing ones
        self.queue_samples = {}  # Only new samples after watching started
        self.pending_files = []  # Initialize pending_files list
        self.config_file = os.path.join(ROOT, "src", "resources", "pipeline_config.json")
        
        # Load configuration first
        self.load_config()
        if not self.config:
            messagebox.showwarning(
                "Configuration Missing",
                "Please configure analysis parameters in the Configuration tab first."
            )
            return
        
        # UI update control
        self.last_update_time = time.time()
        self.needs_update = False
        self.expanded_items = set()
        
        # Create UI
        self.create_widgets()
        self.after(100, self.process_queue)
        
        # Bind to parent window close event
        if isinstance(parent, tk.Tk) or isinstance(parent, tk.Toplevel):
            parent.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        print("FolderWatcherWindow initialized")

    def load_config(self):
        """Load analysis configuration from file"""
        if os.path.exists(self.config_file):
            try:
                with open(self.config_file, 'r') as f:
                    self.config = json.load(f)
                print(f"Loaded configuration: {self.config}")
                return True
            except Exception as e:
                print(f"Error loading configuration: {e}")
                return False
        else:
            print("No configuration file found")
            return False

    def create_widgets(self):
        # Configure self first
        is_dark = ctk.get_appearance_mode() == "dark"
        bg_color = "transparent"  # Use transparent background
        
        self.configure(fg_color=bg_color)
        
        # Main container with padding and rounded corners
        self.main_container = ctk.CTkFrame(self, fg_color=bg_color)
        self.main_container.pack(fill="both", expand=True, pady=10)

        # Add title section with modern typography
        title_frame = ctk.CTkFrame(self.main_container, fg_color=bg_color)
        title_frame.pack(fill="x", pady=(0, 10))
        ctk.CTkLabel(title_frame, text="Folder Watcher", 
                     font=("Helvetica", 24, "bold")).pack(pady=5)
        ctk.CTkLabel(title_frame, text="Monitor folders for new sequencing data",
                     font=("Helvetica", 12)).pack()

        # Create sections with proper spacing
        self.create_top_section(self.main_container)
        self.create_separator(self.main_container)
        self.create_file_trees(self.main_container)
        self.create_separator(self.main_container)
        self.create_bottom_controls(self.main_container)

        # Configure Treeview style
        style = ttk.Style()
        style.configure(
            "Custom.Treeview",
            background="#ffffff",
            fieldbackground="#ffffff",
            foreground="#000000",
            font=("Helvetica", 11)
        )
        style.configure(
            "Custom.Treeview.Heading",
            background="#f0f0f0",
            foreground="#000000",
            font=("Helvetica", 11, "bold")
        )
        # Configure selection colors
        style.map(
            "Custom.Treeview",
            background=[("selected", "#e1e1e1")],
            foreground=[("selected", "#000000")]
        )

    def create_separator(self, parent):
        """Create a visual separator line"""
        separator = ctk.CTkFrame(parent, height=2)
        separator.pack(fill="x", pady=10)

    def create_top_section(self, container):
        top_frame = ctk.CTkFrame(container, fg_color="transparent")
        top_frame.pack(fill="x", pady=(0, 5))

        # Folder selection with improved visual grouping
        folder_group = ctk.CTkFrame(top_frame, fg_color="transparent")
        folder_group.pack(fill="x", pady=3)
        
        folder_label = ctk.CTkLabel(folder_group, text="Folder to watch:", 
                                   font=("Helvetica", 12, "bold"))
        folder_label.pack(side="left")
        
        self.folder_entry = ctk.CTkEntry(folder_group, height=32)
        self.folder_entry.pack(side="left", expand=True, fill="x", padx=5)
        
        button_frame = ctk.CTkFrame(folder_group, fg_color="transparent")
        button_frame.pack(side="left")
        
        # Modern buttons with consistent styling
        self.browse_button = ctk.CTkButton(
            button_frame, 
            text="Browse", 
            command=self.browse_folder,
            height=32,
            font=("Helvetica", 12)
        )
        self.browse_button.pack(side="left", padx=2)
        
        self.csv_button = ctk.CTkButton(
            button_frame, 
            text="Import CSV", 
            command=self.load_from_csv,
            height=32,
            font=("Helvetica", 12)
        )
        self.csv_button.pack(side="left", padx=2)
        
        self.watch_button = ctk.CTkButton(
            button_frame, 
            text="Start Watching", 
            command=self.toggle_watching,
            height=32,
            font=("Helvetica", 12)
        )
        self.watch_button.pack(side="left", padx=2)

        # Project information section
        project_frame = ctk.CTkFrame(top_frame, fg_color="transparent")
        project_frame.pack(fill="x", pady=(10, 0))

        # Project name
        project_group = ctk.CTkFrame(project_frame, fg_color="transparent")
        project_group.pack(fill="x", pady=2)
        
        ctk.CTkLabel(project_group, text="Project:", 
                     font=("Helvetica", 12, "bold"), 
                     width=80).pack(side="left")
        
        self.project_entry = ctk.CTkEntry(project_group, height=32)
        self.project_entry.pack(side="left", expand=True, fill="x", padx=5)

        # Subproject name
        subproject_group = ctk.CTkFrame(project_frame, fg_color="transparent")
        subproject_group.pack(fill="x", pady=2)
        
        ctk.CTkLabel(subproject_group, text="Subproject:", 
                     font=("Helvetica", 12, "bold"), 
                     width=80).pack(side="left")
        
        self.subproject_entry = ctk.CTkEntry(subproject_group, height=32)
        self.subproject_entry.pack(side="left", expand=True, fill="x", padx=5)

        # Date selection section
        date_frame = ctk.CTkFrame(top_frame, fg_color="transparent")
        date_frame.pack(fill="x", pady=(10, 0))

        # Date picker
        date_group = ctk.CTkFrame(date_frame, fg_color="transparent")
        date_group.pack(fill="x", pady=2)
        
        ctk.CTkLabel(date_group, text="Date:", 
                     font=("Helvetica", 12, "bold"), 
                     width=80).pack(side="left")
        
        # Initialize date variables
        self.selected_date = datetime.date.today()
        self.use_file_date = tk.BooleanVar(value=True)
        
        # Date entry and calendar button
        date_entry_frame = ctk.CTkFrame(date_frame, fg_color="transparent")
        date_entry_frame.pack(side="left", expand=True, fill="x", padx=5)
        
        self.date_entry = ctk.CTkEntry(date_entry_frame, height=32)
        self.date_entry.pack(side="left", expand=True, fill="x")
        self.date_entry.insert(0, self.selected_date.strftime('%Y-%m-%d'))
        
        calendar_button = ctk.CTkButton(
            date_entry_frame,
            text="📅",
            width=32,
            height=32,
            command=self.show_calendar
        )
        calendar_button.pack(side="left", padx=(5, 0))

        # Use file date checkbox
        checkbox_frame = ctk.CTkFrame(date_frame, fg_color="transparent")
        checkbox_frame.pack(fill="x", pady=2)
        
        self.use_file_date_checkbox = ctk.CTkCheckBox(
            checkbox_frame,
            text="Use file creation date",
            variable=self.use_file_date,
            font=("Helvetica", 11)
        )
        self.use_file_date_checkbox.pack(side="left", padx=(80, 0))

    def create_file_trees(self, container):
        # Create a frame for both trees with proper spacing
        trees_frame = ctk.CTkFrame(container, fg_color="transparent")
        trees_frame.pack(fill="both", expand=True, pady=5)
        
        # Left tree (70% width)
        left_frame = ctk.CTkFrame(trees_frame, fg_color="transparent")
        left_frame.pack(side="left", fill="both", expand=True, padx=(0, 2))
        
        # Title and count for left tree
        left_header = ctk.CTkFrame(left_frame, fg_color="transparent")
        left_header.pack(fill="x", pady=2)
        ctk.CTkLabel(left_header, text="All Files", 
                     font=("Helvetica", 14, "bold")).pack(side="left")
        
        # Create the file tree with updated style
        self.file_tree = ttk.Treeview(
            left_frame,
            columns=("files", "status"),
            selectmode="extended",
            style="Custom.Treeview",
            height=15
        )
        self.file_tree.pack(fill="both", expand=True)
        
        # Configure file tree columns
        self.file_tree.heading("#0", text="Sample/File")
        self.file_tree.heading("files", text="Files")
        self.file_tree.heading("status", text="Status")
        
        # Right tree (30% width)
        right_frame = ctk.CTkFrame(trees_frame, fg_color="transparent")
        right_frame.pack(side="right", fill="both", padx=(2, 0))
        
        # Title and count for right tree
        right_header = ctk.CTkFrame(right_frame, fg_color="transparent")
        right_header.pack(fill="x", pady=2)
        ctk.CTkLabel(right_header, text="Queue", 
                     font=("Helvetica", 14, "bold")).pack(side="left")
        self.queue_count = ctk.CTkLabel(right_header, text="(0 samples, 0 files)", 
                                       font=("Helvetica", 12))
        self.queue_count.pack(side="left", padx=5)
        
        # Create the queue tree with updated style
        self.queue_tree = ttk.Treeview(
            right_frame,
            columns=("files", "first_seen"),
            selectmode="extended",
            style="Custom.Treeview",
            height=15
        )
        self.queue_tree.pack(fill="both", expand=True)
        
        # Configure queue tree columns
        self.queue_tree.heading("#0", text="Sample")
        self.queue_tree.heading("files", text="Files")
        self.queue_tree.heading("first_seen", text="First Seen")
        
        # Configure column widths
        self.file_tree.column("#0", width=200)
        self.file_tree.column("files", width=100)
        self.file_tree.column("status", width=100)
        
        self.queue_tree.column("#0", width=150)
        self.queue_tree.column("files", width=100)
        self.queue_tree.column("first_seen", width=100)

    def create_bottom_controls(self, container):
        bottom_frame = ctk.CTkFrame(container)
        bottom_frame.pack(fill="x", pady=5)

        # Left controls
        self.auto_analyze_checkbox = ctk.CTkCheckBox(bottom_frame, 
                                                   text="Auto-analyze new samples",
                                                   variable=self.auto_analyze_var,
                                                   font=("Helvetica", 11))
        self.auto_analyze_checkbox.pack(side="left", padx=3)

        # Right controls with smaller buttons
        right_controls = ctk.CTkFrame(bottom_frame)
        right_controls.pack(side="right")
        
        button_params = {
            "height": 28,
            "font": ("Helvetica", 11)
        }
        
        ctk.CTkButton(right_controls, text="Analyze Selected", 
                     command=self.run_selected_sample, 
                     **button_params).pack(side="left", padx=2)
        
        ctk.CTkButton(right_controls, text="Analyze Queue", 
                     command=self.analyze_selected_queue, 
                     **button_params).pack(side="left", padx=2)
        
        ctk.CTkButton(right_controls, text="Clear Queue", 
                     command=self.clear_queue, 
                     **button_params).pack(side="left", padx=2)

    def add_new_file(self, file_path):
        """Add a new file to both trees if appropriate"""
        # Skip concatenated files
        if '_concatenated' in file_path:
            return
        
        if not self.watch_start_time:
            return
        
        # Store currently expanded items before update
        self.store_expanded_state()
        
        print(f"Processing new file: {file_path}")
        sample_name = self.suggest_sample_name(os.path.dirname(file_path))
        current_time = time.strftime("%H:%M:%S", time.localtime(time.time()))
        
        # Add to all_samples (left tree)
        if sample_name not in self.all_samples:
            print(f"Adding new sample {sample_name}")
            self.all_samples[sample_name] = {
                "files": [file_path],
                "path": os.path.dirname(file_path)
            }
            self.samples[sample_name] = [file_path]
        else:
            if file_path not in self.all_samples[sample_name]["files"]:
                print(f"Adding new file to existing sample {sample_name}")
                self.all_samples[sample_name]["files"].append(file_path)
                if sample_name in self.samples:
                    self.samples[sample_name].append(file_path)
                else:
                    self.samples[sample_name] = [file_path]
        
        # Add to queue_samples (right tree)
        if sample_name not in self.queue_samples:
            print(f"Adding sample {sample_name} to queue")
            self.queue_samples[sample_name] = {
                "files": [file_path],
                "first_seen": current_time
            }
        else:
            if file_path not in self.queue_samples[sample_name]["files"]:
                print(f"Adding new file to queued sample {sample_name}")
                self.queue_samples[sample_name]["files"].append(file_path)
        
        # Set needs_update flag instead of forcing immediate updates
        self.needs_update = True
        
        print(f"Current samples in all_samples: {list(self.all_samples.keys())}")
        print(f"Current samples in queue: {list(self.queue_samples.keys())}")

    def schedule_updates(self):
        """Schedule updates to avoid redundant calls"""
        self.after_cancel(self.process_queue) if hasattr(self, '_after_id') else None
        self._after_id = self.after(100, self.process_queue)

    def update_file_tree(self):
        """Update the main file tree view"""
        try:
            # Store expanded state
            self.store_expanded_state()
            
            # Clear existing items
            self.file_tree.delete(*self.file_tree.get_children())
            
            # Add samples and their files
            for sample_name, data in self.all_samples.items():
                # Add sample
                sample_id = self.file_tree.insert("", "end", text=sample_name,
                                    values=(f"{len(data['files'])} files",
                                           data['path']))
                
                # Add files under sample with unique IDs
                for i, file_path in enumerate(data['files']):
                    try:
                        # Use index to ensure unique IDs
                        file_id = f"{sample_name}_file_{i}"
                        self.file_tree.insert(sample_id, "end", 
                                            iid=file_id,
                                            text=os.path.basename(file_path),
                                            values=("", ""))
                    except Exception as e:
                        print(f"Error displaying file {file_path}: {e}")
            
            # Restore expanded state
            self.restore_expanded_state()
                
        except Exception as e:
            print(f"Error updating file tree: {e}")

    def run_analysis_for_sample(self, sample_name):
        """Run analysis for a single sample using MMonitorCMD"""
        if not self.config:
            messagebox.showerror("Error", "No configuration loaded. Please configure analysis parameters first.")
            return

        if sample_name not in self.all_samples:
            print(f"Error: Sample {sample_name} not found in all_samples")
            return

        try:
            # Get analysis parameters from config
            if 'analysis_type' not in self.config:
                messagebox.showerror("Error", "Analysis type not found in configuration.")
                return
                
            analysis_type = self.config['analysis_type']
            print(f"Starting {analysis_type} analysis for {sample_name}")
            
            # Get sample data
            sample_data = self.all_samples[sample_name]
            sample_folder = sample_data["path"]  # This is the actual folder containing the files
            sample_files = sample_data["files"]  # These are the full paths to the files
            
            print(f"Processing sample {sample_name} from folder: {sample_folder}")
            print(f"Found {len(sample_files)} files to analyze")
            
            try:
                # Update tree item status
                self.file_tree.item(sample_name, values=(f"{len(sample_files)} files", "Running..."))
                self.file_tree.update()
                
                # Create MMonitorCMD instance with authentication state
                cmd_runner = MMonitorCMD()
                cmd_runner.offline_mode = self.offline_mode
                cmd_runner.logged_in = self.logged_in
                cmd_runner.current_user = self.current_user
                
                # Create temporary concatenated file
                with tempfile.NamedTemporaryFile(delete=False, suffix='.fastq.gz') as temp_file:
                    print(f"Creating concatenated file: {temp_file.name}")
                    # Concatenate input files
                    cmd_runner.centrifuger_runner.concatenate_gzipped_files(sample_files, temp_file.name)
                    
                    # Run Centrifuger analysis
                    success = cmd_runner.centrifuger_runner.run_centrifuger(
                        input_file=temp_file.name,
                        sample_name=sample_name,
                        db_path=self.config['centrifuger_db']  # Use the configured database path
                    )
                    
                    # Clean up temporary file
                    os.unlink(temp_file.name)
                    
                    # Update status based on success
                    if success:
                        self.file_tree.item(sample_name, values=(f"{len(sample_files)} files", "Completed"))
                        messagebox.showinfo("Success", f"Analysis completed for sample {sample_name}")
                    else:
                        self.file_tree.item(sample_name, values=(f"{len(sample_files)} files", "Failed"))
                        messagebox.showerror("Error", f"Analysis failed for sample {sample_name}")
                    
                    # Remove from queue if present
                    if sample_name in self.queue_samples:
                        del self.queue_samples[sample_name]
                        self.update_queue_display()
                        
            except tk.TclError as e:
                print(f"GUI Error updating tree: {e}")
                # Continue with analysis even if GUI update fails
                
        except Exception as e:
            error_msg = f"Error running analysis for {sample_name}: {str(e)}"
            print(error_msg)
            import traceback
            traceback.print_exc()
            messagebox.showerror("Error", error_msg)

    def check_and_start_analysis(self, sample_name):
        """Check if conditions are met to start analysis"""
        if not self.auto_analyze_var.get():
            return
        
        if sample_name not in self.all_samples:
            return
        
        # Get minimum files threshold from config or use default
        min_files = self.config.get('min_files_for_analysis', 10)
        current_files = len(self.all_samples[sample_name]["files"])
        
        print(f"Checking analysis conditions for {sample_name}: {current_files}/{min_files} files")
        
        if current_files >= min_files:
            print(f"Starting analysis for {sample_name}")
            self.run_analysis_for_sample(sample_name)

    def browse_folder(self):
        folder = filedialog.askdirectory()
        if folder:
            self.folder_entry.delete(0, tk.END)
            self.folder_entry.insert(0, folder)
            self.suggest_project_names(folder)
            # Auto-start watching when folder is selected
            if not self.watching:
                self.start_watching()

    def suggest_project_names(self, folder):
        """Suggest project and subproject names based on the folder path"""
        folder_parts = folder.split(os.sep)
        if len(folder_parts) >= 2:
            project_name = folder_parts[-1]
            subproject_name = folder_parts[-2]
        else:
            project_name = folder_parts[-1]
            subproject_name = ""

        self.project_entry.delete(0, tk.END)
        self.project_entry.insert(0, project_name)
        
        self.subproject_entry.delete(0, tk.END)
        self.subproject_entry.insert(0, subproject_name)

    def toggle_watching(self):
        if not self.watching:
            self.start_watching()
        else:
            self.stop_watching()

    def start_watching(self):
        if not self.check_config():
            return

        folder = self.folder_entry.get()
        if not folder:
            messagebox.showerror("Error", "Please select a folder to watch.")
            return

        # Clear existing data
        self.file_tree.delete(*self.file_tree.get_children())
        self.queue_tree.delete(*self.queue_tree.get_children())
        self.all_samples.clear()
        self.samples.clear()
        self.queue_samples.clear()
        
        # Set watch start time
        self.watch_start_time = time.time()
        
        # Now scan existing files
        self.scan_existing_files(folder)
        self.update_file_tree()  # Force initial update

        # Start monitoring
        if self.folder_monitor:
            self.folder_monitor.stop()

        self.folder_monitor = FolderMonitor([folder], self)
        self.folder_monitor.start()
        self.watching = True
        self.watch_button.configure(text="Stop Watching")
        
        # Disable folder entry and browse button while watching
        self.folder_entry.configure(state="disabled")
        self.browse_button.configure(state="disabled")

        messagebox.showinfo("Folder Watch Started", f"Now watching folder: {folder}")

    def stop_watching(self):
        if self.folder_monitor:
            self.folder_monitor.stop()
            self.folder_monitor = None
        self.watching = False
        self.watch_button.configure(text="Start Watching")
        
        # Re-enable folder entry and browse button
        self.folder_entry.configure(state="normal")
        self.browse_button.configure(state="normal")
        
        messagebox.showinfo("Folder Watch Stopped", "Stopped watching folder.")

    def scan_existing_files(self, folder):
        """Scan and add existing files to the file tree only"""
        print(f"Scanning existing files in {folder}")
        self.file_tree.delete(*self.file_tree.get_children())
        self.all_samples.clear()
        self.samples.clear()
        
        for root, _, files in os.walk(folder):
            # Filter out concatenated files and get only fastq files
            fastq_files = [f for f in files 
                          if (f.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')) and 
                              '_concatenated' not in f)]
            
            if fastq_files:
                sample_name = self.suggest_sample_name(root)
                file_paths = [os.path.join(root, f) for f in fastq_files]
                print(f"Found sample {sample_name} with {len(file_paths)} files")
                
                self.all_samples[sample_name] = {
                    "files": file_paths,
                    "path": root
                }
                self.samples[sample_name] = file_paths
        
        print(f"Total samples found: {len(self.all_samples)}")
        self.needs_update = True  # Set flag instead of immediate update

    """
    Detect if the reads look like 16S or WGS based on read length distribution.
    """
    
    def detect_read_type(self, file_path, length_threshold=2000, fraction_threshold=0.6, max_reads=30):
        read_lengths = []
        
        # Open file (gzipped or not)
        open_func = gzip.open if file_path.endswith('.gz') else open
        try:
            with open_func(file_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    read_lengths.append(len(record.seq))
                    if len(read_lengths) >= max_reads:
                        break
        except EOFError:
            # Handle the case where the file is empty or corrupted
            return "Unknown (File error)"
        
        if not read_lengths:
            return "Unknown (No reads)"
        
        # Handle the case of a single read
        if len(read_lengths) == 1:
            return "Looks like WGS" if read_lengths[0] >= length_threshold else "Looks like 16S"
        
        read_lengths = np.array(read_lengths)
        fraction_below_threshold = np.sum(read_lengths < length_threshold) / len(read_lengths)

        if fraction_below_threshold >= fraction_threshold:
            return "Looks like 16S"
        else:
            return "Looks like WGS"

    def add_sample(self, sample_name, file_paths):
        if sample_name not in self.samples:
            self.samples[sample_name] = file_paths
            read_types = []
            import random
            for file_path in random.sample(file_paths, min(5, len(file_paths))):  # Use 5 random files
                read_types.append(self.detect_read_type(file_path))
            if read_types:
                read_type = max(set(read_types), key=read_types.count)  # Majority vote
            else:
                read_type = "Unknown"
            
            self.file_tree.insert("", "end", sample_name, text=sample_name, values=(f"{len(file_paths)} files", "New", read_type))
        else:
            new_files = [f for f in file_paths if f not in self.samples[sample_name]]
            self.samples[sample_name].extend(new_files)
            current_status = self.file_tree.item(sample_name, "values")[1]
            current_read_type = self.file_tree.item(sample_name, "values")[2]
            if new_files and current_status != "New":
                new_status = "New"
            else:
                new_status = current_status
            self.file_tree.item(sample_name, values=(f"{len(self.samples[sample_name])} files", new_status, current_read_type))

        for file_path in file_paths:
            if self.file_tree.exists(file_path):
                continue
            self.file_tree.insert(sample_name, "end", file_path, text=os.path.basename(file_path), values=("", ""))

        if self.auto_analyze_var.get() and new_status == "New":
            self.run_analysis_for_sample(sample_name)

    def run_all_analyses(self):
        for sample in self.samples:
            self.run_analysis_for_sample(sample)

    def run_analysis_for_sample(self, sample_name):
        with open(self.config_file, 'r') as f:
            config = json.load(f)
        analysis_type = config.get('analysis_type', 'default_analysis_type')
        
        with open(self.config_file, 'r') as f:
            config = json.load(f)
        
        # Add conda environment to PATH
        # conda_env_path = os.path.expanduser("~/miniconda3/envs/fastcat/bin")
        # os.environ["PATH"] = f"{conda_env_path}:{os.environ['PATH']}"
        
        cmd_runner = MMonitorCMD()
        args = [
            "-a", analysis_type,
            "-c", self.config_file,
            "-i"] + self.samples[sample_name] + [
            "-s", sample_name,
            "-d", self.selected_date.strftime("%Y-%m-%d"),
            "-p", self.project_entry.get(),
            "-u", self.subproject_entry.get(),
            "-n", config.get('min_abundance', '0.01'),
            "-t", config.get('threads', multiprocessing.cpu_count()),
            "--overwrite"
        ]

        if analysis_type == "assembly-functional":
            args.extend([
                "--assembly-mode", config.get('assembly_mode', 'nano-raw'),
                "--medaka-model", config.get('medaka_model', 'r1041_e82_400bps_sup_v5.0.0'),
            ])
            if config.get('is_isolate', False):
                args.append("--isolate")

        cmd_runner.initialize_from_args(cmd_runner.parse_arguments(args))
        try:
            cmd_runner.run()
            self.file_tree.item(sample_name, values=(f"{len(self.samples[sample_name])} files", "Completed", self.file_tree.item(sample_name, "values")[2]))
        except Exception as e:
            print(f"Error running analysis for {sample_name}: {str(e)}")
            self.file_tree.item(sample_name, values=(f"{len(self.samples[sample_name])} files", "Error", self.file_tree.item(sample_name, "values")[2]))

    def suggest_sample_name(self, path):
        # Try to extract meaningful sample name from the path
        path_parts = path.split(os.sep)
        for part in reversed(path_parts):
            if re.match(r'barcode\d+', part, re.IGNORECASE):
                return part
            if re.match(r'sample\d+', part, re.IGNORECASE):
                return part
        # If no meaningful name found, use the last directory name
        return os.path.basename(path)

    def show_calendar(self):
        """Show calendar popup for date selection"""
        top = tk.Toplevel(self)
        top.title("Select Date")
        
        cal = Calendar(top, selectmode='day', date_pattern='yyyy-mm-dd',
                      year=self.selected_date.year,
                      month=self.selected_date.month,
                      day=self.selected_date.day)
        cal.pack(padx=10, pady=10)
        
        def set_date():
            date_str = cal.get_date()
            self.selected_date = datetime.datetime.strptime(date_str, '%Y-%m-%d').date()
            self.date_entry.delete(0, tk.END)
            self.date_entry.insert(0, date_str)
            top.destroy()
        
        ttk.Button(top, text="OK", command=set_date).pack(pady=5)

    def toggle_sample_name_entry(self):
        if self.use_suggested_name_var.get():
            self.sample_name_entry.configure(state="disabled")
        else:
            self.sample_name_entry.configure(state="normal")

    def update_appearance(self):
        """Update appearance based on current theme"""
        is_dark = ctk.get_appearance_mode() == "dark"
        
        # Update background colors
        bg_color = "#1a1a1a" if is_dark else "#ffffff"
        frame_color = "#2d2d2d" if is_dark else "#f5f5f5"
        text_color = "#ffffff" if is_dark else "#000000"
        
        # Update main container
        self.configure(fg_color=bg_color)
        
        # Update all frames
        for widget in self.winfo_children():
            if isinstance(widget, ctk.CTkFrame):
                widget.configure(fg_color=frame_color)
        
        # Update labels and buttons
        for widget in self.winfo_children():
            if isinstance(widget, ctk.CTkLabel):
                widget.configure(text_color=text_color)
            elif isinstance(widget, ctk.CTkButton):
                widget.configure(
                    fg_color=frame_color,
                    text_color=text_color,
                    hover_color="#3d3d3d" if is_dark else "#e0e0e0"
                )
        
        # Update treeview colors
        self.configure_treeview_colors()

    def check_config(self):
        if not self.config:
            messagebox.showwarning("Invalid Configuration", "Please set up analysis parameters in the Analysis tab first.")
            return False
        return True

    def update_queue_display(self):
        """Update the queue tree view"""
        try:
            self.queue_tree.delete(*self.queue_tree.get_children())
            
            for sample_name, data in self.queue_samples.items():
                self.queue_tree.insert("", "end", text=sample_name,
                                     values=(f"{len(data['files'])} files",
                                            data['first_seen']))
            
            self.update_queue_count()
        except Exception as e:
            print(f"Error updating queue display: {e}")

    def process_queue(self):
        """Only update trees if there are changes"""
        current_time = time.time()
        if current_time - self.last_update_time > 1:
            if self.needs_update:  # Only update if there are changes
                self.update_file_tree()
                self.update_queue_display()
                self.needs_update = False
            self.last_update_time = current_time
        self.after(100, self.process_queue)

    def run_selected_sample(self):
        """Run analysis for selected samples in the file tree using threading"""
        selected_items = self.file_tree.selection()
        if not selected_items:
            messagebox.showinfo("Info", "Please select a sample to analyze.")
            return
        
        # Start a thread for each selected sample
        for item in selected_items:
            # Get the parent item if this is a file
            parent = self.file_tree.parent(item)
            sample_name = self.file_tree.item(item)["text"] if not parent else self.file_tree.item(parent)["text"]
            
            if sample_name in self.all_samples:
                # Create and start a thread for this sample
                analysis_thread = threading.Thread(
                    target=self._run_analysis_thread,
                    args=(sample_name, item),
                    daemon=True
                )
                analysis_thread.start()
            else:
                print(f"Error: Sample {sample_name} not found in samples list")
                messagebox.showerror("Error", f"Sample {sample_name} not found in samples list")

    def _run_analysis_thread(self, sample_name, tree_item):
        """Thread function to run analysis for a single sample"""
        try:
            print(f"Starting analysis thread for sample: {sample_name}")
            sample_data = self.all_samples[sample_name]
            sample_files = sample_data["files"]
            
            # Update GUI in main thread
            self.after(0, lambda: self.file_tree.item(
                tree_item, 
                values=(f"{len(sample_files)} files", "Running...")
            ))
            
            # Create MMonitorCMD instance
            cmd_runner = MMonitorCMD()
            cmd_runner.offline_mode = self.offline_mode
            cmd_runner.logged_in = self.logged_in
            cmd_runner.current_user = self.current_user
            
            # Initialize Django DB interface
            db_interface = DjangoDBInterface(self.config_file)
            db_interface.offline_mode = self.offline_mode
            db_interface.username = self.current_user
            
            # Get analysis type from config
            analysis_type = self.config.get('analysis_type', 'taxonomy-wgs')
            print(f"Starting {analysis_type} analysis for {sample_name}")
            
            # Run analysis and get output files
            success = False
            try:
                # Create temporary concatenated file
                with tempfile.NamedTemporaryFile(delete=False, suffix='.fastq.gz') as temp_file:
                    print(f"Creating concatenated file: {temp_file.name}")
                    # Concatenate input files
                    cmd_runner.centrifuger_runner.concatenate_gzipped_files(sample_files, temp_file.name)
                    
                    # Run Centrifuger analysis
                    success = cmd_runner.centrifuger_runner.run_centrifuger(
                        input_file=temp_file.name,
                        sample_name=sample_name,
                        db_path=self.config['centrifuger_db']
                    )
                    
                    if success:
                        # Get the report file path
                        report_file = os.path.join(
                            cmd_runner.pipeline_out,
                            sample_name,
                            "centrifuger_report.tsv"
                        )
                        
                        # Send results to Django DB
                        if os.path.exists(report_file):
                            print(f"Uploading results for sample {sample_name} to database")
                            db_interface.send_nanopore_record_centrifuge(
                                kraken_out_path=report_file,
                                sample_name=sample_name,
                                project_id=self.project_entry.get() or "default_project",
                                subproject_id=self.subproject_entry.get() or "default_subproject",
                                date=self.selected_date.strftime('%Y-%m-%d'),
                                overwrite=True  # Since we're running a new analysis
                            )
                        else:
                            print(f"Report file not found: {report_file}")
                            success = False
                    
                    # Clean up temporary file
                    os.unlink(temp_file.name)
            
            except Exception as e:
                print(f"Error during analysis: {e}")
                success = False
            
            # Update GUI in main thread
            if success:
                self.after(0, lambda: [
                    self.file_tree.item(tree_item, values=(f"{len(sample_files)} files", "Completed")),
                    messagebox.showinfo("Success", f"Analysis completed for sample {sample_name}")
                ])
            else:
                self.after(0, lambda: [
                    self.file_tree.item(tree_item, values=(f"{len(sample_files)} files", "Failed")),
                    messagebox.showerror("Error", f"Analysis failed for sample {sample_name}")
                ])
            
            # Remove from queue if present
            if sample_name in self.queue_samples:
                self.after(0, lambda: [
                    self.queue_samples.pop(sample_name, None),
                    self.update_queue_display()
                ])
            
        except Exception as e:
            error_msg = f"Error running analysis for {sample_name}: {str(e)}"
            print(error_msg)
            traceback.print_exc()
            self.after(0, lambda: [
                self.file_tree.item(tree_item, values=(f"{len(sample_files)} files", "Error")),
                messagebox.showerror("Error", error_msg)
            ])

    def analyze_selected_queue(self):
        selected_items = self.queue_tree.selection()
        if not selected_items:
            messagebox.showinfo("Info", "Please select samples from the queue to analyze.")
            return
            
        for item in selected_items:
            sample_name = self.queue_tree.item(item)["text"]
            self.run_analysis_for_sample(sample_name)
            self.queue_tree.delete(item)
            # Remove from pending files
            self.pending_files = [(f, s) for f, s in self.pending_files if s != sample_name]
        
        self.update_queue_count()

    def clear_queue(self):
        if messagebox.askyesno("Clear Queue", "Are you sure you want to clear the analysis queue?"):
            self.queue_tree.delete(*self.queue_tree.get_children())
            self.queue_samples.clear()  # Clear queue samples instead of pending_files
            self.update_queue_count()

    def update_queue_count(self):
        sample_count = len(self.queue_samples)
        file_count = sum(len(data["files"]) for data in self.queue_samples.values())
        self.queue_count.configure(text=f"({sample_count} sample{'s' if sample_count != 1 else ''}, "
                                      f"{file_count} file{'s' if file_count != 1 else ''})")

    def configure_treeview_colors(self):
        """Configure treeview colors for light theme"""
        style = ttk.Style()
        
        # Configure Treeview style
        style.configure(
            "Custom.Treeview",
            background="white",
            fieldbackground="white",
            foreground="black",
            font=("Helvetica", 11)
        )
        
        # Configure Treeview headers
        style.configure(
            "Custom.Treeview.Heading",
            background="#F5F5F5", 
            foreground="black",
            font=("Helvetica", 11, "bold")
        )
        
        # Configure selection colors
        style.map(
            "Custom.Treeview",
            background=[("selected", "#E8F0FE")],  # Light blue selection
            foreground=[("selected", "#000000")]   # Keep text black when selected
        )
        
        # Update both trees
        for tree in [self.file_tree, self.queue_tree]:
            tree.configure(style="Custom.Treeview")
            tree.tag_configure("new", foreground="#28a745")  # Green for new items

    def store_expanded_state(self):
        """Store which items are currently expanded"""
        self.expanded_items = {
            item for item in self.file_tree.get_children() 
            if self.file_tree.item(item, 'open')
        }

    def restore_expanded_state(self):
        """Restore previously expanded items"""
        for item in self.expanded_items:
            if self.file_tree.exists(item):
                self.file_tree.item(item, open=True)

    def on_closing(self):
        """Handle window closing event"""
        if self.watching:
            if messagebox.askyesno("Confirm Exit", "Folder watching is active. Do you want to stop watching and close?"):
                self.stop_watching()
                self.parent.destroy()
        else:
            self.parent.destroy()

    def load_from_csv(self):
        """Load multiple samples from a CSV file"""
        file_path = filedialog.askopenfilename(
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")],
            title="Select CSV file with sample information"
        )

        if not file_path:
            print("No file selected.")
            return

        # Check if CSV is empty
        if os.path.getsize(file_path) == 0:
            print("Selected CSV file is empty.")
            return

        error_messages = []  # list to accumulate error messages

        try:
            with open(file_path, 'r') as file:
                reader = csv.DictReader(file)
                
                # Verify required columns
                required_columns = ["sample_name", "date", "project_name", "subproject_name", "sample folder"]
                missing_columns = [col for col in required_columns if col not in reader.fieldnames]
                if missing_columns:
                    messagebox.showerror("Error", f"Missing required columns in CSV: {', '.join(missing_columns)}")
                    return

                for row in reader:
                    # Check if provided path exists
                    if not os.path.exists(row["sample folder"].strip()):
                        error_message = f"Invalid path from CSV: {row['sample folder'].strip()}"
                        print(error_message)
                        error_messages.append(error_message)
                        continue

                    # Look for the fastq_pass folder in the provided path and its child directories
                    folder_path = None
                    for root, dirs, files in os.walk(row["sample folder"].strip()):
                        if "fastq_pass" in dirs:
                            folder_path = os.path.join(root, "fastq_pass")
                            break

                    if not folder_path:
                        error_message = f"'fastq_pass' directory not found for path: {row['sample folder'].strip()}"
                        print(error_message)
                        error_messages.append(error_message)
                        continue

                    # If multiplexing is selected, navigate further to the barcode_x folder
                    if hasattr(self, 'use_multiplexing') and self.use_multiplexing.get():
                        barcode_id_string = str(row.get('Barcode ID', ''))
                        if barcode_id_string:
                            if len(barcode_id_string) == 1:
                                barcode_id_string = f"0{barcode_id_string}"
                            barcode_folder = f"barcode{barcode_id_string}"
                            folder_path = os.path.join(folder_path, barcode_folder)

                            if not os.path.exists(folder_path):
                                error_message = f"Barcode folder '{barcode_folder}' not found."
                                print(error_message)
                                error_messages.append(error_message)
                                continue

                    # Get files from the folder
                    files = get_files_from_folder(folder_path)
                    if not files:
                        error_message = f"No valid sequence files found in {folder_path}"
                        print(error_message)
                        error_messages.append(error_message)
                        continue

                    # Add to all_samples and queue_samples
                    sample_name = row["sample_name"]
                    self.all_samples[sample_name] = {
                        "files": files,
                        "path": folder_path
                    }
                    
                    # Update project and subproject entries with the first row's values
                    if not self.project_entry.get():
                        self.project_entry.insert(0, row["project_name"])
                    if not self.subproject_entry.get():
                        self.subproject_entry.insert(0, row["subproject_name"])
                    
                    # Add to queue if it's a new sample
                    current_time = time.strftime("%H:%M:%S", time.localtime(time.time()))
                    self.queue_samples[sample_name] = {
                        "files": files,
                        "first_seen": current_time
                    }

            # Update displays
            self.needs_update = True
            
            # Show summary
            num_samples = len(self.all_samples)
            if num_samples > 0:
                messagebox.showinfo(
                    "Import Complete",
                    f"Successfully imported {num_samples} sample{'s' if num_samples != 1 else ''} from CSV."
                )
            
            # Show errors if any
            if error_messages:
                if len(error_messages) > 2:
                    error_text = "\n".join(error_messages[:2]) + f"\n...and {len(error_messages)-2} more errors"
                else:
                    error_text = "\n".join(error_messages)
                messagebox.showwarning("Import Warnings", error_text)

        except Exception as e:
            messagebox.showerror("Error", f"Error reading CSV file:\n{str(e)}")
            print(f"Error reading CSV: {e}")



















