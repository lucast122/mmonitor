import os
import datetime
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import traceback
import customtkinter as ctk

from mmonitor.userside.InputWindow import get_files_from_folder
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
import argparse
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


class FolderWatcherWindow(ctk.CTkFrame):
    def __init__(self, parent, gui_ref):
        print("Initializing FolderWatcherWindow")
        super().__init__(parent)
        self.parent = parent
        self.gui = gui_ref
        
        # Get authentication state from main window
        self.offline_mode = gui_ref.offline_mode
        self.logged_in = gui_ref.logged_in
        self.current_user = gui_ref.current_user
        
        print(f"Authentication state: offline={self.offline_mode}, logged_in={self.logged_in}, user={self.current_user}")
        
        # Initialize variables
        self.folder_var = tk.StringVar()
        self.auto_analyze_var = tk.BooleanVar(value=False)
        self.project_var = tk.StringVar()
        self.subproject_var = tk.StringVar()
        self.date_var = tk.StringVar(value=datetime.datetime.now().strftime("%Y-%m-%d"))
        
        self.samples = {}
        self.watching = False
        self.folder_monitor = None
        
        self.config = {}
        self.watch_start_time = None
        self.all_samples = {}
        self.queue_samples = {}
        self.pending_files = []
        self.config_file = os.path.join(ROOT, "src", "resources", "pipeline_config.json")
        
        # Load configuration first
        self.load_config()
        if not self.config:
            messagebox.showwarning(
                "Configuration Missing",
                "No configuration found. Please configure the pipeline first."
            )
        
        self.create_widgets()
        
        # Add thread tracking
        self.running_analyses = {}  # Keep track of running analysis threads
        self.analysis_lock = threading.Lock()  # Thread safety for analysis operations

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
        """Create all widgets for the folder watcher window"""
        self.main_container = ctk.CTkFrame(self, fg_color="transparent")
        self.main_container.pack(fill="both", expand=True, padx=10, pady=5)
        
        # Create folder selection frame
        folder_frame = ctk.CTkFrame(self.main_container, fg_color="transparent")
        folder_frame.pack(fill="x", pady=5)
        
        # Project info frame
        project_frame = ctk.CTkFrame(self.main_container, fg_color="transparent")
        project_frame.pack(fill="x", pady=5)
        
        # Project fields
        ctk.CTkLabel(project_frame, text="Project:").pack(side="left", padx=5)
        self.project_entry = ctk.CTkEntry(project_frame, width=200)
        self.project_entry.pack(side="left", padx=5)
        
        ctk.CTkLabel(project_frame, text="Subproject:").pack(side="left", padx=5)
        self.subproject_entry = ctk.CTkEntry(project_frame, width=200)
        self.subproject_entry.pack(side="left", padx=5)
        
        ctk.CTkLabel(project_frame, text="Date:").pack(side="left", padx=5)
        self.date_entry = ctk.CTkEntry(project_frame, width=200)
        self.date_entry.pack(side="left", padx=5)
        self.date_entry.insert(0, datetime.datetime.now().strftime("%Y-%m-%d"))
        
        # Folder selection
        ctk.CTkLabel(folder_frame, text="Watch Folder:").pack(side="left", padx=5)
        self.folder_entry = ctk.CTkEntry(folder_frame, textvariable=self.folder_var, width=400)
        self.folder_entry.pack(side="left", padx=5)
        
        # Store reference to browse button
        self.browse_button = ctk.CTkButton(folder_frame, text="Browse", 
                                          command=self.browse_folder)
        self.browse_button.pack(side="left", padx=5)
        
        # Watch control buttons
        self.watch_button = ctk.CTkButton(folder_frame, text="Start Watching",
                                         command=self.toggle_watching)
        self.watch_button.pack(side="left", padx=5)
        
        # Create trees
        self.create_file_trees(self.main_container)
        
        # Create control buttons frame at the bottom
        control_frame = ctk.CTkFrame(self.main_container, fg_color="transparent")
        control_frame.pack(fill="x", pady=5)
        
        # Left side - Analyze Selected button
        left_controls = ctk.CTkFrame(control_frame, fg_color="transparent")
        left_controls.pack(side="left", padx=5)
        
        self.analyze_button = ctk.CTkButton(
            left_controls,
            text="Analyze Selected",
            command=self.analyze_selected
        )
        self.analyze_button.pack(side="left", padx=5)
        
        # Right side - Auto-analyze checkbox and Run Queue button
        right_controls = ctk.CTkFrame(control_frame, fg_color="transparent")
        right_controls.pack(side="right", padx=5)
        
        self.auto_analyze_checkbox = ctk.CTkCheckBox(
            right_controls,
            text="Auto-analyze new samples",
            variable=self.auto_analyze_var,
            onvalue=True,
            offvalue=False
        )
        self.auto_analyze_checkbox.pack(side="left", padx=5)
        
        self.run_queue_button = ctk.CTkButton(
            right_controls,
            text="Run Queue",
            command=self.run_queue
        )
        self.run_queue_button.pack(side="left", padx=5)

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
        # Project, subproject and date in one line
        info_frame = ctk.CTkFrame(top_frame, fg_color="transparent")
        info_frame.pack(fill="x", pady=(10, 0))

        # Initialize date variables
        self.selected_date = datetime.date.today()
        self.use_file_date = tk.BooleanVar(value=True)

        # Project name
        ctk.CTkLabel(info_frame, text="Project:", 
                     font=("Helvetica", 12, "bold"), 
                     width=60).pack(side="left")
        self.project_entry = ctk.CTkEntry(info_frame, height=32, width=150)
        self.project_entry.pack(side="left", padx=5)

        # Subproject
        ctk.CTkLabel(info_frame, text="Subproject:", 
                     font=("Helvetica", 12, "bold"), 
                     width=80).pack(side="left", padx=(10,0))
        self.subproject_entry = ctk.CTkEntry(info_frame, height=32, width=150)
        self.subproject_entry.pack(side="left", padx=5)

        # Date picker
        ctk.CTkLabel(info_frame, text="Date:", 
                     font=("Helvetica", 12, "bold"), 
                     width=40).pack(side="left", padx=(10,0))
        self.date_entry = ctk.CTkEntry(info_frame, height=32, width=100)
        self.date_entry.pack(side="left", padx=5)
        self.date_entry.insert(0, self.selected_date.strftime('%Y-%m-%d'))

        calendar_button = ctk.CTkButton(
            info_frame,
            text="📅", 
            width=32,
            height=32,
            command=self.show_calendar
        )
        calendar_button.pack(side="left")

        self.use_file_date_checkbox = ctk.CTkCheckBox(
            info_frame,
            text="Use file date",
            variable=self.use_file_date,
            font=("Helvetica", 11)
        )
        self.use_file_date_checkbox.pack(side="left", padx=10)

    def create_file_trees(self, container):
        # Create a frame for both trees with proper spacing
        trees_frame = ctk.CTkFrame(container, fg_color="transparent")
        trees_frame.pack(fill="both", expand=True, pady=5, padx=20)
        
        # Left tree (70% width)
        left_frame = ctk.CTkFrame(trees_frame, fg_color="transparent")
        left_frame.pack(side="left", fill="both", expand=True, padx=(0, 10))
        
        # Title for left tree
        left_header = ctk.CTkFrame(left_frame, fg_color="transparent")
        left_header.pack(fill="x", pady=2)
        ctk.CTkLabel(left_header, text="Detected Files", 
                     font=("Helvetica", 14, "bold")).pack(anchor="w")
        
        # Create container frame
        tree_container = ctk.CTkFrame(left_frame, fg_color="#dbdbdb")
        tree_container.pack(fill="both", expand=True)
        
        # Create customtkinter scrollbar
        y_scrollbar = ctk.CTkScrollbar(tree_container)
        y_scrollbar.pack(side="right", fill="y")
        
        # Create the file tree
        self.file_tree = ttk.Treeview(
            tree_container,
            columns=("files", "status"),
            selectmode="extended",
            height=15,
            yscrollcommand=y_scrollbar.set,
            style="Treeview"
        )
        self.file_tree.pack(side="left", fill="both", expand=True)
        
        # Configure scrollbar
        y_scrollbar.configure(command=self.file_tree.yview)
        
        # Configure file tree columns
        self.file_tree.heading("#0", text="Sample Name")
        self.file_tree.heading("files", text="Files")
        self.file_tree.heading("status", text="Status")
        
        self.file_tree.column("#0", width=350, minwidth=350)
        self.file_tree.column("files", width=100, minwidth=100)
        self.file_tree.column("status", width=100, minwidth=100)
        
        # Right tree (30% width)
        right_frame = ctk.CTkFrame(trees_frame, fg_color="transparent")
        right_frame.pack(side="right", fill="both", expand=True)
        
        # Title for right tree
        right_header = ctk.CTkFrame(right_frame, fg_color="transparent")
        right_header.pack(fill="x", pady=2)
        ctk.CTkLabel(right_header, text="Analysis Queue", 
                     font=("Helvetica", 14, "bold")).pack(side="left")
        self.queue_count = ctk.CTkLabel(right_header, text="(0 samples, 0 files)", 
                                       font=("Helvetica", 12))
        self.queue_count.pack(side="left", padx=5)
        
        # Create container frame for right tree
        right_tree_container = ctk.CTkFrame(right_frame, fg_color="#dbdbdb")
        right_tree_container.pack(fill="both", expand=True)
        
        # Create customtkinter scrollbar for right tree
        right_y_scrollbar = ctk.CTkScrollbar(right_tree_container)
        right_y_scrollbar.pack(side="right", fill="y")
        
        # Create queue tree
        self.queue_tree = ttk.Treeview(
            right_tree_container,
            columns=("files", "first_seen"),
            selectmode="extended",
            height=15,
            yscrollcommand=right_y_scrollbar.set,
            style="Treeview"
        )
        self.queue_tree.pack(side="left", fill="both", expand=True)
        
        # Configure right scrollbar
        right_y_scrollbar.configure(command=self.queue_tree.yview)
        
        # Configure queue tree columns
        self.queue_tree.heading("#0", text="Sample Name")
        self.queue_tree.heading("files", text="Files")
        self.queue_tree.heading("first_seen", text="Added")
        
        self.queue_tree.column("#0", width=200, minwidth=150)
        self.queue_tree.column("files", width=100, minwidth=80)
        self.queue_tree.column("first_seen", width=120, minwidth=100)
        
        # Apply color configuration
        self.configure_treeview_colors()

    def configure_file_tree(self, tree):
        """Configure the file tree columns and headings"""
        tree.configure(show="headings")  # Hide the tree view border
        
        # Configure columns
        tree.heading("files", text="Files")
        tree.heading("status", text="Status")
        
        # Set column widths
        tree.column("files", width=100, minwidth=80)
        tree.column("status", width=150, minwidth=100)

    def configure_queue_tree(self, tree):
        """Configure the queue tree columns and headings"""
        tree.configure(show="headings")  # Hide the tree view border
        
        # Configure columns
        tree.heading("files", text="Files")
        tree.heading("first_seen", text="Added")
        
        # Set column widths
        tree.column("files", width=100, minwidth=80)
        tree.column("first_seen", width=120, minwidth=100)

    def add_scrollbars(self, container, tree):
        """Add scrollbars to a treeview"""
        # Create frame for scrollbars
        scroll_frame = ttk.Frame(container)
        scroll_frame.pack(fill="both", expand=True)
        
        # Add vertical scrollbar (removed width parameter)
        vsb = ttk.Scrollbar(scroll_frame, orient="vertical", command=tree.yview)
        vsb.pack(side="right", fill="y")
        
        # Add horizontal scrollbar
        hsb = ttk.Scrollbar(scroll_frame, orient="horizontal", command=tree.xview)
        hsb.pack(side="bottom", fill="x")
        
        # Configure tree to use scrollbars
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.pack(in_=scroll_frame, side="left", fill="both", expand=True)

    def create_bottom_controls(self, container):
        """Create bottom controls with improved spacing"""
        bottom_frame = ctk.CTkFrame(container, fg_color="transparent")
        bottom_frame.pack(fill="x", padx=20, pady=(0, 20))  # Added bottom padding
        
        # Left side controls
        left_controls = ctk.CTkFrame(bottom_frame, fg_color="transparent")
        left_controls.pack(side="left", fill="y")
        
        self.auto_analyze_checkbox = ctk.CTkCheckBox(
            left_controls,
            text="Auto-analyze new samples",
            variable=self.auto_analyze_var,
            font=("Helvetica", 12)
        )
        self.auto_analyze_checkbox.pack(side="left", padx=(0, 10))
        
        # Right side controls with spacing
        right_controls = ctk.CTkFrame(bottom_frame, fg_color="transparent")
        right_controls.pack(side="right", fill="y", padx=(0, 20))  # Added right padding
        
        self.analyze_button = ctk.CTkButton(
            right_controls,
            text="Analyze Selected",
            command=self.analyze_selected_queue,
            height=32,
            font=("Helvetica", 12)
        )
        self.analyze_button.pack(side="left", padx=5)
        
        self.clear_button = ctk.CTkButton(
            right_controls,
            text="Clear Queue",
            command=self.clear_queue,
            height=32,
            font=("Helvetica", 12)
        )
        self.clear_button.pack(side="left", padx=5)

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
        """Update the file tree view with proper error handling"""
        try:
            # Clear existing items
            self.file_tree.delete(*self.file_tree.get_children())
            
            print(f"Updating tree with {len(self.all_samples)} samples")  # Debug print
            
            for sample_name, data in self.all_samples.items():
                try:
                    # Insert sample as parent
                    files = data.get('files', [])
                    status = data.get('status', 'Pending')
                    
                    sample_id = self.file_tree.insert(
                        "", "end",
                        text=sample_name,
                        values=(f"{len(files)} files", status))
                    
                    # Insert each file as child
                    for file_path in files:
                        try:
                            self.file_tree.insert(
                                sample_id, "end",
                                text=os.path.basename(file_path),
                                values=("", "")
                            )
                        except Exception as e:
                            print(f"Error inserting file {file_path}: {e}")
                
                except Exception as e:
                    print(f"Error processing sample {sample_name}: {e}")
            
            # Update queue count
            queue_count = len(self.queue_samples)
            total_files = sum(len(data.get('files', [])) for data in self.queue_samples.values())
            self.queue_count.configure(text=f"({queue_count} samples, {total_files} files)")
            
        except Exception as e:
            print(f"Error updating file tree: {e}")
            import traceback
            traceback.print_exc()

    
    def analysis_completed(self, sample_name, success):
        """Handle analysis completion in main thread"""
        if success:
            self.update_sample_status(sample_name, "Completed")
            print(f"Analysis completed successfully for {sample_name}")
        else:
            self.update_sample_status(sample_name, "Failed")
            print(f"Analysis failed for {sample_name}")
        
        # Remove from queue if present
        if sample_name in self.queue_samples:
            del self.queue_samples[sample_name]
            self.needs_update = True

    def analysis_failed(self, sample_name, error_msg):
        """Handle analysis failure in main thread"""
        self.update_sample_status(sample_name, "Failed")
        messagebox.showerror("Analysis Error", f"Analysis failed for {sample_name}: {error_msg}")

    def run_all_analyses(self):
        """Run analyses for all samples in separate threads"""
        for sample_name in list(self.samples.keys()):
            if sample_name not in self.running_analyses:
                self.run_analysis_for_sample(sample_name)
            else:
                print(f"Analysis already running for {sample_name}")

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
        """Stop watching and clean up running analyses"""
        if self.folder_monitor:
            self.folder_monitor.stop()
            self.folder_monitor = None
        
        # Wait for running analyses to complete
        for thread in self.running_analyses.values():
            thread.join(timeout=0.1)  # Short timeout to avoid blocking
        
        self.watching = False
        self.watch_button.configure(text="Start Watching")
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
        """Process the update queue with debug output"""
        if self.needs_update:
            print("Processing update queue...")  # Debug print
            self.update_file_tree()
            self.update_queue_tree()
            self.needs_update = False
        
        # Schedule next update
        self.after(100, self.process_queue)


    def _run_analysis_thread(self, sample_name, tree_item):
        """Thread function to run analysis for a single sample"""
        try:
            print(f"Starting analysis thread for sample: {sample_name}")
            
            # Update status in GUI thread
            self.after(0, lambda: self.update_sample_status(sample_name, "Analyzing..."))
            
            # Get sample files
            files = self.all_samples[sample_name]['files']
            
            # Get project information
            project_info = {
                'project': self.project_entry.get() or "default_project",
                'subproject': self.subproject_entry.get() or "default_subproject",
                'date': self.date_entry.get()
            }
            
            # Initialize MMonitorCMD
            cmd_runner = MMonitorCMD()
            cmd_runner.offline_mode = self.offline_mode
            cmd_runner.logged_in = self.logged_in
            cmd_runner.current_user = self.current_user
            
            # Load config file
            try:
                with open(self.config_file, 'r') as f:
                    config = json.load(f)
                    print(f"Loaded config: {config}")
            except Exception as e:
                print(f"Error loading config: {e}")
                raise
            
            # Verify database paths exist
            centrifuger_db = config.get('centrifuger_db')
            if not centrifuger_db:
                raise ValueError("Centrifuger database path not found in config")
            
            if not os.path.exists(os.path.dirname(centrifuger_db)):
                raise ValueError(f"Centrifuger database directory not found: {centrifuger_db}")
            
            # Set up arguments with proper database paths
            args = argparse.Namespace(
                analysis=config.get('analysis_type', 'taxonomy-wgs'),
                config=self.config_file,
                threads=config.get('threads', 4),
                sample=sample_name,
                project=project_info['project'],
                subproject=project_info['subproject'],
                date=datetime.datetime.strptime(project_info['date'], '%Y-%m-%d').date(),
                input=files,
                multicsv=None,
                centrifuger_db=centrifuger_db,  # Use verified path
                emu_db=config.get('emu_db'),
                min_length=config.get('min_length', 1000),
                min_quality=config.get('min_quality', 10.0),
                min_abundance=config.get('min_abundance', 0.01),
                overwrite=True,
                qc=True,
                update=False,
                verbose=True,
                loglevel='INFO'
            )
            
            print(f"Analysis arguments:")
            print(f"Centrifuger DB: {args.centrifuger_db}")
            print(f"EMU DB: {args.emu_db}")
            print(f"Sample: {args.sample}")
            print(f"Project: {args.project}")
            print(f"Subproject: {args.subproject}")
            
            # Initialize and run analysis
            cmd_runner.initialize_from_args(args)
            success = cmd_runner.run()
            
            # Update GUI in main thread
            def update_gui():
                if success:
                    self.update_sample_status(sample_name, "Completed")
                    messagebox.showinfo("Success", f"Analysis completed for sample {sample_name}")
                else:
                    self.update_sample_status(sample_name, "Failed")
                    messagebox.showerror("Error", f"Analysis failed for sample {sample_name}")
                
                # Remove from queue if present
                if sample_name in self.queue_samples:
                    del self.queue_samples[sample_name]
                    self.update_queue_display()
            
            self.after(0, update_gui)
            
        except Exception as e:
            error_msg = f"Error running analysis for {sample_name}: {str(e)}"
            print(error_msg)
            traceback.print_exc()
            
            # Update GUI in main thread
            self.after(0, lambda: [
                self.update_sample_status(sample_name, "Error"),
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

    def configure_treeview_colors(self):
        """Configure cross-platform consistent treeview styling"""
        style = ttk.Style()
        
        # Define colors
        bg_color = "#dbdbdb"      # Light grey background
        selected_bg = "#b4b4b4"   # Slightly darker grey for selection
        
        # Configure Treeview
        style.configure(
            "Treeview",
            background=bg_color,
            fieldbackground=bg_color,
            foreground="black",
            borderwidth=0,
            font=("Helvetica", 11),
            rowheight=16
        )
        
        # Configure Treeview headers
        style.configure(
            "Treeview.Heading",
            background=bg_color,
            foreground="black",
            borderwidth=0,
            font=("Helvetica", 11, "bold")
        )
        
        # Configure selection colors
        style.map(
            "Treeview",
            background=[("selected", selected_bg)],
            foreground=[("selected", "black")]
        )

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

    def update_queue_tree(self):
        """Update the queue tree view"""
        try:
            # Clear existing items
            self.queue_tree.delete(*self.queue_tree.get_children())
            
            for sample_name, data in self.queue_samples.items():
                try:
                    # Insert sample with file count and first seen time
                    files = data.get('files', [])
                    first_seen = data.get('first_seen', '')
                    
                    sample_id = self.queue_tree.insert(
                        "", "end",
                        text=sample_name,
                        values=(f"{len(files)} files", first_seen)
                    )
                    
                    # Insert each file as child
                    for file_path in files:
                        try:
                            self.queue_tree.insert(
                                sample_id, "end",
                                text=os.path.basename(file_path),
                                values=("", "")
                            )
                        except Exception as e:
                            print(f"Error inserting queue file {file_path}: {e}")
                
                except Exception as e:
                    print(f"Error processing queue sample {sample_name}: {e}")
            
            # Update queue count label
            queue_count = len(self.queue_samples)
            total_files = sum(len(data.get('files', [])) for data in self.queue_samples.values())
            self.queue_count.configure(text=f"({queue_count} samples, {total_files} files)")
            
        except Exception as e:
            print(f"Error updating queue tree: {e}")
            import traceback
            traceback.print_exc()

    def create_control_buttons(self, container):
        """Create control buttons with all required functionality"""
        button_frame = ctk.CTkFrame(container, fg_color="transparent")
        button_frame.pack(fill="x", pady=5)
        
        # Left side controls
        left_controls = ctk.CTkFrame(button_frame, fg_color="transparent")
        left_controls.pack(side="left", padx=5)
        
        # Analyze Selected button
        self.analyze_button = ctk.CTkButton(
            left_controls,
            text="Analyze Selected",
            command=self.analyze_selected
        )
        self.analyze_button.pack(side="left", padx=5)
        create_tooltip(self.analyze_button, "Analyze selected samples (Shift+Click for multiple)")
        
        # Right side controls
        right_controls = ctk.CTkFrame(button_frame, fg_color="transparent")
        right_controls.pack(side="right", padx=5)
        
        # Auto-analyze checkbox
        self.auto_analyze_checkbox = ctk.CTkCheckBox(
            right_controls,
            text="Auto-analyze new samples",
            variable=self.auto_analyze_var,
            onvalue=True,
            offvalue=False
        )
        self.auto_analyze_checkbox.pack(side="left", padx=5)
        
        # Run Queue button
        self.run_queue_button = ctk.CTkButton(
            right_controls,
            text="Run Queue",
            command=self.run_queue
        )
        self.run_queue_button.pack(side="left", padx=5)

    def analyze_selected(self):
        """Analyze selected samples from the file tree using threading"""
        selected_items = self.file_tree.selection()
        if not selected_items:
            messagebox.showwarning("No Selection", "Please select at least one sample to analyze")
            return
        
        # Create a thread for each selected sample
        for item in selected_items:
            # Get the parent item if a file is selected
            parent = self.file_tree.parent(item)
            sample_item = parent if parent else item
            sample_name = self.file_tree.item(sample_item)['text']
            
            if sample_name in self.all_samples:
                # Create and start analysis thread
                analysis_thread = threading.Thread(
                    target=self._analyze_sample_thread,
                    args=(sample_name, sample_item),
                    daemon=True
                )
                analysis_thread.start()

    def _analyze_sample_thread(self, sample_name, tree_item):
        """Thread function to run analysis for a single sample"""
        try:
            print(f"Starting analysis for sample: {sample_name}")
            
            # Update status in GUI thread
            self.after(0, lambda: self.update_sample_status(sample_name, "Analyzing..."))
            
            # Get sample files
            files = self.all_samples[sample_name]['files']
            
            # Get project information
            project_info = {
                'project': self.project_entry.get() or "default_project",
                'subproject': self.subproject_entry.get() or "default_subproject",
                'date': self.date_entry.get()
            }
            
            # Initialize MMonitorCMD
            cmd_runner = MMonitorCMD()
            cmd_runner.offline_mode = self.offline_mode
            cmd_runner.logged_in = self.logged_in
            cmd_runner.current_user = self.current_user
            
            # Load config file
            with open(self.config_file, 'r') as f:
                config = json.load(f)
            
            analysis_type = config.get('analysis_type', 'taxonomy-wgs')
            print(f"Analysis type from config: {analysis_type}")
            
            # Set up arguments with proper database paths
            args = argparse.Namespace(
                analysis=analysis_type,
                config=self.config_file,
                threads=config.get('threads', 4),
                sample=sample_name,
                project=project_info['project'],
                subproject=project_info['subproject'],
                date=datetime.datetime.strptime(project_info['date'], '%Y-%m-%d').date(),
                input=files,
                multicsv=None,
                centrifuger_db=os.path.abspath(config.get('centrifuger_db')),
                emu_db=os.path.abspath(config.get('emu_db')),
                min_length=config.get('min_length', 1000),
                min_quality=config.get('min_quality', 10.0),
                min_abundance=config.get('min_abundance', 0.01),
                overwrite=True,
                qc=True,
                update=False,
                verbose=True,
                loglevel='INFO'
            )
            
            # Initialize and run analysis
            cmd_runner.initialize_from_args(args)
            success = cmd_runner.run()
            
            if success:
                # Upload results based on analysis type
                if analysis_type == 'taxonomy-wgs':
                    # Path to centrifuger report
                    report_file = os.path.join(
                        cmd_runner.pipeline_out,
                        sample_name,
                        f"{sample_name}_centrifuger_report.tsv"
                    )
                    
                    # Upload centrifuger results
                    if os.path.exists(report_file):
                        cmd_runner.django_db.send_nanopore_record_centrifuge(
                            kraken_out_path=report_file,
                            sample_name=sample_name,
                            project_id=project_info['project'],
                            subproject_id=project_info['subproject'],
                            date=project_info['date'],
                            overwrite=True
                        )
                else:  # taxonomy-16s
                    # Path to EMU output
                    emu_out_path = os.path.join(cmd_runner.pipeline_out, sample_name)
                    
                    # Upload EMU results using the correct method name
                    cmd_runner.django_db.update_django_with_emu_out(
                        emu_out_path=emu_out_path,
                        tax_rank="species",
                        sample_name=sample_name,
                        project_name=project_info['project'],
                        sample_date=project_info['date'],
                        subproject_name=project_info['subproject'],
                        overwrite=True
                    )
            
            # Update GUI in main thread
            def update_gui():
                if success:
                    self.update_sample_status(sample_name, "Completed")
                    messagebox.showinfo("Success", f"Analysis completed for sample {sample_name}")
                else:
                    self.update_sample_status(sample_name, "Failed")
                    messagebox.showerror("Error", f"Analysis failed for sample {sample_name}")
                
                # Remove from queue if present
                if sample_name in self.queue_samples:
                    del self.queue_samples[sample_name]
                    self.update_queue_display()
            
            self.after(0, update_gui)
            
        except Exception as e:
            error_msg = f"Error analyzing sample {sample_name}: {str(e)}"
            print(error_msg)
            traceback.print_exc()
            
            # Update GUI in main thread
            self.after(0, lambda: [
                self.update_sample_status(sample_name, "Error"),
                messagebox.showerror("Error", error_msg)
            ])

    def run_queue(self):
        """Run all samples in the queue"""
        if not self.queue_samples:
            messagebox.showinfo("Queue Empty", "No samples in the queue to analyze")
            return
        
        try:
            # Get project information
            project_info = {
                'project': self.project_entry.get(),
                'subproject': self.subproject_entry.get(),
                'date': self.date_entry.get()
            }
            
            # Process each sample in the queue
            for sample_name in list(self.queue_samples.keys()):
                try:
                    print(f"Processing queued sample: {sample_name}")
                    
                    # Get all files for this sample (including new ones)
                    all_files = self.all_samples[sample_name]['files']
                    
                    # Update status
                    self.update_sample_status(sample_name, "Analyzing...")
                    
                    # Run the analysis
                    cmd_runner = MMonitorCMD(
                        files=all_files,
                        sample_name=sample_name,
                        config=self.config,
                        project_info=project_info
                    )
                    cmd_runner.run()
                    
                    # Update status and remove from queue
                    self.update_sample_status(sample_name, "Analyzed")
                    del self.queue_samples[sample_name]
                    
                except Exception as e:
                    print(f"Error processing queued sample {sample_name}: {e}")
                    self.update_sample_status(sample_name, "Failed")
                    messagebox.showerror("Queue Error", f"Error processing sample {sample_name}: {str(e)}")
            
            self.needs_update = True
            
        except Exception as e:
            print(f"Error in run_queue: {e}")
            messagebox.showerror("Error", f"Error processing queue: {str(e)}")

    def update_sample_status(self, sample_name, status):
        """Update the status of a sample in the file tree"""
        try:
            # Update status in all_samples
            if sample_name in self.all_samples:
                self.all_samples[sample_name]['status'] = status
            
            # Find and update the item in the tree
            for item in self.file_tree.get_children():
                if self.file_tree.item(item)['text'] == sample_name:
                    self.file_tree.set(item, "status", status)
                    break
            
        except Exception as e:
            print(f"Error updating sample status: {e}")

    def on_file_change(self, event):
        """Handle file system events"""
        try:
            if event.is_directory:
                return
                
            file_path = event.src_path
            if not any(file_path.endswith(ext) for ext in ['.fastq', '.fastq.gz', '.fq', '.fq.gz']):
                return
                
            sample_name = self.get_sample_name(file_path)
            if not sample_name:
                return
                
            # Add to all_samples
            if sample_name not in self.all_samples:
                self.all_samples[sample_name] = {'files': [], 'status': 'New'}
            if file_path not in self.all_samples[sample_name]['files']:
                self.all_samples[sample_name]['files'].append(file_path)
                
            # Add to queue if watching
            if self.watching and self.watch_start_time:
                file_time = datetime.datetime.fromtimestamp(os.path.getctime(file_path))
                if file_time > self.watch_start_time:
                    self.add_to_queue(sample_name)
                    if self.auto_analyze_var.get():
                        self.run_queue()
            
            self.needs_update = True
            
        except Exception as e:
            print(f"Error handling file change: {e}")



















