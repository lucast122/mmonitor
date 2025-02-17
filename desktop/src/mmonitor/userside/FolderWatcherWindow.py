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
import multiprocessing
import queue
import threading
from mmonitor.paths import ROOT, RESOURCES_DIR
import tkinter as tk
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
        self.analysis_type = tk.StringVar(value="taxonomy-wgs")
        
        # Initialize trace for analysis_type
        self.analysis_type_trace = self.analysis_type.trace_add('write', lambda *args: self.on_analysis_type_change())
        
        self.samples = {}
        self.watching = False
        self.folder_monitor = None
        
        self.config = {}
        self.watch_start_time = None
        self.all_samples = {}
        self.queue_samples = {}
        self.pending_files = []
        self.config_file = os.path.join(RESOURCES_DIR, "pipeline_config.json")
        
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
        
        # Initialize update flags
        self._update_pending = False
        self._update_after_id = None
        self._last_update = 0
        self.UPDATE_INTERVAL = 500  # ms
        
        # Initialize file monitoring
        self.file_queue = queue.Queue()
        self.event_handler = None
        self.observer = None
        
        # Start processing queue
        self.process_queue()

    def load_config(self):
        """Load analysis configuration from file"""
        # Ensure resources directory exists
        os.makedirs(RESOURCES_DIR, exist_ok=True)
        
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
            print(f"No configuration file found at {self.config_file}")
            return False

    def create_widgets(self):
        print("Creating FolderWatcherWindow widgets")
        
        # Main container
        main_container = ctk.CTkFrame(self)
        main_container.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Create top section
        self.create_top_section(main_container)
        
        # Create separator
        self.create_separator(main_container)
        
        # Create file trees section
        self.create_file_trees(main_container)
        
        # Create bottom controls
        self.create_bottom_controls(main_container)
        
        # Configure cross-platform consistent treeview styling
        self.configure_treeview_colors()
        
        print("FolderWatcherWindow widgets created")

    def create_separator(self, parent):
        """Create a visual separator line"""
        separator = ctk.CTkFrame(parent, height=2)
        separator.pack(fill="x", pady=10)

    def create_top_section(self, container):
        top_frame = ctk.CTkFrame(container, fg_color="transparent")
        top_frame.pack(fill="x", padx=20, pady=(10, 0))
        
        # Folder selection frame
        folder_frame = ctk.CTkFrame(top_frame, fg_color="transparent")
        folder_frame.pack(fill="x", pady=5)
        
        # Folder entry
        self.folder_entry = ctk.CTkEntry(folder_frame, height=32, width=400)
        self.folder_entry.pack(side="left", padx=(0, 5))
        
        # Browse button
        self.browse_button = ctk.CTkButton(
            folder_frame,
            text="Browse",
            width=100,
            height=32,
            command=self.browse_folder
        )
        self.browse_button.pack(side="left", padx=5)
        
        # Import CSV button
        self.import_csv_button = ctk.CTkButton(
            folder_frame,
            text="Import CSV",
            width=100,
            height=32,
            command=self.load_from_csv
        )
        self.import_csv_button.pack(side="left", padx=5)
        
        # Stop button (hidden initially)
        self.stop_button = ctk.CTkButton(
            folder_frame,
            text="Stop Watching",
            width=120,
            height=32,
            command=self.stop_watching
        )
        
        # Auto-analyze checkbox
        self.auto_analyze_checkbox = ctk.CTkCheckBox(
            folder_frame,
            text="Auto-analyze",
            variable=self.auto_analyze_var,
            font=("Helvetica", 11)
        )
        self.auto_analyze_checkbox.pack(side="left", padx=10)

        # Project information section
        # Project, subproject and date in one line
        info_frame = ctk.CTkFrame(top_frame, fg_color="transparent")
        info_frame.pack(fill="x", pady=(10, 0))

        # Initialize date variables
        self.selected_date = datetime.datetime.now().date()
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

        # Pipeline selection
        ctk.CTkLabel(info_frame, text="Pipeline:", 
                     font=("Helvetica", 12, "bold"), 
                     width=60).pack(side="left", padx=(10,0))
        pipeline_options = ["taxonomy-wgs", "taxonomy-16s", "assembly"]
        self.pipeline_menu = ctk.CTkOptionMenu(
            info_frame,
            values=pipeline_options,
            variable=self.analysis_type,
            height=32,
            width=150,
            command=self.on_analysis_type_change
        )
        self.pipeline_menu.pack(side="left", padx=5)

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
        controls_frame = ctk.CTkFrame(container)
        controls_frame.pack(fill="x", pady=(10, 0))
        
        # Left side buttons
        left_buttons = ctk.CTkFrame(controls_frame)
        left_buttons.pack(side="left", fill="x", expand=True)
        
        analyze_button = ctk.CTkButton(
            left_buttons,
            text="Analyze Selected",
            command=self.analyze_selected
        )
        analyze_button.pack(side="left", padx=5)
        
        # Right side buttons
        right_buttons = ctk.CTkFrame(controls_frame)
        right_buttons.pack(side="right", fill="x")
        
        clear_button = ctk.CTkButton(
            right_buttons,
            text="Clear Queue",
            command=self.clear_queue
        )
        clear_button.pack(side="right", padx=5)

    def schedule_update(self):
        """Schedule a tree update with rate limiting"""
        current_time = time.time() * 1000
        if not self._update_pending and (current_time - self._last_update) >= self.UPDATE_INTERVAL:
            self._update_pending = True
            self.after(0, self._do_update)

    def _do_update(self):
        """Perform the actual tree update"""
        try:
            self.update_file_tree()
            self._last_update = time.time() * 1000
        finally:
            self._update_pending = False

    def process_queue(self):
        """Process the file event queue"""
        try:
            while True:
                event = self.file_queue.get_nowait()
                self.handle_file_event(event)
        except queue.Empty:
            pass
        finally:
            # Schedule next queue check
            self.after(100, self.process_queue)

    def handle_file_event(self, event):
        """Handle a single file event"""
        try:
            if event.event_type in ['created', 'modified']:
                if event.src_path.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                    self.add_new_file(event.src_path, is_initial_scan=False)
            elif event.event_type == 'deleted':
                self.remove_file(event.src_path)
            
            # Schedule update if needed
            self.schedule_update()
            
        except Exception as e:
            print(f"Error handling file event: {e}")
            traceback.print_exc()

    def scan_existing_files(self, folder):
        """Scan and add existing files to the file tree"""
        try:
            # Get all files in the folder
            all_files = []
            for root, _, files in os.walk(folder):
                for file in files:
                    if file.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                        all_files.append(os.path.join(root, file))
            
            # Process files in larger batches for better performance
            batch_size = 100
            for i in range(0, len(all_files), batch_size):
                batch = all_files[i:i + batch_size]
                for file_path in batch:
                    # Pass is_initial_scan=True to prevent queueing existing files
                    self.add_new_file(file_path, is_initial_scan=True)
                
                # Update GUI less frequently
                if i % (batch_size * 4) == 0:
                    self.schedule_update()
            
            # Final update
            self.schedule_update()
            
        except Exception as e:
            print(f"Error scanning files: {e}")
            traceback.print_exc()

    def add_new_file(self, file_path, is_initial_scan=False):
        """Add a new file to the samples dict"""
        try:
            if not file_path.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                return
            
            sample_name = self.suggest_sample_name(file_path)
            if not sample_name:
                return
            
            # Add to samples dict
            if sample_name not in self.all_samples:
                self.all_samples[sample_name] = {'files': [], 'status': 'Pending'}
        
            if file_path not in self.all_samples[sample_name]['files']:
                self.all_samples[sample_name]['files'].append(file_path)
            
            # Only add to queue if it's not from initial scan and we're watching
            if not is_initial_scan and self.watching:
                if sample_name not in self.queue_samples:
                    self.queue_samples[sample_name] = {
                        "files": [file_path],
                        "first_seen": datetime.datetime.now().strftime("%H:%M:%S")
                    }
                elif file_path not in self.queue_samples[sample_name]["files"]:
                    self.queue_samples[sample_name]["files"].append(file_path)
                
                # Update queue display
                self.update_queue_display()
                
                # Check if we should auto-analyze
                if self.auto_analyze_var.get():
                    self.check_and_start_analysis(sample_name)
            
        except Exception as e:
            print(f"Error adding file: {e}")

    def update_file_tree(self):
        """Update the file tree with current samples"""
        try:
            # Clear existing items
            for item in self.file_tree.get_children():
                self.file_tree.delete(item)
            
            # Add samples and their files
            for sample_name, sample_data in self.all_samples.items():
                # Add sample
                sample_item = self.file_tree.insert('', 'end', text=sample_name)
                
                # Add status tag if needed
                if 'status' in sample_data:
                    self.file_tree.item(sample_item, tags=(sample_data['status'],))
                
                # Add files
                for file_path in sample_data['files']:
                    self.file_tree.insert(sample_item, 'end', text=os.path.basename(file_path))
            
        except Exception as e:
            print(f"Error updating file tree: {e}")
            traceback.print_exc()

    def toggle_watching(self):
        """Toggle folder watching on/off"""
        if not self.watching:
            folder = self.folder_entry.get()
            if not folder:
                messagebox.showerror("Error", "Please select a folder to watch.")
                return
                
            if not self.check_config():
                return
                
            # Clear existing data
            self.file_tree.delete(*self.file_tree.get_children())
            self.all_samples.clear()
            self.queue_samples.clear()  # Also clear the queue
            
            # Start watching
            self.start_watching(folder)
            
            # Update UI
            self.watching = True
            self.watch_button.configure(
                text="Stop Watching",
                state="normal",  # Ensure button stays enabled
                command=self.toggle_watching  # Ensure command is preserved
            )
            self.folder_entry.configure(state="disabled")
            self.browse_button.configure(state="disabled")
            
            messagebox.showinfo("Folder Watch Started", f"Now watching folder: {folder}")
        else:
            self.stop_watching()
            
            # Update UI
            self.watching = False
            self.watch_button.configure(
                text="Watch Folder",
                state="normal",  # Ensure button stays enabled
                command=self.toggle_watching  # Ensure command is preserved
            )
            self.folder_entry.configure(state="normal")
            self.browse_button.configure(state="normal")

    def start_watching(self, folder):
        """Start watching a folder for changes"""
        try:
            # Stop existing observer if any
            if self.observer:
                self.stop_watching()
            
            # Scan existing files first
            self.scan_existing_files(folder)
            
            # Set up new observer
            self.event_handler = FileSystemEventHandler()
            self.event_handler.on_created = lambda event: self.file_queue.put(event)
            self.event_handler.on_modified = lambda event: self.file_queue.put(event)
            self.event_handler.on_deleted = lambda event: self.file_queue.put(event)
            
            self.observer = Observer()
            self.observer.schedule(self.event_handler, folder, recursive=True)
            self.observer.start()
            
            # Record start time
            self.watch_start_time = datetime.datetime.now()
            
        except Exception as e:
            print(f"Error starting folder watch: {e}")
            traceback.print_exc()
            messagebox.showerror("Error", f"Failed to start watching folder: {str(e)}")

    def stop_watching(self):
        """Stop watching and clean up resources"""
        try:
            if self.observer:
                self.observer.stop()
                self.observer.join()
                self.observer = None
            
            self.event_handler = None
            self.watch_start_time = None
            
            # Clear queue samples but keep all_samples
            self.queue_samples.clear()
            self.update_queue_display()
            
        except Exception as e:
            print(f"Error stopping folder watch: {e}")
            traceback.print_exc()

    def browse_folder(self):
        """Browse for a folder to watch"""
        folder = filedialog.askdirectory()
        if folder:
            self.folder_entry.delete(0, tk.END)
            self.folder_entry.insert(0, folder)
            
            # Suggest project names
            self.suggest_project_names(folder)
            
            # Start watching automatically
            self.start_watching_folder()
            
    def start_watching_folder(self):
        """Start watching the selected folder"""
        folder = self.folder_entry.get()
        if not folder:
            messagebox.showerror("Error", "Please select a folder to watch.")
            return
            
        if not self.check_config():
            return
            
        # Clear existing data
        self.file_tree.delete(*self.file_tree.get_children())
        self.all_samples.clear()
        self.queue_samples.clear()
        
        # Start watching
        self.start_watching(folder)
        
        # Update UI
        self.watching = True
        self.folder_entry.configure(state="disabled")
        self.browse_button.pack_forget()  # Hide browse button
        self.stop_button.pack(side="left", padx=5)  # Show stop button
        
        messagebox.showinfo("Folder Watch Started", f"Now watching folder: {folder}")

    def stop_watching(self):
        """Stop watching and clean up resources"""
        try:
            if self.observer:
                self.observer.stop()
                self.observer.join()
                self.observer = None
            
            self.event_handler = None
            self.watch_start_time = None
            
            # Clear queue samples but keep all_samples
            self.queue_samples.clear()
            self.update_queue_display()
            
            # Update UI
            self.watching = False
            self.folder_entry.configure(state="normal")
            self.stop_button.pack_forget()  # Hide stop button
            self.browse_button.pack(side="left", padx=5)  # Show browse button
            
        except Exception as e:
            print(f"Error stopping folder watch: {e}")
            traceback.print_exc()

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
                      year=datetime.datetime.now().year,
                      month=datetime.datetime.now().month,
                      day=datetime.datetime.now().day)
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
        samples_to_analyze = []  # list to track samples that need to be analyzed

        try:
            with open(file_path, 'r') as file:
                reader = csv.DictReader(file)
                
                # Verify required columns
                required_columns = ["sample_name", "date", "project_name", "subproject_name", "sample folder"]
                missing_columns = [col for col in required_columns if col not in reader.fieldnames]
                if missing_columns:
                    messagebox.showerror("Error", f"Missing required columns in CSV: {', '.join(missing_columns)}")
                    return

                # Get current analysis type
                current_analysis_type = self.analysis_type.get()

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
                        "path": folder_path,
                        "project": row["project_name"],
                        "subproject": row["subproject_name"],
                        "date": row["date"]
                    }
                    
                    # Update project and subproject entries with the first row's values
                    if not self.project_entry.get():
                        self.project_entry.insert(0, row["project_name"])
                    if not self.subproject_entry.get():
                        self.subproject_entry.insert(0, row["subproject_name"])
                    
                    # Add to queue and prepare for analysis
                    current_time = datetime.datetime.now().strftime("%H:%M:%S")
                    self.queue_samples[sample_name] = {
                        "files": files,
                        "first_seen": current_time,
                        "project": row["project_name"],
                        "subproject": row["subproject_name"],
                        "date": row["date"]
                    }
                    samples_to_analyze.append(sample_name)

            # Update displays
            self.schedule_update()
            
            # Show summary
            num_samples = len(samples_to_analyze)
            if num_samples > 0:
                # Start analysis for all samples
                progress_window = tk.Toplevel(self)
                progress_window.title("Analysis Progress")
                
                # Create progress bar
                progress_label = ttk.Label(progress_window, text="Starting analysis...", wraplength=300)
                progress_label.pack(pady=10)
                progress_bar = ttk.Progressbar(progress_window, mode='indeterminate')
                progress_bar.pack(pady=10, padx=20, fill=tk.X)
                progress_bar.start()
                
                # Start analysis for each sample
                for sample_name in samples_to_analyze:
                    sample_data = self.queue_samples[sample_name]
                    progress_label.config(text=f"Analyzing {sample_name}...")
                    
                    # Update project, subproject, and date from CSV
                    self.project_entry.delete(0, tk.END)
                    self.project_entry.insert(0, sample_data["project"])
                    self.subproject_entry.delete(0, tk.END)
                    self.subproject_entry.insert(0, sample_data["subproject"])
                    self.date_entry.delete(0, tk.END)
                    self.date_entry.insert(0, sample_data["date"])
                    
                    # Run analysis
                    self.run_analysis_for_sample(sample_name)
                
                progress_window.destroy()
                messagebox.showinfo(
                    "Import Complete",
                    f"Successfully imported and started analysis for {num_samples} sample{'s' if num_samples != 1 else ''}."
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
            traceback.print_exc()

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
            file_count = sum(len(data["files"]) for data in self.queue_samples.values())
            self.queue_count.configure(text=f"({queue_count} samples, {file_count} files)")
            
        except Exception as e:
            print(f"Error updating queue tree: {e}")
            traceback.print_exc()

    def create_control_buttons(self, container):
        """Create control buttons with all required functionality"""
        button_frame = ctk.CTkFrame(container, fg_color="transparent")
        button_frame.pack(fill="x", pady=5)
        
        # Left side controls
        left_controls = ctk.CTkFrame(button_frame, fg_color="transparent")
        left_controls.pack(side="left", fill="x", expand=True)
        
        # Analyze Selected button
        self.analyze_button = ctk.CTkButton(
            left_controls,
            text="Analyze Selected",
            command=self.analyze_selected
        )
        self.analyze_button.pack(side="left", padx=5)
        create_tooltip(self.analyze_button, "Analyze selected samples (Shift+Click for multiple)")
        
        # Run All button
        self.run_all_button = ctk.CTkButton(
            left_controls,
            text="Run All",
            command=self.run_all_analyses
        )
        self.run_all_button.pack(side="left", padx=5)
        create_tooltip(self.run_all_button, "Run analysis for all samples in the file tree")
        
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
        create_tooltip(self.auto_analyze_checkbox, "Automatically analyze new samples when they are added")
        
        # Run Queue button
        self.run_queue_button = ctk.CTkButton(
            right_controls,
            text="Run Queue",
            command=self.analyze_selected_queue
        )
        self.run_queue_button.pack(side="left", padx=5)
        create_tooltip(self.run_queue_button, "Run analysis for selected samples in the queue")

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

    class ProgressPopup(tk.Toplevel):
        def __init__(self, parent, title="Progress", cancel_callback=None):
            super().__init__(parent)
            self.title(title)
            self.parent = parent
            self.cancel_callback = cancel_callback
            
            # Set up the window
            self.geometry("300x150")
            self.resizable(False, False)
            
            # Progress label
            self.label = ttk.Label(self, text="Processing...", wraplength=280)
            self.label.pack(pady=10)
            
            # Progress bar
            self.progress = ttk.Progressbar(self, mode='indeterminate')
            self.progress.pack(pady=10, padx=20, fill=tk.X)
            self.progress.start()
            
            # Cancel button
            self.cancel_btn = ttk.Button(self, text="Cancel", command=self.cancel)
            self.cancel_btn.pack(pady=10)
            
            # Handle window close button
            self.protocol("WM_DELETE_WINDOW", self.cancel)
            
            # Center the window
            self.center_window()
            
        def cancel(self):
            if self.cancel_callback:
                self.cancel_callback()
            self.destroy()
            
        def center_window(self):
            self.update_idletasks()
            width = self.winfo_width()
            height = self.winfo_height()
            x = (self.winfo_screenwidth() // 2) - (width // 2)
            y = (self.winfo_screenheight() // 2) - (height // 2)
            self.geometry('{}x{}+{}+{}'.format(width, height, x, y))
            
        def update_progress(self, text):
            self.label.config(text=text)
            
    def run_analysis_for_sample(self, sample_name):
        """Run analysis for a single sample"""
        if sample_name in self.running_analyses:
            print(f"Analysis already running for {sample_name}")
            return
            
        # Create progress popup
        self.progress_popup = self.ProgressPopup(
            self, 
            title=f"Running {self.analysis_type.get()} for {sample_name}",
            cancel_callback=lambda: self.cancel_analysis(sample_name)
        )
        
        # Create and start analysis thread
        analysis_thread = threading.Thread(
            target=self._run_analysis_thread,
            args=(sample_name,)
        )
        analysis_thread.daemon = True
        
        # Store thread reference
        self.running_analyses[sample_name] = {
            'thread': analysis_thread,
            'cancel': False
        }
        
        analysis_thread.start()
        
    def cancel_analysis(self, sample_name):
        """Cancel running analysis for a sample"""
        if sample_name in self.running_analyses:
            print(f"Canceling analysis for {sample_name}")
            self.running_analyses[sample_name]['cancel'] = True
            # Remove from running analyses
            del self.running_analyses[sample_name]
            # Close progress popup if it exists
            if hasattr(self, 'progress_popup') and self.progress_popup:
                self.progress_popup.destroy()
                self.progress_popup = None
                
    def _run_analysis_thread(self, sample_name):
        """Thread function to run analysis for a single sample"""
        try:
            # Get files for this sample and filter out any concatenated files
            files = [f for f in self.all_samples[sample_name]['files'] 
                    if not any(x in f for x in ['_concatenated.fastq', '_concat.fastq'])]
            if not files:
                raise Exception("No files found for sample")
            
            # Initialize MMonitorCMD
            cmd_runner = MMonitorCMD()
            cmd_runner.offline_mode = self.offline_mode
            cmd_runner.logged_in = self.logged_in
            cmd_runner.current_user = self.current_user
            
            # Load config file for threads
            try:
                with open(self.config_file, 'r') as f:
                    config = json.load(f)
            except Exception as e:
                print(f"Error loading config: {e}")
                self.after(0, lambda: self.analysis_failed(sample_name))
                return
                
            # Set up args for analysis
            args = argparse.Namespace(
                analysis=self.analysis_type.get(),
                config=self.config_file,
                threads=int(config.get('threads', 12)),
                sample=sample_name,
                project=self.project_var.get(),
                subproject=self.subproject_var.get(),
                date=datetime.datetime.strptime(self.date_var.get(), '%Y-%m-%d').date(),
                input=files,
                multicsv=None,
                centrifuger_db=config.get('centrifuge_db', ''),
                emu_db=os.path.abspath(config.get('emu_db', '')),
                min_length=config.get('min_length', 1000),
                min_quality=config.get('min_quality', 10.0),
                min_abundance=float(config.get('min_abundance', 0.01)),
                overwrite=True,
                qc=True,
                update=False,
                verbose=True,
                loglevel='INFO'
            )
            
            # Initialize and run analysis
            cmd_runner.initialize_from_args(args)
            
            # Check for cancellation before starting
            if sample_name in self.running_analyses and self.running_analyses[sample_name]['cancel']:
                print(f"Analysis canceled for {sample_name}")
                return
                
            success = cmd_runner.run()
            
            # Check for cancellation after run
            if sample_name in self.running_analyses and self.running_analyses[sample_name]['cancel']:
                print(f"Analysis canceled for {sample_name}")
                return
            
            if success:
                # Get the output directory
                output_dir = os.path.join(cmd_runner.pipeline_out, sample_name)
                
                # Update database with results
                if self.analysis_type.get() == "taxonomy-16s":
                    print(f"\nUpdating database with EMU results for {sample_name}...")
                    try:
                        cmd_runner.django_db.update_django_with_emu_out(
                            emu_out_path=output_dir,
                            tax_rank="species",
                            sample_name=sample_name,
                            project_name=self.project_var.get(),
                            subproject_name=self.subproject_var.get(),
                            sample_date=self.date_var.get(),
                            overwrite=True
                        )
                        print("Database update completed")
                    except Exception as e:
                        print(f"Error updating database: {e}")
                        self.after(0, lambda: self.analysis_failed(sample_name))
                        return
                        
                elif self.analysis_type.get() == "taxonomy-wgs":
                    print(f"\nUpdating database with Centrifuger results for {sample_name}...")
                    try:
                        # Construct the path to the Centrifuger report file
                        centrifuger_report = os.path.join(output_dir, f"{sample_name}_centrifuger_report.tsv")
                        if not os.path.isfile(centrifuger_report):
                            raise FileNotFoundError(f"Centrifuger report not found at {centrifuger_report}")
                            
                        cmd_runner.django_db.send_nanopore_record_centrifuger(
                            kraken_out_path=centrifuger_report,
                            sample_name=sample_name,
                            project_id=self.project_var.get(),
                            subproject_id=self.subproject_var.get(),
                            date=datetime.datetime.now().strftime("%Y-%m-%d"),
                            overwrite=True
                        )
                        print("Successfully updated database with Centrifuger results")
                    except Exception as e:
                        print(f"Error updating database: {str(e)}")
                        self.after(0, lambda: self.analysis_failed(sample_name))
                        return
                
                self.after(0, lambda: self.analysis_completed(sample_name))
            else:
                self.after(0, lambda: self.analysis_failed(sample_name))
            
        except Exception as e:
            print(f"Error running analysis for {sample_name}: {e}")
            self.after(0, lambda: self.analysis_failed(sample_name))
        finally:
            # Clean up
            if sample_name in self.running_analyses:
                del self.running_analyses[sample_name]
            
            # Close progress window in main thread
            self.after(0, self.close_progress_popup)
            
    def analysis_completed(self, sample_name):
        """Handle analysis completion in main thread"""
        try:
            if hasattr(self, 'progress_popup') and self.progress_popup.winfo_exists():
                self.progress_popup.update_progress(f"Analysis completed for {sample_name}")
        except Exception as e:
            print(f"Error updating progress window: {e}")
        
        self.update_sample_status(sample_name, "Completed")

    def analysis_failed(self, sample_name):
        """Handle analysis failure in main thread"""
        try:
            if hasattr(self, 'progress_popup') and self.progress_popup.winfo_exists():
                self.progress_popup.update_progress(f"Analysis failed for {sample_name}")
        except Exception as e:
            print(f"Error updating progress window: {e}")
        
        self.update_sample_status(sample_name, "Failed")

    def close_progress_popup(self):
        """Safely close the progress popup"""
        try:
            if hasattr(self, 'progress_popup') and self.progress_popup.winfo_exists():
                self.progress_popup.destroy()
                del self.progress_popup
        except Exception as e:
            print(f"Error closing progress window: {e}")

    def remove_file(self, file_path):
        """Remove a file from the samples dict"""
        try:
            for sample_name, sample_data in self.all_samples.items():
                if file_path in sample_data['files']:
                    self.all_samples[sample_name]['files'].remove(file_path)
                    if not self.all_samples[sample_name]['files']:
                        del self.all_samples[sample_name]
                    break
            
        except Exception as e:
            print(f"Error removing file: {e}")

    def add_to_queue(self, sample_name):
        """Add a sample to the queue"""
        try:
            if sample_name not in self.queue_samples:
                self.queue_samples[sample_name] = {
                    "files": self.all_samples[sample_name]['files'],
                    "first_seen": datetime.datetime.now().strftime("%H:%M:%S")
                }
            
        except Exception as e:
            print(f"Error adding to queue: {e}")

    def on_file_change(self, event):
        """Handle file system events"""
        try:
            if event.is_directory:
                return
                
            file_path = event.src_path
            if not file_path.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
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
            
            self.schedule_update()
            
        except Exception as e:
            print(f"Error handling file change: {e}")

    def on_analysis_type_change(self, *args):
        """Handle analysis type change"""
        analysis_type = self.analysis_type.get()
        print(f"Analysis type changed to: {analysis_type}")
        # Update configuration based on analysis type
        if hasattr(self, 'gui') and hasattr(self.gui, 'pipeline_config'):
            self.gui.pipeline_config['analysis_type'] = analysis_type

    def destroy(self):
        """Clean up before destroying the window"""
        try:
            # Remove the analysis type trace
            if hasattr(self, 'analysis_type') and hasattr(self, 'analysis_type_trace'):
                self.analysis_type.trace_remove('write', self.analysis_type_trace)
            
            # Stop watching if active
            if self.watching:
                self.stop_watching()
            
            # Destroy all widgets
            super().destroy()
            
        except Exception as e:
            print(f"Error during cleanup: {e}")
            traceback.print_exc()
            # Still try to destroy even if cleanup failed
            super().destroy()

    def analyze_selected(self):
        """Analyze selected samples from the file tree"""
        selected_items = self.file_tree.selection()
        if not selected_items:
            messagebox.showwarning("No Selection", "Please select at least one sample to analyze")
            return
        
        for item in selected_items:
            # Get the parent item if a file is selected
            parent = self.file_tree.parent(item)
            sample_item = parent if parent else item
            sample_name = self.file_tree.item(sample_item)['text']
            
            if sample_name in self.all_samples:
                self.run_analysis_for_sample(sample_name)

    def update_sample_status(self, sample_name, status):
        """Update the status of a sample in the file tree"""
        try:
            # Find the sample item
            for item in self.file_tree.get_children():
                if self.file_tree.item(item)['text'] == sample_name:
                    # Update status
                    self.file_tree.item(item, tags=(status,))
                    break
            
        except Exception as e:
            print(f"Error updating sample status: {e}")