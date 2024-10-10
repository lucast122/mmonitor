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
from mmonitor.userside.PipelineWindow import PipelinePopup
from Bio import SeqIO
import numpy as np


class FolderWatcherWindow(ctk.CTkFrame):
    def __init__(self, parent, gui_ref):
        print("Initializing FolderWatcherWindow")
        super().__init__(parent)
        self.parent = parent
        self.gui = gui_ref
        self.samples = {}
        self.watching = False
        self.selected_date = datetime.date.today()
        self.auto_analyze_var = tk.BooleanVar(value=False)
        self.use_suggested_name_var = tk.BooleanVar(value=True)
        self.use_file_date_var = tk.BooleanVar(value=True)
        self.config_file = os.path.join(os.path.expanduser("~"), ".mmonitor_config.json")
        self.load_config()
        self.analysis_type = tk.StringVar(value="taxonomy-wgs")
        self.create_widgets()
        
        self.update_queue = queue.Queue()
        self.parent.after(100, self.process_queue)
        print("FolderWatcherWindow initialized")

        self.update_appearance()

    def load_config(self):
        if os.path.exists(self.config_file):
            with open(self.config_file, 'r') as f:
                self.config = json.load(f)
        else:
            messagebox.showwarning("Configuration Missing", "Please set up analysis parameters in the Analysis tab first.")
            self.config = {}

    def create_widgets(self):
        print("Creating FolderWatcherWindow widgets")
        
        main_frame = ctk.CTkFrame(self)
        main_frame.pack(fill="both", expand=True, padx=20, pady=20)

        title_label = ctk.CTkLabel(main_frame, text="Folder Watcher", font=("Helvetica", 24, "bold"))
        title_label.pack(pady=(0, 10))
        
        explanation = ("Monitor a folder for new sequencing files and automatically analyze samples.\n"
                       "Select multiple samples for batch analysis when 'Use suggested name' is checked.")
        ctk.CTkLabel(main_frame, text=explanation, wraplength=500).pack(pady=(0, 20))

        folder_frame = ctk.CTkFrame(main_frame)
        folder_frame.pack(fill="x", pady=(0, 10))

        ctk.CTkLabel(folder_frame, text="Folder to watch:").pack(side="left", padx=5)
        self.folder_entry = ctk.CTkEntry(folder_frame, width=300)
        self.folder_entry.pack(side="left", padx=5, expand=True, fill="x")
        browse_button = ctk.CTkButton(folder_frame, text="Browse", command=self.browse_folder)
        browse_button.pack(side="left", padx=5)
        create_tooltip(browse_button, "Select a folder to monitor for new sequencing files.")
        
        self.watch_button = ctk.CTkButton(folder_frame, text="Start Watching", command=self.toggle_watching)
        self.watch_button.pack(side="left", padx=5)
        create_tooltip(self.watch_button, "Start or stop monitoring the selected folder.")

        info_frame = ctk.CTkFrame(main_frame)
        info_frame.pack(fill="x", pady=10)

        sample_name_frame = ctk.CTkFrame(info_frame)
        sample_name_frame.pack(fill="x", pady=5)
        ctk.CTkLabel(sample_name_frame, text="Sample Name:", width=120, anchor="e").pack(side="left", padx=5)
        self.sample_name_entry = ctk.CTkEntry(sample_name_frame, width=200)
        self.sample_name_entry.pack(side="left", padx=5)
        self.use_suggested_checkbox = ctk.CTkCheckBox(sample_name_frame, text="Use suggested name", variable=self.use_suggested_name_var, 
                        command=self.toggle_sample_name_entry)
        self.use_suggested_checkbox.pack(side="left", padx=5)
        create_tooltip(self.use_suggested_checkbox, "Automatically use suggested sample names based on file paths.\nEnables multi-sample selection.")

        for label, attr in [("Project Name:", "project_entry"), ("Subproject Name:", "subproject_entry")]:
            frame = ctk.CTkFrame(info_frame)
            frame.pack(fill="x", pady=5)
            ctk.CTkLabel(frame, text=label, width=120, anchor="e").pack(side="left", padx=5)
            entry = ctk.CTkEntry(frame, width=200)
            entry.pack(side="left", padx=5)
            setattr(self, attr, entry)

        date_frame = ctk.CTkFrame(info_frame)
        date_frame.pack(fill="x", pady=5)
        ctk.CTkLabel(date_frame, text="Date:", width=120, anchor="e").pack(side="left", padx=5)
        self.date_btn = ctk.CTkButton(date_frame, text="Select Date", command=self.open_calendar)
        self.date_btn.pack(side="left", padx=5)
        use_file_date_checkbox = ctk.CTkCheckBox(date_frame, text="Use file creation date", variable=self.use_file_date_var)
        use_file_date_checkbox.pack(side="left", padx=5)
        create_tooltip(use_file_date_checkbox, "Use the file creation date instead of manually selecting a date.")

        analysis_frame = ctk.CTkFrame(main_frame)
        analysis_frame.pack(fill="x", pady=10)
        ctk.CTkLabel(analysis_frame, text="Analysis Type:").pack(side="left", padx=5)
        analysis_menu = ctk.CTkOptionMenu(analysis_frame, variable=self.analysis_type, 
                                          values=["taxonomy-wgs", "taxonomy-16s", "assembly-functional"])
        analysis_menu.pack(side="left", padx=5)
        create_tooltip(analysis_menu, "Select the type of analysis to perform")

        self.create_file_treeview(main_frame)

        control_frame = ctk.CTkFrame(main_frame)
        control_frame.pack(fill="x", pady=10)
        
        self.run_selected_button = ctk.CTkButton(control_frame, text="Run Selected Sample", command=self.run_selected_sample)
        self.run_selected_button.pack(side="left", padx=5)
        create_tooltip(self.run_selected_button, "Run analysis for the selected sample")

        self.run_all_button = ctk.CTkButton(control_frame, text="Run All Analyses", command=self.run_all_analyses)
        self.run_all_button.pack(side="left", padx=5)
        create_tooltip(self.run_all_button, "Run analysis for all samples")

        auto_analyze_checkbox = ctk.CTkCheckBox(control_frame, text="Auto-analyze new samples", variable=self.auto_analyze_var)
        auto_analyze_checkbox.pack(side="left", padx=5)
        create_tooltip(auto_analyze_checkbox, "Automatically start analysis when new samples are detected")

        self.configure_treeview_colors()
        print("FolderWatcherWindow widgets created")

    def create_file_treeview(self, parent):
        tree_frame = ttk.Frame(parent)
        tree_frame.pack(fill="both", expand=True, pady=10)

        self.file_tree = ttk.Treeview(tree_frame, columns=("files", "status", "read_type"), selectmode="extended")
        self.file_tree.heading("#0", text="Sample")
        self.file_tree.heading("files", text="Files")
        self.file_tree.heading("status", text="Status")
        self.file_tree.heading("read_type", text="Detected Read Type")
        self.file_tree.pack(side="left", fill="both", expand=True)

        scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.file_tree.yview)
        scrollbar.pack(side="right", fill="y")
        self.file_tree.configure(yscrollcommand=scrollbar.set)

    def configure_treeview_colors(self):
        style = ttk.Style()
        if ctk.get_appearance_mode() == "Dark":
            style.configure("Treeview", 
                            background="gray20", 
                            foreground="white", 
                            fieldbackground="gray20")
            style.map('Treeview', background=[('selected', 'gray30')])
        else:
            style.configure("Treeview", 
                            background="white", 
                            foreground="black", 
                            fieldbackground="white")
            style.map('Treeview', background=[('selected', 'gray70')])

    def browse_folder(self):
        folder = filedialog.askdirectory()
        if folder:
            self.folder_entry.delete(0, tk.END)
            self.folder_entry.insert(0, folder)
            self.suggest_project_names(folder)
            self.scan_existing_files(folder)

    def suggest_project_names(self, folder):
        project_name = os.path.basename(folder)
        subproject_name = os.path.basename(os.path.dirname(folder))
        
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
        if os.path.exists(self.config_file):
            with open(self.config_file, 'r') as f:
                pipeline_config = json.load(f)
        else:
            messagebox.showerror("Error", "Configuration file not found. Please set up analysis parameters in the Analysis tab first.")
            return

        pipeline_popup = PipelinePopup(self, self.gui)
        missing_params = pipeline_popup.check_missing_parameters(pipeline_config)
        if missing_params:
            messagebox.showerror("Error", f"Cannot start watching. The following parameters are missing in the pipeline configuration:\n{', '.join(missing_params)}")
            return

        folder = self.folder_entry.get()
        if not folder:
            messagebox.showerror("Error", "Please select a folder to watch.")
            return

        self.samples.clear()
        self.file_tree.delete(*self.file_tree.get_children())

        self.scan_existing_files(folder)

        if hasattr(self, 'folder_monitor'):
            self.folder_monitor.stop()

        self.folder_monitor = FolderMonitor([folder], self)
        self.folder_monitor.start()
        self.watching = True
        self.watch_button.configure(text="Stop Watching")

        messagebox.showinfo("Folder Watch Started", f"Now watching folder: {folder}")

    def stop_watching(self):
        if hasattr(self, 'folder_monitor'):
            self.folder_monitor.stop()
        self.watching = False
        self.watch_button.configure(text="Start Watching")

    def scan_existing_files(self, folder):
        for root, dirs, files in os.walk(folder):
            fastq_files = [f for f in files if f.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz'))]
            if fastq_files:
                sample_name = self.suggest_sample_name(root)
                file_paths = [os.path.join(root, f) for f in fastq_files]
                self.add_sample(sample_name, file_paths)
    """
    Detect if the reads look like 16S or WGS based on read length distribution.
    """
    
    def detect_read_type(self, file_path, length_threshold=2000, fraction_threshold=0.6, max_reads=100):
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

    def run_selected_sample(self):
        selected_items = self.file_tree.selection()
        if not selected_items:
            messagebox.showinfo("Info", "Please select a sample to analyze.")
            return
        
        for item in selected_items:
            if self.file_tree.parent(item) == "":  # Only process parent items (samples)
                self.run_analysis_for_sample(item)

    def run_all_analyses(self):
        for sample in self.samples:
            self.run_analysis_for_sample(sample)

    def run_analysis_for_sample(self, sample_name):
        analysis_type = self.analysis_type.get()
        
        with open(self.config_file, 'r') as f:
            config = json.load(f)
        
        # Add conda environment to PATH
        conda_env_path = os.path.expanduser("~/miniconda3/envs/fastcat/bin")
        os.environ["PATH"] = f"{conda_env_path}:{os.environ['PATH']}"
        
        cmd_runner = MMonitorCMD()
        args = [
            "-a", analysis_type,
            "-c", self.gui.db_path,
            "-i"] + self.samples[sample_name] + [
            "-s", sample_name,
            "-d", self.selected_date.strftime("%Y-%m-%d"),
            "-p", self.project_entry.get(),
            "-u", self.subproject_entry.get(),
            "-n", config.get('min_abundance', '0.01'),
            "-t", config.get('threads', '1'),
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

    def process_queue(self):
        try:
            while True:
                method, args = self.update_queue.get_nowait()
                if method == 'add_sample':
                    self.add_sample(*args)
        except queue.Empty:
            pass
        self.after(100, self.process_queue)

    def open_calendar(self):
        def on_close():
            selected_date = cal.selection_get()
            self.date_btn.configure(text=selected_date.strftime("%Y-%m-%d"))
            self.selected_date = selected_date
            date_win.destroy()

        date_win = tk.Toplevel(self)
        date_win.title("Select a Date")

        cal = Calendar(date_win, selectmode='day')
        cal.pack(pady=20, padx=20)

        ctk.CTkButton(date_win, text="OK", command=on_close).pack(pady=20)

    def toggle_sample_name_entry(self):
        if self.use_suggested_name_var.get():
            self.sample_name_entry.configure(state="disabled")
        else:
            self.sample_name_entry.configure(state="normal")

    def update_appearance(self):
        self.configure_treeview_colors()
        
        # Update colors based on the current theme
        bg_color = self.cget("fg_color")
        text_color = "white" if ctk.get_appearance_mode() == "Dark" else "black"

        for widget in self.winfo_children():
            if isinstance(widget, ctk.CTkFrame):
                widget.configure(fg_color=bg_color)
            elif isinstance(widget, (ctk.CTkLabel, ctk.CTkButton, ctk.CTkCheckBox, ctk.CTkRadioButton)):
                widget.configure(text_color=text_color)

        # Update treeview colors
        self.configure_treeview_colors()