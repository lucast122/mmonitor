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
import queue  # Add this import
import threading  # Make sure this is imported as well
from build_mmonitor_pyinstaller import ROOT
from tkcalendar import Calendar

class FolderWatcherWindow(ctk.CTkFrame):
    def __init__(self, parent):
        print("Initializing FolderWatcherWindow")
        super().__init__(parent)
        self.parent = parent
        self.samples = {}
        self.watching = False
        self.selected_date = datetime.date.today()
        self.auto_analyze_var = tk.BooleanVar(value=False)
        self.use_suggested_name_var = tk.BooleanVar(value=False)
        self.create_widgets()
        
        self.queue = queue.Queue()
        self.after(100, self.process_queue)
        print("FolderWatcherWindow initialized")

    def create_widgets(self):
        print("Creating FolderWatcherWindow widgets")
        
        # Main frame
        main_frame = ctk.CTkFrame(self)
        main_frame.pack(fill="both", expand=True, padx=20, pady=20)

        # Title and explanation
        title_label = ctk.CTkLabel(main_frame, text="Folder Watcher", font=("Helvetica", 24, "bold"))
        title_label.pack(pady=(0, 10))
        
        explanation = ("This tool allows you to monitor a folder for new sequencing files.\n"
                       "You can set up automatic analysis for new samples as they arrive.")
        ctk.CTkLabel(main_frame, text=explanation, wraplength=500).pack(pady=(0, 20))

        # Folder selection frame
        folder_frame = ctk.CTkFrame(main_frame)
        folder_frame.pack(fill="x", pady=(0, 10))

        ctk.CTkLabel(folder_frame, text="Folder to watch:").pack(side="left", padx=5)
        self.folder_entry = ctk.CTkEntry(folder_frame, width=300)
        self.folder_entry.pack(side="left", padx=5, expand=True, fill="x")
        ctk.CTkButton(folder_frame, text="Browse", command=self.browse_folder).pack(side="left", padx=5)
        self.watch_button = ctk.CTkButton(folder_frame, text="Start Watching", command=self.toggle_watching)
        self.watch_button.pack(side="left", padx=5)

        # Sample information frame
        info_frame = ctk.CTkFrame(main_frame)
        info_frame.pack(fill="x", pady=10)

        # Sample name
        sample_name_frame = ctk.CTkFrame(info_frame)
        sample_name_frame.pack(fill="x", pady=5)
        ctk.CTkLabel(sample_name_frame, text="Sample Name:").pack(side="left", padx=5)
        self.sample_name_entry = ctk.CTkEntry(sample_name_frame, width=200)
        self.sample_name_entry.pack(side="left", padx=5)
        ctk.CTkCheckBox(sample_name_frame, text="Use suggested name", variable=self.use_suggested_name_var, 
                        command=self.toggle_sample_name_entry).pack(side="left", padx=5)

        # Project and subproject
        for label, attr in [("Project Name:", "project_entry"), ("Subproject Name:", "subproject_entry")]:
            frame = ctk.CTkFrame(info_frame)
            frame.pack(fill="x", pady=5)
            ctk.CTkLabel(frame, text=label).pack(side="left", padx=5)
            setattr(self, attr, ctk.CTkEntry(frame, width=200))
            getattr(self, attr).pack(side="left", padx=5)

        # Date selection
        date_frame = ctk.CTkFrame(info_frame)
        date_frame.pack(fill="x", pady=5)
        ctk.CTkLabel(date_frame, text="Date:").pack(side="left", padx=5)
        self.date_btn = ctk.CTkButton(date_frame, text="Select Date", command=self.open_calendar)
        self.date_btn.pack(side="left", padx=5)
        ctk.CTkCheckBox(date_frame, text="Use file creation date").pack(side="left", padx=5)

        # Analysis type selection
        analysis_frame = ctk.CTkFrame(main_frame)
        analysis_frame.pack(fill="x", pady=10)
        self.analysis_var = tk.StringVar(value="taxonomy-wgs")
        ctk.CTkRadioButton(analysis_frame, text="WGS Analysis", variable=self.analysis_var, value="taxonomy-wgs").pack(side="left", padx=10)
        ctk.CTkRadioButton(analysis_frame, text="16S Analysis", variable=self.analysis_var, value="taxonomy-16s").pack(side="left", padx=10)

        # Treeview
        tree_frame = ttk.Frame(main_frame)
        tree_frame.pack(fill="both", expand=True, pady=10)
        self.tree = ttk.Treeview(tree_frame, columns=("status",))
        self.tree.heading("#0", text="Sample/File")
        self.tree.heading("status", text="Status")
        self.tree.pack(side="left", fill="both", expand=True)
        scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        scrollbar.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=scrollbar.set)

        # Control buttons
        control_frame = ctk.CTkFrame(main_frame)
        control_frame.pack(fill="x", pady=10)
        self.start_analysis_button = ctk.CTkButton(control_frame, text="Start Analysis", command=self.start_analysis)
        self.start_analysis_button.pack(side="left", padx=5)
        ctk.CTkCheckBox(control_frame, text="Auto-analyze new samples", variable=self.auto_analyze_var).pack(side="left", padx=5)

        print("FolderWatcherWindow widgets created")

    def browse_folder(self):
        folder = filedialog.askdirectory()
        if folder:
            self.folder_entry.delete(0, tk.END)
            self.folder_entry.insert(0, folder)

    def toggle_watching(self):
        if not self.watching:
            self.start_watching()
        else:
            self.stop_watching()

    def start_watching(self):
        folder = self.folder_entry.get()
        if not folder:
            messagebox.showerror("Error", "Please select a folder to watch.")
            return

        self.samples.clear()
        self.tree.delete(*self.tree.get_children())

        threading.Thread(target=self.scan_existing_files, args=(folder,), daemon=True).start()

        self.folder_monitor = FolderMonitor([folder], self)  # Changed from self.gui.folder_monitor
        self.folder_monitor.start()
        self.watching = True
        self.watch_button.configure(text="Stop Watching")

    def scan_existing_files(self, folder):
        for root, dirs, files in os.walk(folder):
            if "fastq_pass" in dirs:
                fastq_pass_dir = os.path.join(root, "fastq_pass")
                barcode_dirs = [d for d in os.listdir(fastq_pass_dir) if d.startswith("barcode")]
                
                if barcode_dirs:
                    for barcode_dir in barcode_dirs:
                        parent_folder = os.path.basename(os.path.dirname(fastq_pass_dir))
                        self.process_sample_folder(os.path.join(fastq_pass_dir, barcode_dir), f"{parent_folder}_{barcode_dir}")
                else:
                    parent_folder = os.path.basename(os.path.dirname(fastq_pass_dir))
                    self.process_sample_folder(fastq_pass_dir, parent_folder)

    def process_sample_folder(self, folder, sample_name):
        files = [f for f in os.listdir(folder) if f.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')) and not f.endswith('_concatenated.fastq.gz')]
        if files:
            self.queue.put(('add_sample', sample_name, folder, files))

    def add_new_file(self, file_path):
        if "_concatenated" in file_path:
            return  # Ignore concatenated files

        parent_folder = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file_path))))
        barcode_folder = os.path.basename(os.path.dirname(file_path))
        suggested_sample_name = f"{parent_folder}_{barcode_folder}"
        self.queue.put(('add_file', suggested_sample_name, file_path))

    def process_queue(self):
        try:
            while True:
                action, *args = self.queue.get_nowait()
                if action == 'add_sample':
                    self._add_sample(*args)
                elif action == 'add_file':
                    self._add_file(*args)
        except queue.Empty:
            pass
        finally:
            self.after(100, self.process_queue)

    def _add_sample(self, sample_name, folder, files):
        if sample_name not in self.samples:
            self.samples[sample_name] = []
            self.tree.insert("", "end", sample_name, text=sample_name)
        
        for file in files:
            file_path = os.path.join(folder, file)
            if file_path not in self.samples[sample_name]:
                self.samples[sample_name].append(file_path)
                self.tree.insert(sample_name, "end", text=file, values=("Existing",))

    def _add_file(self, suggested_sample_name, file_path):
        if suggested_sample_name not in self.samples:
            self.samples[suggested_sample_name] = []
            self.tree.insert("", "end", suggested_sample_name, text=suggested_sample_name)
        
        if file_path not in self.samples[suggested_sample_name]:
            self.samples[suggested_sample_name].append(file_path)
            self.tree.insert(suggested_sample_name, "end", text=os.path.basename(file_path), values=("New",))

        if self.auto_analyze_var.get() and self.should_start_analysis(suggested_sample_name):
            threading.Thread(target=self.start_analysis_for_sample, args=(suggested_sample_name,), daemon=True).start()

    def stop_watching(self):
        if hasattr(self, 'folder_monitor'):  # Changed from self.gui.folder_monitor
            self.folder_monitor.stop()
        self.watching = False
        self.watch_button.configure(text="Start Watching")

    def get_sample_name(self, folder_path):
        if "barcode" in folder_path:
            return os.path.basename(folder_path)
        return os.path.basename(os.path.dirname(folder_path))

    def should_start_analysis(self, sample_name):
        # Add your conditions here, e.g., number of files, time since last analysis, etc.
        return len(self.samples[sample_name]) >= 5  # Start analysis when 5 or more files are available

    def start_analysis_for_sample(self, suggested_sample_name):
        analysis_type = self.analysis_var.get()
        sample_date = self.selected_date.strftime("%Y-%m-%d") if not self.use_file_date.get() else datetime.datetime.fromtimestamp(os.path.getctime(self.samples[suggested_sample_name][0])).strftime("%Y-%m-%d")
        
        if self.use_suggested_name_var.get():
            sample_name = suggested_sample_name
        else:
            sample_name = self.sample_name_entry.get()

        input_files = self.samples[suggested_sample_name]  # Use all files for the sample
        
        cmd = [
            sys.executable,
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", analysis_type,
            "-c", self.gui.db_path or os.path.join(ROOT, "src", "resources", "db_config.json"),
            "-s", sample_name,
            "-d", sample_date,
            "-p", self.project_entry.get(),
            "-u", self.subproject_entry.get(),
            "--overwrite"
        ]
        
        # Add input files to the command
        cmd.extend(["-i"] + input_files)
        
        env = os.environ.copy()
        mmonitor_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        env["PYTHONPATH"] = f"{mmonitor_path}:{env.get('PYTHONPATH', '')}"

        subprocess.Popen(cmd, env=env)

    def get_or_create_concat_file(self, sample_name):
        concat_file = os.path.join(os.path.dirname(self.samples[sample_name][0]), f"{sample_name}_concatenated.fastq.gz")
        
        if not os.path.exists(concat_file):
            # Create new concatenated file
            with gzip.open(concat_file, 'wb') as outfile:
                for file in self.samples[sample_name]:
                    with gzip.open(file, 'rb') as infile:
                        outfile.write(infile.read())
        else:
            # Append new files to existing concatenated file
            with gzip.open(concat_file, 'ab') as outfile:
                for file in self.samples[sample_name]:
                    if os.path.getmtime(file) > os.path.getmtime(concat_file):
                        with gzip.open(file, 'rb') as infile:
                            outfile.write(infile.read())
        
        return concat_file

    def start_analysis(self):
        selected_items = self.tree.selection()
        if not selected_items:
            messagebox.showinfo("Info", "Please select samples to analyze.")
            return

        for item in selected_items:
            suggested_sample_name = self.tree.item(item, "text")
            if suggested_sample_name in self.samples:
                self.start_analysis_for_sample(suggested_sample_name)

    
    def set_date(self, date):
        self.selected_date = date
        self.date_button.configure(text=date.strftime("%Y-%m-%d"))

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