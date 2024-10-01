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

class FolderWatcherWindow(ctk.CTkToplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.title("Folder Watcher")
        self.geometry("800x700")
        self.samples = {}
        self.watching = False
        self.selected_date = datetime.date.today()
        self.auto_analyze_var = tk.BooleanVar(value=False)
        self.use_suggested_name_var = tk.BooleanVar(value=False)
        self.create_widgets()
        
        self.queue = queue.Queue()
        self.after(100, self.process_queue)

    def create_widgets(self):
        self.folder_frame = ctk.CTkFrame(self)
        self.folder_frame.pack(pady=10, padx=10, fill="x")

        self.folder_label = ctk.CTkLabel(self.folder_frame, text="Folder to watch:")
        self.folder_label.pack(side="left", padx=5)

        self.folder_entry = ctk.CTkEntry(self.folder_frame, width=400)
        self.folder_entry.pack(side="left", padx=5)

        self.browse_button = ctk.CTkButton(self.folder_frame, text="Browse", command=self.browse_folder)
        self.browse_button.pack(side="left", padx=5)

        self.watch_button = ctk.CTkButton(self.folder_frame, text="Start Watching", command=self.toggle_watching)
        self.watch_button.pack(side="left", padx=5)

        # Sample information frame
        self.info_frame = ctk.CTkFrame(self)
        self.info_frame.pack(pady=10, padx=10, fill="x")

        # Modify the sample name entry section
        self.sample_name_frame = ctk.CTkFrame(self.info_frame)
        self.sample_name_frame.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky="w")

        self.sample_name_label = ctk.CTkLabel(self.sample_name_frame, text="Sample Name:")
        self.sample_name_label.pack(side="left", padx=5)

        self.sample_name_entry = ctk.CTkEntry(self.sample_name_frame, width=200)
        self.sample_name_entry.pack(side="left", padx=5)

        self.use_suggested_name_check = ctk.CTkCheckBox(self.sample_name_frame, text="Use suggested name", 
                                                        variable=self.use_suggested_name_var, 
                                                        command=self.toggle_sample_name_entry)
        self.use_suggested_name_check.pack(side="left", padx=5)

        self.project_label = ctk.CTkLabel(self.info_frame, text="Project Name:")
        self.project_label.grid(row=1, column=0, padx=5, pady=5)
        self.project_entry = ctk.CTkEntry(self.info_frame, width=200)
        self.project_entry.grid(row=1, column=1, padx=5, pady=5)

        self.subproject_label = ctk.CTkLabel(self.info_frame, text="Subproject Name:")
        self.subproject_label.grid(row=2, column=0, padx=5, pady=5)
        self.subproject_entry = ctk.CTkEntry(self.info_frame, width=200)
        self.subproject_entry.grid(row=2, column=1, padx=5, pady=5)

        self.date_label = ctk.CTkLabel(self.info_frame, text="Date:")
        self.date_label.grid(row=3, column=0, padx=5, pady=5)
        self.date_btn = ctk.CTkButton(self.info_frame, text="Select Date", command=self.open_calendar)
        self.date_btn.grid(row=3, column=1, padx=5, pady=5)

        self.use_file_date = ctk.CTkCheckBox(self.info_frame, text="Use file creation date")
        self.use_file_date.grid(row=3, column=2, padx=5, pady=5)

        # Analysis type selection
        self.analysis_frame = ctk.CTkFrame(self)
        self.analysis_frame.pack(pady=10, padx=10, fill="x")

        self.analysis_var = tk.StringVar(value="taxonomy-wgs")
        ctk.CTkRadioButton(self.analysis_frame, text="WGS Analysis", variable=self.analysis_var, value="taxonomy-wgs").pack(side="left", padx=10)
        ctk.CTkRadioButton(self.analysis_frame, text="16S Analysis", variable=self.analysis_var, value="taxonomy-16s").pack(side="left", padx=10)

        # Create a frame for the Treeview
        tree_frame = ttk.Frame(self)
        tree_frame.pack(pady=10, padx=10, fill="both", expand=True)

        # Create the Treeview widget
        self.tree = ttk.Treeview(tree_frame, columns=("status",))
        self.tree.heading("#0", text="Sample/File")
        self.tree.heading("status", text="Status")
        self.tree.pack(side="left", fill="both", expand=True)

        # Add a scrollbar
        scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        scrollbar.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=scrollbar.set)

        # Manual start analysis button
        self.start_analysis_button = ctk.CTkButton(self, text="Start Analysis", command=self.start_analysis)
        self.start_analysis_button.pack(pady=10)

        # Add auto-analyze checkbox
        self.auto_analyze_check = ctk.CTkCheckBox(self, text="Auto-analyze new samples", variable=self.auto_analyze_var)
        self.auto_analyze_check.pack(pady=10)

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

        self.parent.folder_monitor = FolderMonitor([folder], self)
        self.parent.folder_monitor.start()
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
        if self.parent.folder_monitor:
            self.parent.folder_monitor.stop()
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
            "-c", self.parent.db_path or os.path.join(ROOT, "src", "resources", "db_config.json"),
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