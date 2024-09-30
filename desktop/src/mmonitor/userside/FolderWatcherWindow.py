import os
import datetime
import tkinter as tk
from tkinter import filedialog, ttk
import customtkinter as ctk
from mmonitor.userside.FolderMonitor import FolderMonitor
import subprocess
import sys

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
        self.create_widgets()

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

        self.sample_name_label = ctk.CTkLabel(self.info_frame, text="Sample Name:")
        self.sample_name_label.grid(row=0, column=0, padx=5, pady=5)
        self.sample_name_entry = ctk.CTkEntry(self.info_frame, width=200)
        self.sample_name_entry.grid(row=0, column=1, padx=5, pady=5)

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
        self.date_button = ctk.CTkButton(self.info_frame, text="Select Date", command=self.open_date_picker)
        self.date_button.grid(row=3, column=1, padx=5, pady=5)

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

        self.parent.folder_monitor = FolderMonitor([folder], self)
        self.parent.folder_monitor.start()
        self.watching = True
        self.watch_button.configure(text="Stop Watching")

    def stop_watching(self):
        if self.parent.folder_monitor:
            self.parent.folder_monitor.stop()
        self.watching = False
        self.watch_button.configure(text="Start Watching")

    def add_new_file(self, file_path):
        if "_concatenated" in file_path:
            return  # Ignore concatenated files

        parent_folder = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file_path))))
        barcode_folder = os.path.basename(os.path.dirname(file_path))
        sample_name = f"{parent_folder}_{barcode_folder}"

        if sample_name not in self.samples:
            self.samples[sample_name] = []
            self.tree.insert("", "end", sample_name, text=sample_name)
        
        if file_path not in self.samples[sample_name]:
            self.samples[sample_name].append(file_path)
            self.tree.insert(sample_name, "end", text=os.path.basename(file_path), values=("New",))

        if self.auto_analyze_var.get() and self.should_start_analysis(sample_name):
            self.start_analysis_for_sample(sample_name)

    def should_start_analysis(self, sample_name):
        return len(self.samples[sample_name]) >= 5  # Start analysis when 5 or more files are available

    def start_analysis_for_sample(self, sample_name):
        self.parent.run_analysis_for_sample(sample_name)

    def start_analysis(self):
        selected_items = self.tree.selection()
        if not selected_items:
            messagebox.showinfo("Info", "Please select samples to analyze.")
            return

        for item in selected_items:
            sample_name = self.tree.item(item, "text")
            if sample_name in self.samples:
                self.start_analysis_for_sample(sample_name)

    def open_date_picker(self):
        # Implement date picker functionality
        pass

    def get_project_name(self):
        return self.project_entry.get()

    def get_sample_date(self):
        return self.selected_date.strftime("%Y-%m-%d") if not self.use_file_date.get() else None

    def get_subproject_name(self):
        return self.subproject_entry.get()

    def run_analysis_for_sample(self, sample_name):
        analysis_type = self.analysis_var.get()
        sample_date = self.selected_date.strftime("%Y-%m-%d") if not self.use_file_date.get() else datetime.datetime.fromtimestamp(os.path.getctime(self.samples[sample_name][0])).strftime("%Y-%m-%d")
        
        concat_file = self.get_or_create_concat_file(sample_name)
        
        cmd = [
            sys.executable,
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", analysis_type,
            "-c", self.parent.db_path or os.path.join(self.parent.ROOT, "src", "resources", "db_config.json"),
            "-i", concat_file,
            "-s", sample_name,
            "-d", sample_date,
            "-p", self.project_entry.get(),
            "-u", self.subproject_entry.get(),
            "--overwrite",
            "--update"  # Add this to ensure database update
        ]
        
        env = os.environ.copy()
        mmonitor_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        env["PYTHONPATH"] = f"{mmonitor_path}:{env.get('PYTHONPATH', '')}"

        try:
            print("Executing analysis command...")
            result = subprocess.run(cmd, env=env, text=True, capture_output=True)
            print(result.stdout)
            print(result.stderr, file=sys.stderr)
            print("Analysis command executed successfully.")
        except Exception as e:
            print(f"Error running analysis: {e}")