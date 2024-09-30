import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk  # Add ttk import
from threading import Thread
from webbrowser import open_new
from datetime import datetime
from PIL import Image
import customtkinter as ctk
from CTkMessagebox import CTkMessagebox
import subprocess
import csv
import io
import threading
import gzip
from tkcalendar import DateEntry
import datetime  # Make sure this import is at the top of your file
import queue

from build_mmonitor_pyinstaller import ROOT, IMAGES_PATH
from mmonitor.dashapp.index import Index
from mmonitor.database.DBConfigForm import DataBaseConfigForm
from mmonitor.database.django_db_interface import DjangoDBInterface
from mmonitor.database.mmonitor_db import MMonitorDBInterface
from mmonitor.userside.CentrifugeRunner import CentrifugeRunner
from mmonitor.userside.EmuRunner import EmuRunner
from mmonitor.userside.FastqStatistics import FastqStatistics
from mmonitor.userside.InputWindow import InputWindow
from mmonitor.userside.PipelineWindow import PipelinePopup
from mmonitor.userside.FunctionalRunner import FunctionalRunner
from mmonitor.userside.FolderMonitor import FolderMonitor
from mmonitor.userside.MMonitorCMD import MMonitorCMD

VERSION = "v1.0.0"
MAIN_WINDOW_X, MAIN_WINDOW_Y = 260, 300

class ConsoleWindow(ctk.CTkToplevel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title("Console Output")
        self.geometry("600x400")
        self.console_text = scrolledtext.ScrolledText(self, wrap=tk.WORD)
        self.console_text.pack(expand=True, fill='both')
        self.console_text.config(state="disabled")

    def write(self, text):
        self.console_text.config(state="normal")
        self.console_text.insert(tk.END, text)
        self.console_text.see(tk.END)
        self.console_text.config(state="disabled")
        self.lift()
        self.focus_force()

    def flush(self):
        pass

class StdoutRedirector(io.StringIO):
    def __init__(self, gui):
        super().__init__()
        self.gui = gui

    def write(self, string):
        super().write(string)
        self.gui.update_console(string)

class GUI(ctk.CTk):
    def __init__(self):
        print("Initializing GUI...")
        super().__init__()
        self.title(f"MMonitor {VERSION}")
        self.geometry(f"{MAIN_WINDOW_X}x{MAIN_WINDOW_Y}")
        self.minsize(MAIN_WINDOW_X, MAIN_WINDOW_Y)

        self.console_window = None
        self.console_visible = False
        self.folder_monitor = None
        self.folder_watcher_window = None
        self.setup_variables()
        self.setup_runners()
        self.init_layout()

        # Redirect stdout and stderr
        sys.stdout = StdoutRedirector(self)
        sys.stderr = StdoutRedirector(self)

        print("GUI initialization complete.")

    def setup_variables(self):
        print("Setting up variables...")
        self.centrifuge_index_path = f"{ROOT}/src/resources/dec22"
        self.analysis_var = tk.StringVar(value="taxonomy-wgs")
        self.sample_files = {}
        self.pipeline_popup = None
        self.django_db = DjangoDBInterface(f"{ROOT}/src/resources/db_config.json")
        self.db = None
        self.db_path = None
        self.dashapp = None
        self.monitor_thread = None
        print("Variables setup complete.")

    def setup_runners(self):
        print("Setting up runners...")
        self.centrifuge_runner = CentrifugeRunner()
        self.emu_runner = EmuRunner()
        self.functional_analysis_runner = FunctionalRunner()
        self.cmd_runner = MMonitorCMD()
        print("Runners setup complete.")

    def init_layout(self):
        print("Initializing layout...")
        ctk.set_default_color_theme("blue")
        self.create_header()
        self.create_buttons()
        self.create_console_toggle()
        print("Layout initialization complete.")

    def create_header(self):
        header_label = ctk.CTkLabel(self, text=f"Metagenome Monitor {VERSION}", font=("Helvetica", 18))
        header_label.pack(pady=10)

    def create_buttons(self):
        buttons = [
            ("Select User", self.open_db_config_form, "mmonitor_button4_authenticate.png"),
            ("Run Analysis", self.checkbox_popup, "button_add_data2.png"),
            ("Toggle Console", self.toggle_console, "mmonitor_button_console.png"),
            ("Folder Watcher", self.open_folder_watcher, "mmonitor-folder-watcher.png"),
            # ("Select Database", self.select_database, "mmonitor_button_select_db.png"),
            ("Quit", self.stop_app, "mmonitor_button_quit.png")
        ]

        for text, command, icon_name in buttons:
            icon = self.load_icon(icon_name)
            btn = ctk.CTkButton(self, text=text, command=command, image=icon, width=210, height=40)
            btn.pack(pady=5)

    def create_console_toggle(self):
        pass  # This method is no longer needed as the console toggle is now part of the main buttons

    def load_icon(self, icon_name):
        icon_path = os.path.join(IMAGES_PATH, icon_name)
        if os.path.exists(icon_path):
            icon = Image.open(icon_path)
            return ctk.CTkImage(icon, size=(35, 35))
        else:
            print(f"Warning: Icon {icon_name} not found. Using default icon.")
            return None  # Or return a default icon if you have one

    def toggle_console(self):
        if self.console_visible:
            self.close_console()
        else:
            self.open_console()

    def open_console(self):
        if self.console_window is None or not self.console_window.winfo_exists():
            self.console_window = ConsoleWindow(self)
            self.console_window.protocol("WM_DELETE_WINDOW", self.close_console)
        self.console_window.deiconify()
        self.console_visible = True

    def close_console(self):
        if self.console_window and self.console_window.winfo_exists():
            self.console_window.withdraw()
        self.console_visible = False

    def update_console(self, text):
        if not self.console_visible:
            self.open_console()
        if self.console_window:
            self.console_window.write(text)

    def open_db_config_form(self):
        print("Opening database configuration form...")
        db_config_form = DataBaseConfigForm(master=self)
        print(f"Database configuration updated: {db_config_form.last_config}")

    def checkbox_popup(self):
        print("Opening analysis pipeline configuration...")
        self.pipeline_popup = PipelinePopup(self, self)
        self.wait_window(self.pipeline_popup)

    def run_pipeline(self, analysis_type):
        print(f"Starting pipeline for analysis type: {analysis_type}")
        input_window = InputWindow(self, self.emu_runner)
        self.wait_window(input_window)

        if input_window.do_quit:
            print("Pipeline configuration cancelled by user.")
            return

        # Create a new thread for running the pipeline
        if input_window.process_multiple_samples:
            print("Processing multiple samples...")
            thread = threading.Thread(target=self.run_multi_sample_pipeline, args=(analysis_type, input_window.multi_sample_input))
        else:
            print("Processing single sample...")
            thread = threading.Thread(target=self.run_single_sample_pipeline, args=(analysis_type, input_window))

        # Start the thread
        thread.start()

        # Optionally, you can add a loading indicator or progress bar here
        self.show_loading_indicator()

    def show_loading_indicator(self):
        # Create a new window for the loading indicator
        loading_window = ctk.CTkToplevel(self)
        loading_window.title("Processing")
        loading_window.geometry("300x100")

        # Add a label
        label = ctk.CTkLabel(loading_window, text="Analysis is running...\nThis may take a while.")
        label.pack(pady=20)

        # Add a progress bar (indeterminate mode)
        progress = ctk.CTkProgressBar(loading_window, mode="indeterminate")
        progress.pack(pady=10)
        progress.start()

        # Function to check if the thread is still running
        def check_thread():
            if threading.active_count() > 1:  # Main thread + worker thread
                loading_window.after(100, check_thread)
            else:
                loading_window.destroy()
                self.show_info("Analysis complete!")

        # Start checking
        loading_window.after(100, check_thread)

    def show_info(self, message):
        CTkMessagebox(message=message, icon="check", title="Info")

    def run_single_sample_pipeline(self, analysis_type, input_window):
        sample_name = input_window.sample_name
        project_name = input_window.project_name
        subproject_name = input_window.subproject_name
        sample_date = input_window.selected_date.strftime('%Y-%m-%d')
        files = input_window.file_paths_single_sample

        print(f"Running pipeline for sample: {sample_name}")
        print(f"Project: {project_name}")
        print(f"Subproject: {subproject_name}")
        print(f"Sample date: {sample_date}")
        print(f"Number of input files: {len(files)}")

        config_path = self.db_path if self.db_path else os.path.join(ROOT, "src", "resources", "db_config.json")
        
        cmd = [
            sys.executable,
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", analysis_type,
            "-c", config_path,
            "-i", *files,
            "-s", sample_name,
            "-d", sample_date,
            "-p", project_name,
            "-u", subproject_name
        ]
        
        env = os.environ.copy()
        mmonitor_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        env["PYTHONPATH"] = f"{mmonitor_path}:{env.get('PYTHONPATH', '')}"

        try:
            print("Executing pipeline command...")
            result = subprocess.run(cmd, env=env, text=True, capture_output=True)
            print(result.stdout)
            print(result.stderr, file=sys.stderr)
            print("Pipeline command executed successfully.")
        except Exception as e:
            print(f"Error running pipeline: {e}")

    def run_multi_sample_pipeline(self, analysis_type, multi_sample_input):
        print("Preparing multi-sample pipeline...")
        temp_csv_path = "temp_multi_sample.csv"
        with open(temp_csv_path, 'w', newline='') as csvfile:
            fieldnames = ['sample_name', 'date', 'project_name', 'subproject_name', 'sample_folder']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for i in range(len(multi_sample_input['sample_names'])):
                writer.writerow({
                    'sample_name': multi_sample_input['sample_names'][i],
                    'date': multi_sample_input['dates'][i],
                    'project_name': multi_sample_input['project_names'][i],
                    'subproject_name': multi_sample_input['subproject_names'][i],
                    'sample_folder': os.path.dirname(multi_sample_input['file_paths_lists'][i][0])
                })
        print(f"Temporary CSV file created: {temp_csv_path}")

        config_path = self.db_path if self.db_path else os.path.join(ROOT, "src", "resources", "db_config.json")

        cmd = [
            sys.executable,
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", analysis_type,
            "-c", config_path,
            "-m", temp_csv_path
        ]

        env = os.environ.copy()
        mmonitor_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        env["PYTHONPATH"] = f"{mmonitor_path}:{env.get('PYTHONPATH', '')}"

        try:
            print("Executing multi-sample pipeline command...")
            result = subprocess.run(cmd, env=env, text=True, capture_output=True)
            print(result.stdout)
            print(result.stderr, file=sys.stderr)
            print("Multi-sample pipeline command executed successfully.")
        except Exception as e:
            print(f"Error running multi-sample pipeline: {e}")

    def start_app(self):
        print("Starting MMonitor application...")
        self.mainloop()

    def stop_app(self):
        print("Stopping MMonitor application...")
        if self.monitor_thread is not None and self.monitor_thread.is_alive():
            print("Shutting down monitoring server...")
            post('http://localhost:8050/shutdown')
            self.monitor_thread.join()
        if self.folder_monitor:
            print("Stopping folder monitor...")
            self.folder_monitor.stop()
        print("Application shutdown complete.")
        self.destroy()

    def start_monitoring(self):
        try:
            self.dashapp = Index(self.db)
            self.monitor_thread = Thread(target=self.dashapp.run_server, args=(False,))
            self.monitor_thread.start()
        except IndexError:
            self.show_info(
                "No data found in database. Please first run analysis pipeline to fill DB with data.")
            return

        sleep(1)
        open_new('http://localhost:8050')

    def check_and_start_analysis(self, sample_name):
        # Example condition: if the sample has accumulated 10 files
        if len(self.sample_files[sample_name]) >= 10:
            print(f"Starting analysis for sample {sample_name}")
            self.run_analysis_for_sample(sample_name)

    def run_analysis_for_sample(self, sample_name):
        analysis_type = self.analysis_var.get()
        sample_date = self.selected_date.strftime("%Y-%m-%d") if not self.use_file_date.get() else datetime.datetime.fromtimestamp(os.path.getctime(self.samples[sample_name][0])).strftime("%Y-%m-%d")
        
        concat_file = self.get_or_create_concat_file(sample_name)
        
        cmd = [
            sys.executable,
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", analysis_type,
            "-c", self.db_path or os.path.join(ROOT, "src", "resources", "db_config.json"),
            "-i", concat_file,
            "-s", sample_name,
            "-d", sample_date,
            "-p", self.project_entry.get(),
            "-u", self.subproject_entry.get(),
            "--overwrite"
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

    def start_analysis(self):
        selected_items = self.file_tree.selection()
        for item in selected_items:
            sample_name = self.file_tree.item(item, "text")
            self.run_analysis_for_sample(sample_name)

    def get_selected_analysis_type(self):
        return self.analysis_var.get()

    def refresh_file_tree(self):
        # Clear and rebuild the treeview
        self.file_tree.delete(*self.file_tree.get_children())
        for sample_name, files in self.sample_files.items():
            self.file_tree.insert("", tk.END, iid=sample_name, text=sample_name, values=("", ""))
            for file_path, is_new in files:
                status = "New" if is_new else "Existing"
                item_id = self.file_tree.insert(sample_name, tk.END, text=os.path.basename(file_path), values=(file_path, status))
                if is_new:
                    self.file_tree.item(item_id, tags=("new",))

        # Configure tag colors
        self.file_tree.tag_configure("new", foreground="green")

    def create_file_list_window(self):
        self.file_list_window = ctk.CTkToplevel(self)
        self.file_list_window.title("Sequencer Files")
        self.file_list_window.geometry("600x400")
        self.file_list_window.protocol("WM_DELETE_WINDOW", self.hide_file_list_window)
        self.file_list_window.withdraw()  # Hide the window initially

        self.file_tree = ttk.Treeview(self.file_list_window)
        self.file_tree.pack(fill=tk.BOTH, expand=True)
        self.file_tree["columns"] = ("filepath", "status")
        self.file_tree.column("#0", width=200, minwidth=200)
        self.file_tree.column("filepath", width=300, minwidth=200)
        self.file_tree.column("status", width=100, minwidth=100)
        self.file_tree.heading("#0", text="Sample Name")
        self.file_tree.heading("filepath", text="File Path")
        self.file_tree.heading("status", text="Status")

        # Add radio buttons for selecting analysis type
        analysis_frame = ctk.CTkFrame(self.file_list_window)
        analysis_frame.pack(pady=10)
        
        ctk.CTkRadioButton(analysis_frame, text="Taxonomy WGS", variable=self.analysis_var, value="taxonomy-wgs").pack(side=tk.LEFT, padx=10)
        ctk.CTkRadioButton(analysis_frame, text="Taxonomy 16S", variable=self.analysis_var, value="taxonomy-16s").pack(side=tk.LEFT, padx=10)

        # Add a button to start analysis
        self.start_analysis_button = ctk.CTkButton(self.file_list_window, text="Start Analysis", command=self.start_analysis)
        self.start_analysis_button.pack(pady=10)

    def show_file_list_window(self):
        self.file_list_window.deiconify()

    def hide_file_list_window(self):
        self.file_list_window.withdraw()

    def add_new_file(self, file_path, sample_name, is_new):
        if sample_name not in self.sample_files:
            self.sample_files[sample_name] = []
            # Add a new item to the treeview
            self.file_tree.insert("", tk.END, iid=sample_name, text=sample_name, values=("", ""))
        
        self.sample_files[sample_name].append((file_path, is_new))
        # Add the file under the sample in the treeview
        status = "New" if is_new else "Existing"
        item_id = self.file_tree.insert(sample_name, tk.END, text=os.path.basename(file_path), values=(file_path, status))
        
        if is_new:
            self.file_tree.item(item_id, tags=("new",))
        
        # Configure tag colors
        self.file_tree.tag_configure("new", foreground="green")

        # Optionally, check if conditions are met to start analysis
        self.check_and_start_analysis(sample_name)

    def update_db_config_path(self):
        self.django_db = DjangoDBInterface(f"{ROOT}/src/resources/db_config.json")

    def open_folder_watcher(self):
        if self.folder_watcher_window is None or not self.folder_watcher_window.winfo_exists():
            self.folder_watcher_window = FolderWatcherWindow(self)
        else:
            self.folder_watcher_window.lift()

    def select_database(self):
        database_path = filedialog.askdirectory(title="Select Database Directory")
        if database_path:
            self.centrifuge_index_path = database_path
            messagebox.showinfo("Database Selected", f"Database directory set to: {database_path}")

class ToolTip:
    def __init__(self, widget, tip_text):
        self.widget = widget
        self.tip_text = tip_text
        self.tip_window = None

    def show_tip(self, event=None):
        x, y, _, _ = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 25

        self.tip_window = tk.Toplevel(self.widget)
        self.tip_window.wm_overrideredirect(True)
        self.tip_window.wm_geometry(f"+{x}+{y}")

        label = tk.Label(self.tip_window, text=self.tip_text, foreground="black", background="white", relief="solid", borderwidth=1,
                         font=("Helvetica", "14", "normal"))
        label.pack(ipadx=1)

    def hide_tip(self, event=None):
        if self.tip_window:
            self.tip_window.destroy()
            self.tip_window = None

class CustomDatePicker(ctk.CTkToplevel):
    def __init__(self, parent, callback):
        super().__init__(parent)
        self.callback = callback
        self.title("Select Date")
        self.geometry("300x250")

        self.year_var = tk.StringVar(value=str(datetime.datetime.now().year))
        self.month_var = tk.StringVar(value=str(datetime.datetime.now().month))
        self.day_var = tk.StringVar(value=str(datetime.datetime.now().day))

        ctk.CTkLabel(self, text="Year:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.year_var, width=100).pack()

        ctk.CTkLabel(self, text="Month:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.month_var, width=100).pack()

        ctk.CTkLabel(self, text="Day:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.day_var, width=100).pack()

        ctk.CTkButton(self, text="Select", command=self.on_select).pack(pady=20)

    def on_select(self):
        try:
            selected_date = datetime.date(int(self.year_var.get()),
                                          int(self.month_var.get()),
                                          int(self.day_var.get()))
            self.callback(selected_date)
            self.destroy()
        except ValueError:
            messagebox.showerror("Invalid Date", "Please enter a valid date.")

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
        sample_name = f"{parent_folder}_{barcode_folder}"
        self.queue.put(('add_file', sample_name, file_path))

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

    def _add_file(self, sample_name, file_path):
        if sample_name not in self.samples:
            self.samples[sample_name] = []
            self.tree.insert("", "end", sample_name, text=sample_name)
        
        if file_path not in self.samples[sample_name]:
            self.samples[sample_name].append(file_path)
            self.tree.insert(sample_name, "end", text=os.path.basename(file_path), values=("New",))

        if self.auto_analyze_var.get() and self.should_start_analysis(sample_name):
            threading.Thread(target=self.start_analysis_for_sample, args=(sample_name,), daemon=True).start()

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

    def start_analysis_for_sample(self, sample_name):
        analysis_type = self.analysis_var.get()
        sample_date = self.selected_date.strftime("%Y-%m-%d") if not self.use_file_date.get() else datetime.datetime.fromtimestamp(os.path.getctime(self.samples[sample_name][0])).strftime("%Y-%m-%d")
        
        concat_file = self.get_or_create_concat_file(sample_name)
        
        cmd = [
            sys.executable,
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", analysis_type,
            "-c", self.parent.db_path or os.path.join(ROOT, "src", "resources", "db_config.json"),
            "-i", concat_file,
            "-s", self.sample_name_entry.get(),
            "-d", sample_date,
            "-p", self.project_entry.get(),
            "-u", self.subproject_entry.get(),
            "--overwrite"
        ]
        
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
            sample_name = self.tree.item(item, "text")
            if sample_name in self.samples:
                self.start_analysis_for_sample(sample_name)

    def open_date_picker(self):
        CustomDatePicker(self, self.set_date)

    def set_date(self, date):
        self.selected_date = date
        self.date_button.configure(text=date.strftime("%Y-%m-%d"))
