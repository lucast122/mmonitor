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
from requests import post
from tkcalendar import DateEntry
import datetime  # Make sure this import is at the top of your file
import queue
from .FolderWatcherWindow import FolderWatcherWindow
from build_mmonitor_pyinstaller import ROOT, IMAGES_PATH
from mmonitor.dashapp.index import Index
from mmonitor.database.DBConfigForm import DataBaseConfigForm
from mmonitor.database.django_db_interface import DjangoDBInterface
from mmonitor.database.mmonitor_db import MMonitorDBInterface
from .CentrifugeRunner import CentrifugeRunner
from .EmuRunner import EmuRunner
from .FastqStatistics import FastqStatistics
from .InputWindow import InputWindow
from .PipelineWindow import PipelinePopup
from .FunctionalRunner import FunctionalRunner
from .MMonitorCMD import MMonitorCMD

VERSION = "v1.0.0"
MAIN_WINDOW_X, MAIN_WINDOW_Y = 260, 320

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
            ("Watch Folder", self.open_folder_watcher, "mmonitor-folder-watcher.png"),
            ("Toggle Console", self.toggle_console, "mmonitor_button_console.png"),
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

    def run_pipeline(self, analysis_type,emu_db,centriduge_db):
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

    # def start_monitoring(self):
    #     try:
    #         self.dashapp = Index(self.db)
    #         self.monitor_thread = Thread(target=self.dashapp.run_server, args=(False,))
    #         self.monitor_thread.start()
    #     except IndexError:
    #         self.show_info(
    #             "No data found in database. Please first run analysis pipeline to fill DB with data.")
    #         return

    #     # sleep(1)
    #     open_new('http://localhost:8050')

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