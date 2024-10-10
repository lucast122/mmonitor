import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
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
import datetime
import queue
from mmonitor.userside.FolderWatcherWindow import FolderWatcherWindow
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
from mmonitor.userside.MMonitorCMD import MMonitorCMD
from mmonitor.userside.DatabaseWindow import DatabaseWindow
from mmonitor.userside.LoginWindow import LoginWindow
import multiprocessing
from mmonitor.userside.utils import create_tooltip
import json

VERSION = "v1.1.0"
MAIN_WINDOW_X, MAIN_WINDOW_Y = 1100, 950
CONSOLE_WIDTH = 300

class TeeOutput(io.StringIO):
    def __init__(self, original_stream, gui):
        super().__init__()
        self.original_stream = original_stream
        self.gui = gui

    def write(self, s):
        self.original_stream.write(s)
        self.gui.update_console(s)
        with open("database_build.log", "a") as log_file:
            log_file.write(s)

    def flush(self):
        self.original_stream.flush()

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
        self.minsize(MAIN_WINDOW_X, MAIN_WINDOW_Y)  # Set minimum window size
        self.resizable(False, False)

        # Define FRAME_PADDING here
        self.FRAME_PADDING = 20

        # Set initial appearance mode
        self.appearance_mode = "light"
        ctk.set_appearance_mode(self.appearance_mode)

        self.console_expanded = False
        self.console_auto_opened = False
        self.folder_monitor = None
        self.folder_watcher_window = None
        self.database_window = None
        self.db_path = os.path.join(ROOT, "src", "resources", "db_config.json")
        self.logged_in = False
        self.offline_mode = False
        self.current_user = None
        self.current_window = None
        self.setup_variables()
        self.setup_runners()
        self.init_layout()

        # Create console before redirecting stdout and stderr
        self.create_console()

        # Redirect stdout and stderr
        sys.stdout = TeeOutput(sys.stdout, self)
        sys.stderr = TeeOutput(sys.stderr, self)

        print("GUI initialization complete.")
        print("Starting MMonitor application...")

        self.show_login()  # Show login screen by default

    def setup_variables(self):
        print("Setting up variables...")
        self.centrifuge_index_path = f"{ROOT}/src/resources/dec22"
        self.analysis_var = tk.StringVar(value="taxonomy-wgs")
        self.sample_files = {}
        self.pipeline_popup = None
        self.django_db = DjangoDBInterface(f"{ROOT}/src/resources/db_config.json")
        self.db = None
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
        
        # Create main content area
        self.content_frame = ctk.CTkFrame(self)
        self.content_frame.pack(side="right", fill="both", expand=True)

        # Create sidebar
        self.sidebar = ctk.CTkFrame(self, width=200, corner_radius=0)
        self.sidebar.pack(side="left", fill="y", padx=0, pady=0)

        # Add logo or app name to sidebar
        logo_label = ctk.CTkLabel(self.sidebar, text="MMonitor", font=("Helvetica", 20, "bold"))
        logo_label.pack(pady=20)

        self.create_sidebar_buttons()

        # Add Toggle Console button
        self.toggle_console_button = ctk.CTkButton(self.sidebar, text="Toggle Console", command=self.toggle_console)
        self.toggle_console_button.pack(side="bottom", pady=10, padx=10, fill="x")

        # Add Toggle Dark Mode button below Toggle Console
        self.toggle_dark_mode_button = ctk.CTkButton(self.sidebar, text="Toggle Dark Mode", command=self.toggle_dark_mode)
        self.toggle_dark_mode_button.pack(side="bottom", pady=10, padx=10, fill="x")

        print("Layout initialization complete.")

    def create_sidebar_buttons(self):
        buttons = [
            ("Home", self.show_home, "home_icon.png"),
            ("Analysis", self.show_analysis, "analysis_icon.png"),
            ("Watch Folder", self.show_folder_watcher, "folder_icon.png"),
            ("Manage Databases", self.show_database_management, "database_icon.png"),
        ]

        for text, command, icon_name in buttons:
            icon = self.load_icon(icon_name, size=(24, 24))  # Slightly larger icons
            btn = ctk.CTkButton(self.sidebar, text=text, command=lambda cmd=command: self.check_login_and_execute(cmd), 
                                image=icon, compound="left", anchor="w", height=40,
                                fg_color="transparent", hover_color=("gray80", "gray30"),
                                font=("Helvetica", 12))
            btn.pack(pady=5, padx=10, fill="x")
            create_tooltip(btn, f"Click to {text.lower()}")

        # Add login status at the top
        self.login_status_label = ctk.CTkLabel(self.sidebar, text="Not logged in", anchor="w", height=40)
        self.login_status_label.pack(pady=10, padx=10, fill="x", side="top")

        # Update appearance for initial mode
        self.update_appearance()

    def load_icon(self, icon_name, size=(20, 20)):
        icon_path = os.path.join(IMAGES_PATH, icon_name)
        if os.path.exists(icon_path):
            return ctk.CTkImage(Image.open(icon_path), size=size)
        else:
            print(f"Warning: Icon {icon_name} not found.")
            return None

    def toggle_console(self):
        if self.console_expanded:
            self.console_frame.pack_forget()
            self.geometry(f"{MAIN_WINDOW_X}x{MAIN_WINDOW_Y}")
        else:
            self.console_frame.pack(side="right", fill="y", expand=False)
            self.geometry(f"{MAIN_WINDOW_X + CONSOLE_WIDTH}x{MAIN_WINDOW_Y}")
        self.console_expanded = not self.console_expanded

    def update_console(self, text):
        self.console_text.insert("end", text)
        self.console_text.see("end")
        if not self.console_expanded and not self.console_auto_opened:
            self.toggle_console()
            self.console_auto_opened = True

    def open_db_config_form(self):
        print("Opening database configuration form...")
        db_config_form = DataBaseConfigForm(master=self)
        self.wait_window(db_config_form)
        print(f"Database configuration updated: {db_config_form.last_config}")

    def checkbox_popup(self):
        print("Opening analysis pipeline configuration...")
        self.pipeline_popup = PipelinePopup(self, self)
        self.wait_window(self.pipeline_popup)

    def run_pipeline(self, analysis_type, emu_db, centrifuge_db):
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

        # Show loading indicator
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
        if self.check_valid_config():
            self.clear_content_frame()
            self.folder_watcher_window = FolderWatcherWindow(self.content_frame, self)
            self.folder_watcher_window.pack(fill="both", expand=True)
        else:
            messagebox.showwarning("Invalid Configuration", "Please set up analysis parameters in the Analysis tab first.")
            self.show_analysis()

    def check_valid_config(self):
        config_file = os.path.join(os.path.expanduser("~"), ".mmonitor_config.json")
        if os.path.exists(config_file):
            with open(config_file, 'r') as f:
                config = json.load(f)
            # Check for essential parameters
            essential_params = ['analysis_type', 'threads', 'min_length', 'min_quality', 'emu_db', 'centrifuge_db', 'min_abundance']
            missing_params = [param for param in essential_params if param not in config or not config[param]]
            if missing_params:
                messagebox.showwarning("Invalid Configuration", f"The following parameters are missing or empty in the configuration:\n{', '.join(missing_params)}\n\nPlease set these parameters in the Analysis tab.")
                return False
            return True
        else:
            # If in offline mode and config doesn't exist, create a default one
            if self.offline_mode:
                default_config = {
                    'analysis_type': 'taxonomy-wgs',
                    'threads': str(multiprocessing.cpu_count()),
                    'min_length': '1000',
                    'min_quality': '10',
                    'emu_db': os.path.join(ROOT, "src", "resources", "emu_db"),
                    'centrifuge_db': os.path.join(ROOT, "src", "resources", "centrifuge_db"),
                    'min_abundance': '0.01'
                }
                with open(config_file, 'w') as f:
                    json.dump(default_config, f, indent=4)
                return True
            else:
                messagebox.showwarning("Configuration Missing", "Configuration file not found. Please set up analysis parameters in the Analysis tab first.")
                return False

    def select_database(self):
        database_path = filedialog.askdirectory(title="Select Database Directory")
        if database_path:
            self.centrifuge_index_path = database_path
            messagebox.showinfo("Database Selected", f"Database directory set to: {database_path}")

    def open_database_window(self):
        self.clear_content_frame()
        self.database_window = DatabaseWindow(self.content_frame)
        self.database_window.pack(fill="both", expand=True)
        self.update()  # Force update of the GUI

    def show_login(self):
        self.clear_content_frame()
        
        login_frame = ctk.CTkFrame(self.content_frame)
        login_frame.pack(padx=20, pady=20, fill="both", expand=True)

        ctk.CTkLabel(login_frame, text="Login to MMonitor", font=("Helvetica", 24, "bold")).pack(pady=(0, 20))

        # Server address
        server_frame = ctk.CTkFrame(login_frame)
        server_frame.pack(fill="x", pady=10)
        ctk.CTkLabel(server_frame, text="Server:", width=100).pack(side="left", padx=(0, 10))
        self.server_entry = ctk.CTkEntry(server_frame)
        self.server_entry.pack(side="left", expand=True, fill="x")
        self.server_entry.insert(0, "mmonitor.org")

        # Port (disabled by default)
        port_frame = ctk.CTkFrame(login_frame)
        port_frame.pack(fill="x", pady=10)
        ctk.CTkLabel(port_frame, text="Port:", width=100).pack(side="left", padx=(0, 10))
        self.port_entry = ctk.CTkEntry(port_frame, state="disabled")
        self.port_entry.pack(side="left", expand=True, fill="x")
        self.port_entry.insert(0, "443")  # Default HTTPS port
        self.port_checkbox = ctk.CTkCheckBox(port_frame, text="Custom Port", command=self.toggle_port_entry)
        self.port_checkbox.pack(side="left", padx=(10, 0))

        # Username
        username_frame = ctk.CTkFrame(login_frame)
        username_frame.pack(fill="x", pady=10)
        ctk.CTkLabel(username_frame, text="Username:", width=100).pack(side="left", padx=(0, 10))
        self.username_entry = ctk.CTkEntry(username_frame)
        self.username_entry.pack(side="left", expand=True, fill="x")

        # Password
        password_frame = ctk.CTkFrame(login_frame)
        password_frame.pack(fill="x", pady=10)
        ctk.CTkLabel(password_frame, text="Password:", width=100).pack(side="left", padx=(0, 10))
        self.password_entry = ctk.CTkEntry(password_frame, show="*")
        self.password_entry.pack(side="left", expand=True, fill="x")

        # Buttons
        button_frame = ctk.CTkFrame(login_frame)
        button_frame.pack(fill="x", pady=(20, 0))
        ctk.CTkButton(button_frame, text="Login", command=self.perform_login).pack(side="left", padx=5, expand=True, fill="x")
        ctk.CTkButton(button_frame, text="Offline Mode", command=self.enter_offline_mode).pack(side="right", padx=5, expand=True, fill="x")

    def toggle_port_entry(self):
        if self.port_checkbox.get():
            self.port_entry.configure(state="normal")
            if not self.port_entry.get():
                self.port_entry.insert(0, "443")  # Default HTTPS port
        else:
            self.port_entry.configure(state="disabled")

    def perform_login(self):
        server = self.server_entry.get()
        port = self.port_entry.get() if self.port_checkbox.get() else "443"
        username = self.username_entry.get()
        password = self.password_entry.get()
        
        # Here you would typically validate the credentials against your backend
        # For this example, we'll just check if both fields are non-empty
        if server and username and password:
            self.logged_in = True
            self.current_user = username
            self.offline_mode = False
            self.update_login_status()
            self.show_home()
        else:
            messagebox.showerror("Login Failed", "Invalid server, username or password")

    def enter_offline_mode(self):
        self.offline_mode = True
        self.update_login_status()
        self.show_home()

    def update_login_status(self):
        if self.logged_in:
            self.login_status_label.configure(text=f"Logged in as {self.current_user}")
        elif self.offline_mode:
            self.login_status_label.configure(text="Offline Mode")
        else:
            self.login_status_label.configure(text="Not logged in")

    def logout(self):
        self.logged_in = False
        self.current_user = None
        self.offline_mode = False
        self.update_login_status()
        self.show_login()

    def toggle_dark_mode(self):
        if self.appearance_mode == "light":
            self.appearance_mode = "dark"
            ctk.set_appearance_mode("dark")
        else:
            self.appearance_mode = "light"
            ctk.set_appearance_mode("light")
        
        self.update_appearance()
        if self.folder_watcher_window:
            self.folder_watcher_window.update_appearance()

    def update_appearance(self):
        if self.appearance_mode == "light":
            self.configure(fg_color="white")
            for widget in self.sidebar.winfo_children():
                if isinstance(widget, ctk.CTkButton):
                    widget.configure(fg_color="#E0E0E0", text_color="black", hover_color="#CCCCCC")
                elif isinstance(widget, ctk.CTkLabel):
                    widget.configure(text_color="black")
        else:
            self.configure(fg_color="gray10")
            for widget in self.sidebar.winfo_children():
                if isinstance(widget, ctk.CTkButton):
                    widget.configure(fg_color="gray20", text_color="white", hover_color="gray30")
                elif isinstance(widget, ctk.CTkLabel):
                    widget.configure(text_color="white")

    def show_home(self):
        self.clear_content_frame()
        
        content_scroll = ctk.CTkScrollableFrame(self.content_frame)
        content_scroll.pack(fill="both", expand=True, padx=self.FRAME_PADDING, pady=self.FRAME_PADDING)
        
        ctk.CTkLabel(content_scroll, text="Welcome to MMonitor", font=("Helvetica", 24, "bold")).pack(pady=20)
        
        if self.logged_in or self.offline_mode:
            explanation = ("MMonitor is a comprehensive tool for metagenomic analysis.\n"
                           "Use the sidebar to navigate between different functionalities.")
            ctk.CTkLabel(content_scroll, text=explanation, wraplength=500).pack(pady=(0, 20))

            actions = [
                ("Start New Analysis", self.show_analysis, "Configure and run a new analysis pipeline"),
                ("View Recent Results", self.view_results, "Check your latest analysis results"),
                ("Manage Databases", self.open_database_window, "Update or configure your databases"),
                ("Watch Folder", self.open_folder_watcher, "Set up automatic analysis for new files"),
            ]

            for title, command, description in actions:
                card = ctk.CTkFrame(content_scroll)
                card.pack(fill="x", padx=10, pady=10)
                
                ctk.CTkLabel(card, text=title, font=("Helvetica", 18, "bold")).pack(pady=5)
                ctk.CTkLabel(card, text=description, wraplength=300).pack(pady=5)
                ctk.CTkButton(card, text="Go", command=command, width=100).pack(pady=10)

            # Add logout button
            ctk.CTkButton(content_scroll, text="Logout", command=self.logout).pack(pady=20)
        else:
            ctk.CTkLabel(content_scroll, text="Please log in to access MMonitor features.").pack(pady=20)
            ctk.CTkButton(content_scroll, text="Go to Login", command=self.show_login).pack(pady=10)

    def view_results(self):
        if self.logged_in:
            # Construct the URL with login credentials
            base_url = f"https://{self.django_db.host}:{self.django_db.port}"
            results_url = f"{base_url}/results/"  # Adjust this URL as needed
            
            # Open the results page in the default browser
            import webbrowser
            webbrowser.open(results_url)
            
            print(f"Opened results page: {results_url}")
        else:
            print("Please log in first to view results.")
            self.open_login_window()

    def show_analysis(self):
        self.clear_content_frame()
        self.pipeline_popup = PipelinePopup(self.content_frame, self)
        self.pipeline_popup.pack(fill="both", expand=True)

    def show_folder_watcher(self):
        if self.check_valid_config():
            self.clear_content_frame()
            self.folder_watcher_window = FolderWatcherWindow(self.content_frame, self)
            self.folder_watcher_window.pack(fill="both", expand=True)
        else:
            messagebox.showwarning("Invalid Configuration", "Please set up analysis parameters in the Analysis tab first.")
            self.show_analysis()

    def show_database_management(self):
        print("Showing Database Management")
        self.clear_content_frame()
        print("Content frame cleared")
        self.database_window = DatabaseWindow(self.content_frame)
        print("DatabaseWindow created")
        self.database_window.pack(fill="both", expand=True)
        print("DatabaseWindow packed")
        self.update()  # Force update of the GUI

    def clear_content_frame(self):
        if self.content_frame.winfo_children():
            print("Clearing content frame")
            for widget in self.content_frame.winfo_children():
                widget.destroy()
            print("Content frame cleared")
        else:
            print("Content frame is already empty")

    def create_console(self):
        self.console_frame = ctk.CTkFrame(self)
        self.console_frame.pack(side="right", fill="y", expand=False)
        self.console_frame.pack_forget()  # Initially hide the console

        self.console_text = ctk.CTkTextbox(self.console_frame, wrap="word", width=300)
        self.console_text.pack(fill="both", expand=True)

    def check_login_and_execute(self, command):
        if self.logged_in or self.offline_mode:
            command()
        else:
            messagebox.showinfo("Login Required", "Please log in or select offline mode to access this feature.")
            self.show_login()

    def run_pipeline(self, analysis_type, params):
        if self.check_valid_config():
            # Existing pipeline execution code
            pass
        else:
            messagebox.showwarning("Invalid Configuration", "Please set up analysis parameters in the Analysis tab first.")
            self.show_analysis()

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