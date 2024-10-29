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
from mmonitor.userside.CentrifugerRunner import CentrifugerRunner
from mmonitor.userside.EmuRunner import EmuRunner
from mmonitor.userside.FastqStatistics import FastqStatistics
from mmonitor.userside.InputWindow import InputWindow
from mmonitor.userside.PipelineConfig import PipelineConfig
from mmonitor.userside.FunctionalRunner import FunctionalRunner
from mmonitor.userside.MMonitorCMD import MMonitorCMD
from mmonitor.userside.DatabaseWindow import DatabaseWindow
from mmonitor.userside.LoginWindow import LoginWindow
import multiprocessing
from mmonitor.userside.utils import create_tooltip
import json

VERSION = "v0.1.0"
MAIN_WINDOW_X, MAIN_WINDOW_Y = 1500, 1000  # Standard window size
CONSOLE_WIDTH = 300  # Fixed console width
SIDEBAR_WIDTH = 220  # Fixed sidebar width
CONTENT_WIDTH = MAIN_WINDOW_X - SIDEBAR_WIDTH - CONSOLE_WIDTH

# Add color scheme constants
COLORS = {
    "primary": "#2B6CB0",        # Deep blue
    "primary_hover": "#2C5282",  # Darker blue
    "secondary": "#718096",      # Slate gray
    "background": "#FFFFFF",     # White
    "surface": "#F7FAFC",        # Light gray background
    "border": "#E2E8F0",        # Light border color
    "text": "#2D3748",          # Dark gray text
    "text_secondary": "#4A5568", # Secondary text color
    "success": "#38A169",       # Green
    "error": "#E53E3E",         # Red
    "warning": "#D69E2E",       # Yellow
    "console_bg": "#F8FAFC",    # Very light gray for console
    "console_text": "#1A202C",  # Very dark gray for console text
    "button_text": "#FFFFFF"    # White text for buttons
}

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
        
        # Center window on screen
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        x = (screen_width - MAIN_WINDOW_X) // 2
        y = (screen_height - MAIN_WINDOW_Y) // 2
        
        # Set window size and position
        self.geometry(f"{MAIN_WINDOW_X}x{MAIN_WINDOW_Y}+{x}+{y}")
        self.minsize(MAIN_WINDOW_X, MAIN_WINDOW_Y)
        
        # Define padding and spacing
        self.FRAME_PADDING = 15
        self.WIDGET_SPACING = 8
        
        # Define font sizes
        self.default_font_size = 14
        self.title_font_size = 20
        self.small_font_size = 12

        # Load appearance mode from system and config
        self.config_file = os.path.join(ROOT, "src", "resources", "pipeline_config.json")
        self.load_appearance_mode()
        
        # Define theme colors based on mode
        self.update_theme_colors()
        
        # Initialize variables
        self.console_expanded = False
        self.console_auto_opened = False  # Start with console closed

        # Rest of initialization...
        self.folder_monitor = None
        self.folder_watcher_window = None
        self.database_window = None
        self.db_path = os.path.join(ROOT, "src", "resources", "pipeline_config.json")
        self.logged_in = False
        self.offline_mode = False
        self.current_user = None
        self.current_window = None
        self.setup_variables()
        self.setup_runners()
        self.init_layout()

        # Create console but don't show it
        # self.create_console()

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
        self.centrifugr_runner = CentrifugerRunner()
        self.emu_runner = EmuRunner()
        self.functional_analysis_runner = FunctionalRunner()
        self.cmd_runner = MMonitorCMD()
        print("Runners setup complete.")

    def init_layout(self):
        print("Initializing layout...")
        ctk.set_default_color_theme("blue")
        
        # Create main frame with fixed minimum size
        self.main_frame = ctk.CTkFrame(self)
        self.main_frame.pack(fill="both", expand=True)
        
        # Create a container frame for sidebar and content
        self.inner_frame = ctk.CTkFrame(self.main_frame)
        self.inner_frame.pack(side="left", fill="both", expand=True)
        
        # Create sidebar with fixed width
        self.sidebar = ctk.CTkFrame(
            self.inner_frame, 
            width=SIDEBAR_WIDTH, 
            corner_radius=0
        )
        self.sidebar.pack(side="left", fill="y")
        self.sidebar.pack_propagate(False)
        
        # Create content frame with minimum width
        content_min_width = MAIN_WINDOW_X - SIDEBAR_WIDTH - CONSOLE_WIDTH
        self.content_frame = ctk.CTkFrame(
            self.inner_frame,
            width=content_min_width
        )
        self.content_frame.pack(side="left", fill="both", expand=True)
        
        # Create console frame (always on right)
        self.create_console()
        
        # Create sidebar buttons
        self.create_sidebar_buttons()

    def create_console(self):
        """Create the console frame with consistent size"""
        self.console_frame = ctk.CTkFrame(
            self.main_frame, 
            width=CONSOLE_WIDTH,
            height=MAIN_WINDOW_Y,
            fg_color=COLORS["console_bg"]
        )
        self.console_frame.pack(side="right", fill="y")
        self.console_frame.pack_propagate(False)
        
        # Modern title bar
        console_title = ctk.CTkFrame(
            self.console_frame, 
            height=40, 
            fg_color=COLORS["primary"]
        )
        console_title.pack(fill="x", padx=1, pady=(1, 0))
        console_title.pack_propagate(False)
        
        ctk.CTkLabel(
            console_title, 
            text="Console Output", 
            font=("Helvetica", 12, "bold"),
            text_color=COLORS["button_text"]
        ).pack(side="left", padx=10)
        
        # Improved toggle button
        self.minimize_btn = ctk.CTkButton(
            console_title,
            text="−",
            width=30,
            height=30,
            corner_radius=5,
            command=self.toggle_console,
            fg_color=COLORS["primary"],
            hover_color=COLORS["primary_hover"],
            text_color=COLORS["button_text"]
        )
        self.minimize_btn.pack(side="right", padx=5, pady=5)
        
        # Console text area with improved styling
        self.console_text = ctk.CTkTextbox(
            self.console_frame, 
            wrap="word",
            width=CONSOLE_WIDTH - 20,
            height=MAIN_WINDOW_Y - 60,
            fg_color="white",
            text_color=COLORS["console_text"],
            font=("Consolas", 11)
        )
        self.console_text.pack(expand=True, fill="both", padx=10, pady=10)
        
        self.console_expanded = True
        self.console_auto_opened = True

    def toggle_console(self):
        """Toggle the console visibility with consistent size"""
        if self.console_expanded:
            self.console_frame.pack_forget()
            self.minimize_btn.configure(text="□")
        else:
            self.console_frame.pack(side="right", fill="y")
            self.minimize_btn.configure(text="−")
        self.console_expanded = not self.console_expanded

    def create_sidebar_buttons(self):
        buttons = [
            ("Home", self.show_home, "🏠", "Navigate to the home screen"),
            ("Configuration", self.show_analysis, "⚙️", "Configure analysis parameters"),
            ("Watch Folder", self.show_folder_watcher, "📁", "Monitor folders for new sequences"),
            ("Manage Databases", self.show_database_management, "🗄️", "Manage reference databases"),
        ]

        # Create a container for the entire sidebar content
        sidebar_content = ctk.CTkFrame(self.sidebar, fg_color="transparent")
        sidebar_content.pack(fill="both", expand=True, pady=20)

        # Logo/Title section
        logo_frame = ctk.CTkFrame(sidebar_content, fg_color="transparent")
        logo_frame.pack(fill="x", pady=(0, 30), padx=20)
        
        logo_label = ctk.CTkLabel(logo_frame, text="Metagenome\nMonitor",
                                 font=("Helvetica", 28, "bold"), justify="left")
        logo_label.pack(pady=(0, 5))
        
        version_label = ctk.CTkLabel(logo_frame, text=VERSION,
                                    font=("Helvetica", 12))
        version_label.pack()

        # Navigation buttons with modern styling
        nav_frame = ctk.CTkFrame(sidebar_content, fg_color="transparent")
        nav_frame.pack(fill="x", expand=True, padx=10)

        for text, command, icon, tooltip in buttons:
            btn_container = ctk.CTkFrame(nav_frame, fg_color="transparent")
            btn_container.pack(fill="x", pady=5)
            
            # Create modern button with hover effect
            btn = ctk.CTkButton(
                btn_container, 
                text=f"{icon}  {text}",
                command=lambda cmd=command: self.check_login_and_execute(cmd),
                height=45,
                corner_radius=10,
                font=("Helvetica", 14),
                anchor="w",
                fg_color="transparent",
                text_color=self.theme_colors["text"],
                hover_color=self.theme_colors["button_hover"],
                border_width=2,
                border_color=self.theme_colors["button_secondary"],
                compound="left"
            )
            btn.pack(fill="x", padx=5)
            
            # Add hover animation with mode-aware text color
            def on_enter(e, b=btn):
                is_dark = self.appearance_mode == "dark"
                b.configure(
                    fg_color=self.theme_colors["button"],
                    text_color="white" if is_dark else "black"  # Mode-specific text color
                )
            
            def on_leave(e, b=btn):
                b.configure(
                    fg_color="transparent",
                    text_color=self.theme_colors["text"]
                )
            
            btn.bind("<Enter>", on_enter)
            btn.bind("<Leave>", on_leave)
            create_tooltip(btn, tooltip)

        # Status section at bottom
        status_frame = ctk.CTkFrame(sidebar_content, fg_color="transparent")
        status_frame.pack(fill="x", side="bottom", pady=20, padx=15)
        
        self.login_status_label = ctk.CTkLabel(
            status_frame, 
            text="Not logged in",
            font=("Helvetica", 12)
        )
        self.login_status_label.pack(fill="x")
        
        # Modern toggle console button
        self.toggle_console_button = ctk.CTkButton(
            status_frame, 
            text="Toggle Console",
            command=self.toggle_console,
            height=35,
            corner_radius=8,
            font=("Helvetica", 12),
            fg_color=self.theme_colors["button"],
            hover_color=self.theme_colors["button_hover"],
            border_width=1,
            border_color=self.theme_colors["button"]
        )
        self.toggle_console_button.pack(fill="x", pady=(10, 0))

        # Update theme colors for dark/light mode
        def update_colors():
            is_dark = self.appearance_mode == "dark"
            bg_color = "gray20" if is_dark else "white"
            text_color = "white" if is_dark else "black"
            hover_color = "gray30" if is_dark else "gray80"
            
            for widget in nav_frame.winfo_children():
                if isinstance(widget, ctk.CTkButton):
                    widget.configure(
                        text_color=text_color,
                        hover_color=hover_color,
                        border_color="gray40" if is_dark else "gray70"
                    )
        
        # Call initially and bind to theme changes
        update_colors()
        self.bind("<<ThemeChanged>>", lambda _: update_colors())

    def load_icon(self, icon_name, size=(20, 20)):
        icon_path = os.path.join(IMAGES_PATH, icon_name)
        if os.path.exists(icon_path):
            return ctk.CTkImage(Image.open(icon_path), size=size)
        else:
            print(f"Warning: Icon {icon_name} not found.")
            return None

    def update_console(self, text):
        """Update console content and auto-show if needed"""
        self.console_text.insert("end", text)
        self.console_text.see("end")
        
        # Auto-show console on first update if configured
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
        """Show login window using the LoginWindow class"""
        self.clear_content_frame()
        
        # Create login window embedded in content frame
        login_frame = LoginWindow(self.content_frame, self.db_path)
        login_frame.pack(fill="both", expand=True)
        
        # Ensure console stays on right if expanded
        if self.console_expanded:
            self.console_frame.pack(side="right", fill="y")

    def enter_offline_mode(self):
        self.logged_in = True
        self.offline_mode = True
        self.current_user = "Offline"
        
        # Start the local Django server
        self.start_local_server()
        
        # Update the configuration for the local server
        self.django_db.set_offline_mode(True)
        
        self.update_login_status()
        self.show_home()
        print("Entered offline mode")

    def start_local_server(self):
        def run_server():
            django_dir = '/Users/timo/Downloads/home/minion-computer/mmonitor_production/MMonitor/server'
            manage_py_path = os.path.join(django_dir, 'manage.py')
            python_interpreter = '/Users/timo/miniconda3/bin/python'
            
            if not os.path.exists(manage_py_path):
                print(f"Error: manage.py not found at {manage_py_path}")
                return

            if not os.path.exists(python_interpreter):
                print(f"Error: Python interpreter not found at {python_interpreter}")
                return

            os.chdir(django_dir)
            try:
                result = subprocess.run(
                    [python_interpreter, "manage.py", "runserver", "127.0.0.1:8000"],
                    check=True,
                    capture_output=True,
                    text=True
                )
                print("Server output:", result.stdout)
            except subprocess.CalledProcessError as e:
                print(f"Error starting Django server: {e}")
                print("Error output:", e.output)
                print("Standard output:", e.stdout)
                print("Standard error:", e.stderr)
            except Exception as e:
                print(f"Unexpected error starting Django server: {e}")

        # Run the server in a separate thread
        server_thread = threading.Thread(target=run_server, daemon=True)
        server_thread.start()
        print("Local server thread started")

    def update_login_status(self):
        if self.offline_mode:
            self.login_status_label.configure(text="Offline Mode")
        elif self.logged_in:
            self.login_status_label.configure(text=f"Logged in as {self.current_user}")
        else:
            self.login_status_label.configure(text="Not logged in")

    def perform_login(self):
        self.login_button.configure(state="disabled", text="Logging in...")
        self.offline_button.configure(state="disabled")
        self.login_status.configure(text="Attempting to log in...")
        
        threading.Thread(target=self._login_thread, daemon=True).start()

    def _login_thread(self):
        try:
            success = self.django_db.login(
                self.username_entry.get(),
                self.password_entry.get(),
                self.server_entry.get(),
                self.port_entry.get(),
                self.remember_var.get()
            )
            self.after(0, self._login_callback, success)
        except Exception as e:
            print(f"Exception during login: {e}")
            self.after(0, self._login_callback, False)

    def _login_callback(self, success):
        if success:
            self.logged_in = True
            self.offline_mode = False
            self.current_user = self.username_entry.get()
            self.update_login_status()
            self.login_status.configure(text="Login successful.")
            self.after(1000, self.show_home)
        else:
            self.login_status.configure(text="Login failed. Please try again.")
            self.login_button.configure(state="normal", text="Login")
            self.offline_button.configure(state="normal")

        print(f"Login {'succeeded' if success else 'failed'}")

    def logout(self):
        self.logged_in = False
        self.current_user = None
        self.offline_mode = False
        self.update_login_status()
        self.db_interface = DjangoDBInterface(self.db_path)  # Reset the DB interface
        self.show_login()

    def load_appearance_mode(self):
        """Load appearance mode from system and config"""
        # Always use light mode for now
        self.appearance_mode = "light"
        ctk.set_appearance_mode("light")
        
        # Commented out dark mode detection
        # try:
        #     with open(self.config_file, 'r') as f:
        #         config = json.load(f)
        #         saved_mode = config.get('appearance_mode')
        # except (FileNotFoundError, json.JSONDecodeError):
        #     saved_mode = None
        # 
        # if saved_mode:
        #     self.appearance_mode = saved_mode
        # else:
        #     # Get system preference
        #     import darkdetect
        #     self.appearance_mode = "dark" if darkdetect.isDark() else "light"
        # 
        # ctk.set_appearance_mode(self.appearance_mode)

    def update_theme_colors(self):
        """Update theme colors based on current mode"""
        # Always use light mode colors
        self.theme_colors = {
            "text": "#1A1A1A",
            "secondary_text": "#666666",
            "button": "#2B7DE9",
            "button_hover": "#1E63C4",
            "button_text": "#FFFFFF",
            "button_secondary": "#E0E0E0",
            "button_secondary_hover": "#CCCCCC",
            "background": "#FFFFFF",
            "sidebar": "#F5F5F5",
            "frame": "#FFFFFF",
            "border": "#E0E0E0",
            "input": "#F8F8F8",
            "hover": "#EEEEEE"
        }

    def toggle_dark_mode(self):
        """Toggle between light and dark mode"""
        # Commented out dark mode toggle functionality
        # new_mode = "light" if self.appearance_mode == "dark" else "dark"
        # self.appearance_mode = new_mode
        # ctk.set_appearance_mode(new_mode)
        # 
        # # Update theme colors
        # self.update_theme_colors()
        # 
        # # Update all windows
        # self.update_all_windows()
        # 
        # # Save to config
        # try:
        #     with open(self.config_file, 'r') as f:
        #         config = json.load(f)
        #     config['appearance_mode'] = new_mode
        #     with open(self.config_file, 'w') as f:
        #         json.dump(config, f, indent=4)
        # except Exception as e:
        #     print(f"Error saving appearance mode: {e}")
        pass

    def update_appearance(self):
        # Always use light mode appearance
        self.configure(fg_color="white")
        for widget in self.sidebar.winfo_children():
            if isinstance(widget, ctk.CTkButton):
                widget.configure(fg_color="#E0E0E0", text_color="black", hover_color="#CCCCCC")
            elif isinstance(widget, ctk.CTkLabel):
                widget.configure(text_color="black")

    def update_all_windows(self):
        """Update appearance of all windows"""
        # Update main window
        self.update_appearance()
        
        # Update folder watcher
        if hasattr(self, 'folder_watcher_window') and self.folder_watcher_window:
            self.folder_watcher_window.update_appearance()
        
        # Update database window
        if hasattr(self, 'database_window') and self.database_window:
            self.database_window.update_appearance()
        
        # Update console
        if hasattr(self, 'console_frame'):
            self.update_console_appearance()

    def update_console_appearance(self):
        """Update console appearance based on current theme"""
        if self.console_expanded:
            self.console_frame.configure(
                fg_color=self.theme_colors["frame"],
                border_color=self.theme_colors["border"]
            )
            self.console_text.configure(
                fg_color=self.theme_colors["input"],
                text_color=self.theme_colors["text"]
            )

    def show_home(self):
        """Show the home screen with vertical layout and clear user guidance"""
        self.clear_content_frame()
        
        # Create main container with padding
        home_frame = ctk.CTkFrame(self.content_frame, fg_color="transparent")
        home_frame.pack(fill="both", expand=True, padx=40, pady=20)
        
        # Title and welcome message
        title_frame = ctk.CTkFrame(home_frame, fg_color="transparent")
        title_frame.pack(fill="x", pady=(0, 30))
        
        ctk.CTkLabel(
            title_frame,
            text="Welcome to MMonitor",
            font=("Helvetica", 28, "bold")
        ).pack(anchor="w")
        
        ctk.CTkLabel(
            title_frame,
            text="Real-time Monitoring and Analysis of Nanopore Metagenome Sequencing Data",
            font=("Helvetica", 14),
            text_color=self.theme_colors["secondary_text"]
        ).pack(anchor="w")
        
        # Create steps container
        steps_frame = ctk.CTkFrame(home_frame, fg_color="transparent")
        steps_frame.pack(fill="x", pady=10)
        
        # Step 1: Analysis Configuration
        step1_frame = self.create_step_frame(
            steps_frame,
            "1. Configure Analysis Parameters",
            "Set up your analysis pipeline parameters including taxonomy settings, "
            "quality thresholds, and processing options.",
            "mmonitor_button_dna.png",
            lambda: self.check_login_and_execute(self.show_analysis),
            "Configure Analysis"
        )
        step1_frame.pack(fill="x", pady=(0, 20))
        
        # Step 2: Database Management
        step2_frame = self.create_step_frame(
            steps_frame,
            "2. Manage Taxonomy Databases",
            "Optional: Build or select custom taxonomy databases for your analysis. "
            "Pre-built databases are available by default.",
            "mmonitor_button_dna.png",
            lambda: self.check_login_and_execute(self.show_database_management),
            "Manage Databases"
        )
        step2_frame.pack(fill="x", pady=(0, 20))
        
        # Step 3: Sample Analysis
        step3_frame = self.create_step_frame(
            steps_frame,
            "3. Analyze Sequencing Data",
            "Process your sequencing data in real-time during runs or analyze "
            "completed sequencing data. Supports both single-sample and batch processing.",
            "mmonitor_button_dna.png",
            lambda: self.check_login_and_execute(self.show_sequencing_monitor),
            "Start Analysis"
        )
        step3_frame.pack(fill="x", pady=(0, 20))

    def create_step_frame(self, parent, title, description, icon_path, command, button_text):
        """Create a consistent frame for each step"""
        frame = ctk.CTkFrame(parent)
        frame.configure(fg_color=self.theme_colors.get("surface", "#F7FAFC"))  # Use default if not found
        
        # Content container
        content = ctk.CTkFrame(frame, fg_color="transparent")
        content.pack(fill="x", padx=20, pady=15)
        
        # Title with step number
        ctk.CTkLabel(
            content,
            text=title,
            font=("Helvetica", 16, "bold"),
            text_color=self.theme_colors.get("text", "#2D3748")
        ).pack(anchor="w")
        
        # Description
        ctk.CTkLabel(
            content,
            text=description,
            font=("Helvetica", 12),
            text_color=self.theme_colors.get("secondary_text", "#4A5568"),
            wraplength=800,
            justify="left"
        ).pack(anchor="w", pady=(5, 10))
        
        # Action button
        ctk.CTkButton(
            content,
            text=button_text,
            command=command,
            height=32,
            font=("Helvetica", 12),
            fg_color=self.theme_colors.get("primary", "#2B6CB0"),
            hover_color=self.theme_colors.get("primary_hover", "#2C5282"),
            text_color="#FFFFFF"
        ).pack(anchor="w")
        
        return frame

    def show_sequencing_monitor(self):
        """Show the sequencing monitor (renamed from folder_watcher)"""
        self.clear_content_frame()
        self.folder_watcher_window = FolderWatcherWindow(self.content_frame, self)
        self.folder_watcher_window.pack(fill="both", expand=True)

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
        
        # Create main configuration frame with modern styling
        config_frame = ctk.CTkFrame(self.content_frame)
        config_frame.pack(fill="both", expand=True, padx=20, pady=20)
        
        # Add title with modern typography
        title_frame = ctk.CTkFrame(config_frame, fg_color="transparent")
        title_frame.pack(fill="x", pady=(0, 20))
        
        ctk.CTkLabel(
            title_frame, 
            text="Analysis Configuration",
            font=("Helvetica", 24, "bold")
        ).pack(anchor="w")
        
        ctk.CTkLabel(
            title_frame,
            text="Configure your analysis parameters and pipeline settings",
            font=("Helvetica", 12),
            text_color=self.theme_colors["secondary_text"]
        ).pack(anchor="w")
        
        # Add pipeline configuration with reference to main window
        self.pipeline_popup = PipelineConfig(config_frame, self)
        self.pipeline_popup.pack(fill="both", expand=True)
        
        # Ensure console stays on right
        if self.console_expanded:
            self.console_frame.pack(side="right", fill="y")

    def show_folder_watcher(self):
        if self.check_valid_config():
            # Clear only the content frame contents
            self.clear_content_frame()
            
            # Ensure main frame and sidebar maintain their sizes
            self.main_frame.pack_propagate(False)
            self.sidebar.pack_propagate(False)
            
            # Create and pack the folder watcher window in the content frame
            self.folder_watcher_window = FolderWatcherWindow(self.content_frame, self)
            self.folder_watcher_window.pack(fill="both", expand=True)
            
            # Ensure proper layout
            self.content_frame.pack(side="right", fill="both", expand=True)
            self.sidebar.pack(side="left", fill="y")
            
            # Update geometry to ensure everything fits
            self.update_idletasks()
            
        else:
            messagebox.showwarning("Invalid Configuration", 
                                 "Please set up analysis parameters in the Analysis tab first.")
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
        """Clear only the content frame while preserving the sidebar and console"""
        # Ensure main frame and sidebar maintain their sizes
        self.main_frame.pack_propagate(False)
        self.sidebar.pack_propagate(False)
        
        # Clear content frame widgets
        for widget in self.content_frame.winfo_children():
            widget.destroy()
        
        # Ensure console stays on the right
        if self.console_expanded:
            self.console_frame.pack(side="right", fill="y")

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






























