import os
import sys
import json
import subprocess
import tkinter as tk
from tkinter import ttk, messagebox
import customtkinter as ctk
from datetime import datetime
from tkcalendar import Calendar
import time
import queue
import datetime
import customtkinter as ctk
from tkinter import messagebox
import tkinter as tk
import multiprocessing
import io
import webbrowser

from mmonitor.userside.utils import create_tooltip
from mmonitor.userside.LoginWindow import LoginWindow
from mmonitor.userside.FolderWatcherWindow import FolderWatcherWindow
from mmonitor.userside.PipelineConfig import PipelineConfig, YAML_CONFIG_PATH
from mmonitor.userside.DatabaseWindow import DatabaseWindow
from mmonitor.userside.FastqStatistics import FastqStatistics
from mmonitor.userside.InputWindow import InputWindow
from mmonitor.pipeline.runner import SnakemakeRunner

from mmonitor.database.DBConfigForm import DataBaseConfigForm
from mmonitor.database.django_db_interface import DjangoDBInterface
from mmonitor.database.mmonitor_db import MMonitorDBInterface

from mmonitor.dashapp.index import Index
from mmonitor.config import _MMONITOR_ROOT as ROOT, _RESOURCES
from ..paths import IMAGES_DIR, RESOURCES_DIR, ROOT as _PROJECT_ROOT


import sys
VERSION = "v0.2.1"
MAIN_WINDOW_X, MAIN_WINDOW_Y = 1500, 1000  # Standard window size
CONSOLE_WIDTH = 300  # Fixed console width
SIDEBAR_WIDTH = 220  # Fixed sidebar width
CONTENT_WIDTH = MAIN_WINDOW_X - SIDEBAR_WIDTH - CONSOLE_WIDTH

# Add color scheme constants
COLORS = {
    "primary": "#2563eb",        # Modern royal blue
    "primary_hover": "#1d4ed8",  # Darker blue for hover
    "secondary": "#64748b",      # Slate gray
    "background": "#ffffff",     # White
    "surface": "#f8fafc",        # Light gray background
    "border": "#e2e8f0",        # Light border color
    "text": "#1e293b",          # Dark blue-gray text
    "text_secondary": "#475569", # Secondary text color
    "success": "#059669",       # Green
    "error": "#dc2626",         # Red
    "warning": "#d97706",       # Amber
    "console_bg": "#f8fafc",    # Very light gray for console
    "console_text": "#1e293b",  # Very dark blue-gray for console text
    "button_text": "#ffffff"    # White text for buttons
}

# Modern styling constants
STYLING = {
    "button_height": 36,        # Taller buttons for better touch targets
    "corner_radius": 8,         # Slightly rounded corners
    "border_width": 1,          # Thin borders
    "font_family": "Helvetica",
    "padding": 16,              # Consistent padding
    "shadow": "0 2px 4px rgba(0, 0, 0, 0.1)"  # Subtle shadow for depth
}

def load_config():
    """Load configuration from config.json"""
    config_path = os.path.join(_RESOURCES, "config.json")
    try:
        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                return json.load(f)
        return {
            "host": "mmonitor.org",
            "port": 443,
            "offline_mode": False
        }
    except Exception as e:
        print(f"Error loading config: {e}")
        return {
            "host": "mmonitor.org",
            "port": 443,
            "offline_mode": False
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
        try:
            # Create a lock file to prevent multiple instances
            self._create_lock_file()
            
            print("Initializing GUI...")
            # Initialize tkinter first
            super().__init__()
            
            # Setup basic logging
            import logging
            logger = logging.getLogger('MMonitor.GUI')
            
            logger.info("Starting GUI initialization")
            
            # Check if we're already running - more aggressive check
            # if self._is_already_running():
            #     logger.info("Application is already running, preventing restart...")
            #     sys.exit(0)
                
            # Set environment variable to indicate we're running
            os.environ['MMONITOR_RUNNING'] = '1'
            
            # Load configuration first
            logger.info("Loading configuration")
            self.config = load_config()
            logger.info(f"Loaded config: {self.config}")
            
            # Set default pipeline config (nested Snakemake-compatible schema)
            self.pipeline_config = {
                'analysis_type': 'taxonomy-wgs',
                'threads': 12,
                'memory_gb': 16,
            }
            
            logger.info("Setting up window properties")
            self.title(f"MMonitor {VERSION}")
            
            # Center window on screen
            screen_width = self.winfo_screenwidth()
            screen_height = self.winfo_screenheight()
            x = (screen_width - MAIN_WINDOW_X) // 2
            y = (screen_height - MAIN_WINDOW_Y) // 2
            
            # Set window size and position
            self.geometry(f"{MAIN_WINDOW_X}x{MAIN_WINDOW_Y}+{x}+{y}")
            self.minsize(MAIN_WINDOW_X, MAIN_WINDOW_Y)
            
            logger.info("Setting up UI constants")
            # Define padding and spacing
            self.FRAME_PADDING = 15
            self.WIDGET_SPACING = 8
            
            # Define font sizes
            self.default_font_size = 14
            self.title_font_size = 20
            self.small_font_size = 12
            
            logger.info("Loading appearance mode")
            # Load appearance mode from system and config
            self.config_file = YAML_CONFIG_PATH
            self.load_appearance_mode()
            
            # Define theme colors based on mode
            self.update_theme_colors()
            
            logger.info("Initializing UI state")
            # Initialize variables
            self.console_expanded = False
            self.console_auto_opened = False
            self.folder_monitor = None
            self.folder_watcher_window = None
            self.database_window = None
            self.db_path = os.path.join(_RESOURCES, "db_config.json")
            self.logged_in = False
            self.offline_mode = False
            self.current_user = None
            self.current_window = None
            
            logger.info("Setting up UI components")
            self.setup_variables()
            self.setup_runners()
            self.init_layout()
            
            logger.info("Setting up output redirection")
            # Redirect stdout and stderr
            sys.stdout = TeeOutput(sys.stdout, self)
            sys.stderr = TeeOutput(sys.stderr, self)
            
            logger.info("GUI initialization complete")
            print("Starting MMonitor application...")
            
            logger.info("Showing login screen")
            self.show_login()
            
            # Bind cleanup to window close
            self.protocol("WM_DELETE_WINDOW", self.on_closing)
            
        except Exception as e:
            import traceback
            error_msg = f"Error during GUI initialization: {str(e)}\n{traceback.format_exc()}"
            print(error_msg)
            try:
                logger.error(error_msg)
            except:
                pass
            
            # Try to show error in GUI if possible
            try:
                import tkinter.messagebox as mb
                mb.showerror("Initialization Error", 
                           f"Failed to start MMonitor:\n{str(e)}\n\nCheck the log file for details.")
            except:
                pass
            
            raise  # Re-raise the exception for the main error handler

    def _create_lock_file(self):
        """Create a lock file to prevent multiple instances from running"""
        import os
        import atexit
        import time
        import random
        
        # Define lock file location in user's home directory
        lock_dir = os.path.expanduser("~/.mmonitor")
        os.makedirs(lock_dir, exist_ok=True)
        self.lock_file = os.path.join(lock_dir, "mmonitor.lock")
        
        # Check if the lock file already exists
        if os.path.exists(self.lock_file):
            # Read the PID from the lock file
            try:
                with open(self.lock_file, 'r') as f:
                    content = f.read().strip()
                    parts = content.split(':')
                    pid = int(parts[0]) if parts else None
                    timestamp = float(parts[1]) if len(parts) > 1 else 0
                    
                    # If lock is older than 60 seconds, assume it's stale
                    if time.time() - timestamp > 60:
                        print(f"Removing stale lock file (age: {time.time() - timestamp:.1f}s)")
                        os.remove(self.lock_file)
                    elif self._is_process_running(pid):
                        print(f"Another instance is already running with PID {pid}")
                        # Small delay to prevent race conditions
                        time.sleep(random.uniform(0.1, 0.5))
                        # Double-check the lock file still exists (maybe it was just removed)
                        if os.path.exists(self.lock_file):
                            print("Exiting due to another running instance")
                            sys.exit(0)
                    else:
                        print(f"Found lock file but process {pid} is not running, removing stale lock")
                        os.remove(self.lock_file)
            except Exception as e:
                print(f"Error checking lock file: {e}")
                # If there's an error reading the lock file, assume it's invalid and remove it
                try:
                    os.remove(self.lock_file)
                except:
                    pass
        
        # Create our lock file with PID and timestamp
        try:
            with open(self.lock_file, 'w') as f:
                f.write(f"{os.getpid()}:{time.time()}")
                
            # Make sure the lock file is removed when the program exits
            atexit.register(self._remove_lock_file)
        except Exception as e:
            print(f"Error creating lock file: {e}")
    
    def _remove_lock_file(self):
        """Remove the lock file when the application exits"""
        try:
            if hasattr(self, 'lock_file') and os.path.exists(self.lock_file):
                os.remove(self.lock_file)
                print("Removed lock file")
        except Exception as e:
            print(f"Error removing lock file: {e}")
    
    def _is_process_running(self, pid):
        """Check if a process with the given PID is running"""
        if pid is None:
            return False
            
        if sys.platform == 'darwin' or sys.platform.startswith('linux'):
            try:
                # Unix-like systems: send signal 0 to check if process exists
                import signal
                os.kill(pid, 0)
                return True
            except OSError:
                return False
        elif sys.platform == 'win32':
            try:
                # Windows
                import ctypes
                kernel32 = ctypes.windll.kernel32
                handle = kernel32.OpenProcess(1, 0, pid)
                if handle == 0:
                    return False
                kernel32.CloseHandle(handle)
                return True
            except:
                return False
        else:
            # Unknown platform
            return False
    
    def _is_already_running(self):
        """Check if the application is already running using multiple methods"""
        # Check if the environment variable is set
        if os.environ.get('MMONITOR_RUNNING') == '1':
            print("MMONITOR_RUNNING environment variable is set")
            return True
        
        # Check the process name for duplicate Python processes with 'MMonitor' in the command line
        import psutil
        current_pid = os.getpid()
        try:
            for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
                try:
                    # Skip if it's the current process
                    if proc.info['pid'] == current_pid:
                        continue
                        
                    # Check the command line for 'MMonitor'
                    cmdline = proc.info.get('cmdline')
                    if cmdline:  # Only check if cmdline is not None
                        if any('MMonitor' in str(cmd) for cmd in cmdline):
                            # Additional check to ensure it's not just a temporary process
                            try:
                                # Try to get process creation time
                                create_time = proc.create_time()
                                # If process is less than 2 seconds old, ignore it
                                if time.time() - create_time < 2:
                                    continue
                                print(f"Found existing MMonitor process: {proc.info['pid']}")
                                return True
                            except (psutil.NoSuchProcess, psutil.AccessDenied):
                                continue
                except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess, KeyError, AttributeError):
                    continue
        except Exception as e:
            print(f"Error checking processes: {e}")
        
        # Check for GUI processes (macOS specific)
        if sys.platform == 'darwin':
            try:
                # Get our own process info
                own_process = psutil.Process(current_pid)
                own_create_time = own_process.create_time()
                
                # Try to use AppleScript to check for running GUI applications
                import subprocess
                proc = subprocess.run([
                    'osascript', 
                    '-e', 
                    'tell application "System Events" to get every process whose name contains "MMonitor"'
                ], capture_output=True, text=True)
                
                if proc.returncode == 0 and proc.stdout.strip():
                    # Parse the output to get PIDs
                    for line in proc.stdout.strip().split('\n'):
                        try:
                            # Try to get the process
                            other_proc = psutil.Process(int(line))
                            # Skip if it's our own process or if it's too new
                            if other_proc.pid == current_pid or time.time() - other_proc.create_time() < 2:
                                continue
                            print(f"Found existing MMonitor GUI process: {other_proc.pid}")
                            return True
                        except (psutil.NoSuchProcess, psutil.AccessDenied, ValueError):
                            continue
            except Exception as e:
                print(f"Error checking GUI processes: {e}")
            
        return False

    def on_closing(self):
        """Clean up resources before closing"""
        try:
            if hasattr(self, 'django_process'):
                print("Shutting down Django server...")
                self.django_process.terminate()
                self.django_process.wait(timeout=5)  # Wait up to 5 seconds for clean shutdown
        except Exception as e:
            print(f"Error during cleanup: {e}")
        finally:
            self.quit()

    def setup_variables(self):
        """Initialize variables and configurations"""
        # Use separate config files for different purposes
        self.db_config_file = os.path.join(_RESOURCES, "db_config.json")
        self.pipeline_config_file = YAML_CONFIG_PATH
        
        # Initialize database interface with login config
        self.django_db = DjangoDBInterface(self.db_config_file)
        
        # Initialize other variables
        self.logged_in = False
        self.offline_mode = False
        self.current_user = None

    def setup_runners(self):
        print("Setting up runners...")
        self.snakemake_runner = SnakemakeRunner()
        print("Runners setup complete.")

    def init_layout(self):
        print("Initializing layout...")
        
        # Set the grey theme using our custom theme file
        if getattr(sys, 'frozen', False):
            # If we're in the bundled app
            theme_path = os.path.join(sys._MEIPASS, "mmonitor", "resources", "grey_theme.json")
        else:
            # Use package-internal resources path
            theme_path = os.path.join(_RESOURCES, "themes", "grey_theme.json")
            if not os.path.exists(theme_path):
                # Fallback: grey_theme.json directly in resources dir
                theme_path = os.path.join(_RESOURCES, "grey_theme.json")
            
        print(f"Loading theme from: {theme_path}")
        
        try:
            # Load and parse the theme file
            with open(theme_path, 'r') as f:
                theme_data = json.load(f)
            
            # Set the theme using the parsed data
            ctk.set_default_color_theme(theme_data)
        except Exception as e:
            print(f"Error loading theme: {e}")
            # Fallback to default theme if loading fails
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
        """Create a console window for output"""
        # Create console frame with fixed width
        self.console_frame = ctk.CTkFrame(
            self.main_frame,
            width=CONSOLE_WIDTH,
            height=MAIN_WINDOW_Y,
            border_width=0,
            fg_color="#E0E0E0"
        )
        self.console_frame.pack(side="right", fill="y")
        self.console_frame.pack_propagate(False)

        # Create title bar frame
        title_bar = ctk.CTkFrame(
            self.console_frame,
            border_width=0,
            fg_color="#D0D0D0"
        )
        title_bar.pack(fill="x", padx=0, pady=0)

        # Add title label
        title_label = ctk.CTkLabel(
            title_bar,
            text="Console",
            font=("Helvetica", 13),
            text_color="#202020"
        )
        title_label.pack(side="left", padx=10)

        # Add minimize button
        self.minimize_btn = ctk.CTkButton(
            title_bar,
            text="−",
            width=20,
            height=20,
            command=self.toggle_console,
            font=("Helvetica", 13),
            fg_color="#D0D0D0",
            hover_color="#C0C0C0",
            text_color="#202020",
            border_width=0
        )
        self.minimize_btn.pack(side="right", padx=5)

        # Create console text widget
        self.console = tk.Text(
            self.console_frame,
            wrap=tk.WORD,
            width=40,
            height=20,
            font=("Courier", 10),
            bg="#E0E0E0",
            fg="#202020",
            insertbackground="#202020",
            borderwidth=0,
            highlightthickness=0
        )
        self.console.pack(fill="both", expand=True, padx=0, pady=0)

        # Add scrollbar
        scrollbar = ctk.CTkScrollbar(
            self.console_frame,
            command=self.console.yview,
            border_spacing=0
        )
        scrollbar.pack(side="right", fill="y")
        self.console.configure(yscrollcommand=scrollbar.set)

        # Configure tag for error messages
        self.console.tag_configure("error", foreground="red")
        
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

        nav_button_args = {
            "fg_color": COLORS["primary"],
            "hover_color": COLORS["primary_hover"],
            "text_color": COLORS["button_text"],
            "border_color": "#000000",
            "border_width": STYLING["border_width"],
            "height": STYLING["button_height"],
            "corner_radius": STYLING["corner_radius"],
            "font": (STYLING["font_family"], 14),
            "anchor": "w",
            "compound": "left"
        }
        
        for text, command, icon, tooltip in buttons:
            btn_container = ctk.CTkFrame(nav_frame, fg_color="transparent")
            btn_container.pack(fill="x", pady=5)
            
            # Create modern button with hover effect
            btn = ctk.CTkButton(
                btn_container, 
                text=f"{icon}  {text}",
                command=lambda cmd=command: self.check_login_and_execute(cmd),
                **nav_button_args
            )
            btn.pack(fill="x", padx=5)
            
            # Add hover animation with mode-aware text color
            def on_enter(e, b=btn):
                b.configure(
                    fg_color="#0056b3",
                    text_color="#FFFFFF"
                )
            
            def on_leave(e, b=btn):
                b.configure(
                    fg_color="#007bff",
                    text_color="#FFFFFF"
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
        
        # Add logout button
        self.logout_button = ctk.CTkButton(
            status_frame,
            text="Logout",
            command=self.logout,
            height=35,
            corner_radius=8,
            font=("Helvetica", 12),
            fg_color="#007bff",
            hover_color="#0056b3",
            text_color="#FFFFFF",
            border_color="#000000",
            border_width=1
        )
        self.logout_button.pack(fill="x", pady=(10, 0))
        
        # Initially hide logout button
        if not self.logged_in:
            self.logout_button.pack_forget()
        
        self.toggle_console_button = ctk.CTkButton(
            status_frame, 
            text="Toggle Console",
            command=self.toggle_console,
            height=35,
            corner_radius=8,
            font=("Helvetica", 12),
            fg_color="#007bff",
            hover_color="#0056b3",
            text_color="#FFFFFF",
            border_color="#000000",
            border_width=1
        )
        self.toggle_console_button.pack(fill="x", pady=(10, 0))

        # Update theme colors for dark/light mode
        def update_colors():
            is_dark = self.appearance_mode == "dark"
            bg_color = "#dbdbdb" if is_dark else "white"
            text_color = "white" if is_dark else "black"
            hover_color = "#b4b4b4" if is_dark else "gray80"
            
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
        self.console.insert("end", text)
        self.console.see("end")
        
        # Auto-show console on first update if configured
        if not self.console_expanded and not self.console_auto_opened:
            self.toggle_console()
            self.console_auto_opened = True

    def open_db_config_form(self):
        print("Opening database configuration form...")
        db_config_form = DataBaseConfigForm(master=self)
        self.wait_window(db_config_form)
        print(f"Database configuration updated: {db_config_form.last_config}")

    
    def run_pipeline(self, params):
        """Run the pipeline with the given parameters"""
        try:
            # Get the analysis type from the FolderWatcher if it exists
            analysis_type = None
            if hasattr(self, 'folder_watcher') and self.folder_watcher:
                analysis_type = self.folder_watcher.analysis_type.get()
            else:
                # Default to WGS if no folder watcher
                analysis_type = 'taxonomy-wgs'
                
            print(f"Starting pipeline for analysis type: {analysis_type}")
            
            # Create input window
            input_window = InputWindow(self)
            
            # If user cancels input, return
            if not input_window.result:
                return
                
            # Create thread for pipeline
            if input_window.multi_sample_mode:
                thread = threading.Thread(target=self.run_multi_sample_pipeline, 
                                       args=(analysis_type, input_window.multi_sample_input))
            else:
                thread = threading.Thread(target=self.run_single_sample_pipeline, 
                                       args=(analysis_type, input_window))
            
            thread.start()
            
        except Exception as e:
            print(f"Error starting pipeline: {e}")
            traceback.print_exc()

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

    def _get_snakemake_target(self, analysis_type: str) -> str:
        """Map GUI analysis type string to Snakemake target rule name."""
        mapping = {
            'taxonomy-16s': 'taxonomy_16s',
            'taxonomy-wgs': 'taxonomy_wgs',
            'assembly': 'assembly_full',   # Full pipeline: assembly+binning+annotation+functional+upload
            'functional': 'functional_full',
        }
        return mapping.get(analysis_type, analysis_type)

    def _build_snakemake_config(self, analysis_type, sample_name, project_name,
                                subproject_name, sample_date, files):
        """Build Snakemake config dict from GUI parameters.

        The PipelineConfig widget now stores a nested Snakemake-compatible dict,
        so we start with that and overlay sample-specific fields.
        """
        import yaml as _yaml

        # Load pipeline config from PipelineConfig widget (already nested)
        pipeline_cfg = {}
        if hasattr(self, 'pipeline_config_frame') and self.pipeline_config_frame:
            pipeline_cfg = getattr(self.pipeline_config_frame, 'pipeline_config', {})

        # If the widget config is empty, load from YAML file directly
        if not pipeline_cfg and os.path.exists(YAML_CONFIG_PATH):
            try:
                with open(YAML_CONFIG_PATH) as f:
                    pipeline_cfg = _yaml.safe_load(f) or {}
            except Exception:
                pass

        # Start with a copy of the nested config
        cfg = dict(pipeline_cfg)

        # Add sample-specific fields
        cfg['sample_name'] = sample_name
        cfg['project_name'] = project_name
        cfg['subproject_name'] = subproject_name
        cfg['sample_date'] = sample_date
        cfg['input_fastq'] = list(files)
        cfg['output_dir'] = os.path.join(os.path.expanduser('~'), 'mmonitor_results')

        # Server config for result upload
        db_config = {}
        db_config_path = self.db_path if self.db_path else os.path.join(_RESOURCES, "db_config.json")
        if os.path.exists(db_config_path):
            try:
                with open(db_config_path) as f:
                    db_config = json.load(f)
            except Exception:
                pass

        cfg['server'] = {
            'url': db_config.get('server_url', 'http://localhost:8000'),
            'username': db_config.get('user', ''),
            'upload_results': not getattr(self, 'offline_mode', False),
        }

        # Set per-tool thread counts to match global threads
        t = cfg.get('threads', 4)
        for tool in ('emu', 'centrifuger', 'flye', 'medaka', 'checkm2', 'gtdbtk', 'bakta'):
            if tool in cfg and isinstance(cfg[tool], dict):
                cfg[tool]['threads'] = t

        return cfg

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

        target = self._get_snakemake_target(analysis_type)
        cfg = self._build_snakemake_config(
            analysis_type, sample_name, project_name,
            subproject_name, sample_date, files
        )

        try:
            print("Executing Snakemake pipeline...")
            runner = SnakemakeRunner()
            exit_code = runner.run(
                target, cfg,
                threads=cfg.get('threads', 4),
                callback=lambda line: print(line),
            )
            if exit_code == 0:
                print("Pipeline completed successfully.")
            else:
                print(f"Pipeline failed with exit code {exit_code}")
        except Exception as e:
            print(f"Error running pipeline: {e}")

    def run_multi_sample_pipeline(self, analysis_type, multi_sample_input):
        print("Preparing multi-sample pipeline...")
        target = self._get_snakemake_target(analysis_type)
        runner = SnakemakeRunner()

        num_samples = len(multi_sample_input['sample_names'])
        for i in range(num_samples):
            sample_name = multi_sample_input['sample_names'][i]
            sample_date = multi_sample_input['dates'][i]
            project_name = multi_sample_input['project_names'][i]
            subproject_name = multi_sample_input['subproject_names'][i]
            files = multi_sample_input['file_paths_lists'][i]

            print(f"\n--- Running sample {i+1}/{num_samples}: {sample_name} ---")

            cfg = self._build_snakemake_config(
                analysis_type, sample_name, project_name,
                subproject_name, sample_date, files
            )

            try:
                exit_code = runner.run(
                    target, cfg,
                    threads=cfg.get('threads', 4),
                    callback=lambda line: print(line),
                )
                if exit_code == 0:
                    print(f"Sample {sample_name} completed successfully.")
                else:
                    print(f"Sample {sample_name} failed with exit code {exit_code}")
            except Exception as e:
                print(f"Error running pipeline for {sample_name}: {e}")

        print("Multi-sample pipeline finished.")

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
        project_name = self.project_entry.get() if hasattr(self, 'project_entry') else ''
        subproject_name = self.subproject_entry.get() if hasattr(self, 'subproject_entry') else ''

        target = self._get_snakemake_target(analysis_type)
        cfg = self._build_snakemake_config(
            analysis_type, sample_name, project_name,
            subproject_name, sample_date, [concat_file]
        )

        try:
            print("Executing Snakemake analysis...")
            runner = SnakemakeRunner()
            exit_code = runner.run(
                target, cfg,
                threads=cfg.get('threads', 4),
                callback=lambda line: print(line),
            )
            if exit_code == 0:
                print("Analysis completed successfully.")
            else:
                print(f"Analysis failed with exit code {exit_code}")
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
        self.django_db = DjangoDBInterface(os.path.join(_RESOURCES, "db_config.json"))


    def create_browser_button(self):
        self.browser_button = ctk.CTkButton(
                    self.login_status_label.master,
                    
                    text="Open Dashboard",
                    height=35,
                    corner_radius=8,
                    font=("Helvetica", 12),
                    fg_color="#007bff",
                    hover_color="#0056b3",
                    text_color="#FFFFFF",
                    border_color="#000000",
                    border_width=1,
                    command=lambda: webbrowser.open("http://127.0.0.1:8000/dash/"),
                    
                    
                )

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
                    'emu_db': os.path.join(_RESOURCES, "emu_db"),
                    'centrifuge_db': os.path.join(_RESOURCES, "centrifuge_db"),
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
        self.database_window = DatabaseWindow(self.content_frame, self)
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
        self.theme_colors = {
            "text": "#202020",              # Dark grey text
            "secondary_text": "#404040",     # Medium grey text
            "button": "#007bff",            # Blue button
            "button_hover": "#0056b3",      # Darker blue for hover
            "button_text": "#FFFFFF",       # White button text
            "button_secondary": "#007bff",  # Blue secondary button
            "button_secondary_hover": "#0056b3",
            "background": "#E0E0E0",        # Light grey background
            "sidebar": "#E0E0E0",           # Light grey sidebar
            "frame": "#E0E0E0",             # Light grey frame
            "border": "#D0D0D0",            # Medium grey border
            "input": "#E0E0E0",             # Light grey input
            "hover": "#D0D0D0"              # Light grey hover
        }


    def update_appearance(self):
        # Always use light mode appearance
        self.configure(fg_color="#E0E0E0")
        for widget in self.sidebar.winfo_children():
            if isinstance(widget, ctk.CTkButton):
                widget.configure(fg_color="#007bff", text_color="#FFFFFF", hover_color="#0056b3", border_color="#000000", border_width=1)
            elif isinstance(widget, ctk.CTkLabel):
                widget.configure(text_color="#202020")

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
                fg_color="#E0E0E0",
                border_color="#D0D0D0"
            )
            self.console.configure(
                fg_color="#E0E0E0",
                text_color="#202020"
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
            text_color="#404040"
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
        frame.configure(fg_color="#E0E0E0")  # Use default if not found
        
        # Content container
        content = ctk.CTkFrame(frame, fg_color="transparent")
        content.pack(fill="x", padx=20, pady=15)
        
        # Title with step number
        ctk.CTkLabel(
            content,
            text=title,
            font=("Helvetica", 16, "bold"),
            text_color="#202020"
        ).pack(anchor="w")
        
        # Description
        ctk.CTkLabel(
            content,
            text=description,
            font=("Helvetica", 12),
            text_color="#404040",
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
            fg_color="#007bff",
            hover_color="#0056b3",
            text_color="#FFFFFF",
            border_color="#000000",
            border_width=1
        ).pack(anchor="w")
        
        return frame

    def show_sequencing_monitor(self):
        
        self.clear_content_frame()
        self.folder_watcher_window = FolderWatcherWindow(self.content_frame, self)
        self.folder_watcher_window.pack(fill="both", expand=True)

    def view_results(self):
        if self.logged_in:
            # Construct the URL with login credentials
            base_url = f"https://{self.django_db.host}:{self.django_db.port}"
            results_url = f"{base_url}/dash/"  
            
            # Open the results page in the default browser
            
            webbrowser.open(results_url)
            
            print(f"Opened results page: {results_url}")
        else:
            print("Please log in first to view results.")
            self.open_login_window()

    def show_analysis(self):
        """Show the analysis configuration window"""
        print("Loading pipeline config:", self.pipeline_config)
        self.clear_content_frame()
        
        # Create a frame for the pipeline configuration
        config_frame = ctk.CTkFrame(self.content_frame)
        config_frame.pack(fill="both", expand=True)
        
        # Create the pipeline configuration window
        self.pipeline_popup = PipelineConfig(config_frame, self)
        self.pipeline_popup.pack(fill="both", expand=True)

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
        self.database_window = DatabaseWindow(self.content_frame, self)
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

    def update_login_status(self):
        """Update login status display"""
        if self.offline_mode:
            self.login_status_label.configure(text="Offline Mode")
            self.logout_button.configure(state="normal")
            self.logout_button.pack(side="bottom", pady=(10, 0))
            
            
            self.create_browser_button()    
            self.browser_button.pack(fill="x", pady=(10,0))
            create_tooltip(self.browser_button, "Open Dashboard in Browser")
        elif self.logged_in:
            self.login_status_label.configure(text=f"Logged in as {self.current_user}")
            self.create_browser_button()    
            
            self.browser_button.configure(command=lambda: self.view_results())
            self.browser_button.pack(fill="x", pady=(10,0))
            
            self.logout_button.configure(state="normal")
            self.logout_button.pack(side="bottom",  pady=(10, 0))
        else:
            self.login_status_label.configure(text="Not logged in")
            self.logout_button.pack_forget()
            self.browser_button.pack_forget()

            
            
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
