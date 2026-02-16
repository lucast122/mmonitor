import customtkinter as ctk
from mmonitor.database.django_db_interface import DjangoDBInterface
from CTkMessagebox import CTkMessagebox
import threading
import time
import os
import keyring
import sys
import subprocess
import requests
import webbrowser
import io
import fcntl
import os
import tkinter.messagebox as messagebox

from mmonitor.paths import ROOT, SERVER_DIR, RESOURCES_DIR

# Define colors here as well for consistency
COLORS = {
    "primary": "#2B6CB0",        # Deep blue
    "primary_hover": "#2C5282",  # Darker blue
    "background": "#dbdbdb",     # Light grey background
    "text": "#1A202C",          # Dark text
    "button_text": "#FFFFFF"     # White text for buttons
}

class LoginWindow(ctk.CTkFrame):
    def __init__(self, parent, db_path):
        super().__init__(parent)
        # Use db_config.json for login settings
        self.db_path = os.path.join(RESOURCES_DIR, "db_config.json")
        self.parent = parent
        
        
        # Get the main window dimensions from parent
        main_window = self.winfo_toplevel()
        
        # Modern styling
        self.configure(fg_color=COLORS["background"])
        
        # Variables
        self.username_var = ctk.StringVar()
        self.password_var = ctk.StringVar()
        self.server_var = ctk.StringVar(value="mmonitor.org")
        self.port_var = ctk.StringVar(value="443")
        self.remember_var = ctk.BooleanVar(value=False)
        self.show_password_var = ctk.BooleanVar(value=False)

        self.logged_in = False
        self.db_interface = DjangoDBInterface(self.db_path)

        # Store reference to main window
        self.main_window = self.winfo_toplevel()

        self.create_widgets()

    def create_widgets(self):
        # Single clean container
        container = ctk.CTkFrame(self, fg_color=COLORS["background"])
        container.pack(expand=True, padx=60, pady=40)

        # Title
        title_label = ctk.CTkLabel(container, 
                                 text="MMonitor Login",
                                 font=("Helvetica", 24, "bold"))
        title_label.pack(pady=(0, 30))

        # Server details
        ctk.CTkLabel(container, text="Server:", 
                    font=("Helvetica", 12)).pack(anchor="w")
        server_entry = ctk.CTkEntry(container, textvariable=self.server_var,
                                  height=32, placeholder_text="mmonitor.org",
                                  width=300)
        server_entry.pack(fill="x", pady=(2, 15))

        ctk.CTkLabel(container, text="Port:", 
                    font=("Helvetica", 12)).pack(anchor="w")
        port_entry = ctk.CTkEntry(container, textvariable=self.port_var,
                                height=32, placeholder_text="443",
                                width=300)
        port_entry.pack(fill="x", pady=(2, 15))

        # Login details
        ctk.CTkLabel(container, text="Username:", 
                    font=("Helvetica", 12)).pack(anchor="w")
        username_entry = ctk.CTkEntry(container, textvariable=self.username_var,
                                    height=32, placeholder_text="Enter username",
                                    width=300)
        username_entry.pack(fill="x", pady=(2, 15))

        ctk.CTkLabel(container, text="Password:", 
                    font=("Helvetica", 12)).pack(anchor="w")
        self.password_entry = ctk.CTkEntry(container, textvariable=self.password_var,
                                         show="*", height=32, 
                                         placeholder_text="Enter password",
                                         width=300)
        self.password_entry.pack(fill="x", pady=(2, 15))

        # Options
        options_frame = ctk.CTkFrame(container, fg_color="transparent")
        options_frame.pack(fill="x", pady=15)
        
        # Show password and Remember me in same row
        checkbox_frame = ctk.CTkFrame(options_frame, fg_color="transparent")
        checkbox_frame.pack(fill="x")
        
        ctk.CTkCheckBox(checkbox_frame, text="Show Password",
                       variable=self.show_password_var,
                       command=self.toggle_password_visibility,
                       font=("Helvetica", 11),
                       height=20).pack(side="left", padx=5)
                       
        # ctk.CTkCheckBox(checkbox_frame, text="Remember Me",
        #                variable=self.remember_var,
        #                font=("Helvetica", 11),
        #                height=20).pack(side="right", padx=5)

        # Buttons
        button_frame = ctk.CTkFrame(container, fg_color="transparent")
        button_frame.pack(fill="x", pady=(20, 0))
        
        # Button container for equal spacing
        btn_container = ctk.CTkFrame(button_frame, fg_color="transparent")
        btn_container.pack(expand=True)
        
        self.login_button = ctk.CTkButton(btn_container, text="Login",
                                        command=self.login,
                                        font=("Helvetica", 12),
                                        height=32,
                                        width=120)
        self.login_button.pack(side="left", padx=10)
        
        self.offline_button = ctk.CTkButton(btn_container, text="Offline Mode",
                                          command=self.enter_offline_mode,
                                          font=("Helvetica", 12),
                                          height=32,
                                          width=120)
        self.offline_button.pack(side="left", padx=10)

        # Status label
        self.status_label = ctk.CTkLabel(container, text="",
                                       font=("Helvetica", 11))
        self.status_label.pack(pady=(15, 0))

    def show_error(self, message):
        """Show error message in a dialog"""
        messagebox.showerror("Error", message)

    def start_local_server(self):
        """Start the local Django server"""
        try:
            # Get server path based on environment
            if getattr(sys, 'frozen', False):
                # If we're in the bundled app
                server_path = os.path.join(sys._MEIPASS, "server")
            else:
                # If we're in development
                server_path = os.path.join(SERVER_DIR)
            
            # Verify server directory exists
            if not os.path.exists(server_path):
                error_msg = f"Server directory not found at {server_path}"
                print(error_msg)
                raise FileNotFoundError(error_msg)
            
            # Verify manage.py exists
            manage_py_path = os.path.join(server_path, "manage.py")
            if not os.path.exists(manage_py_path):
                error_msg = f"Django manage.py not found at {manage_py_path}"
                print(error_msg)
                raise FileNotFoundError(error_msg)
            
            print(f"Starting Django server at {server_path}")
            print(f"Using Python: {sys.executable}")
            
            # Prepare environment variables
            env = os.environ.copy()
            env["PYTHONUNBUFFERED"] = "1"  # Ensure Python output is unbuffered
            env["DJANGO_SETTINGS_MODULE"] = "MMonitor.settings"  # Ensure settings module is set
            
            # Kill any existing processes on port 8000
            self._kill_processes_on_port(8000)
            
            # Change to server directory
            original_dir = os.getcwd()
            os.chdir(server_path)
            
            # Start server process with the enhanced environment
            cmd = [sys.executable, "manage.py", "runserver", "127.0.0.1:8000"]
            print(f"Executing command: {' '.join(cmd)}")
            
            self.server_process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                env=env
            )
            
            # Start thread to monitor server output
            self.monitor_thread = threading.Thread(target=self._monitor_server_output, daemon=True)
            self.monitor_thread.start()
            
            # Wait for server to start
            self._wait_for_server()
            
            # Return to original directory
            os.chdir(original_dir)
            
        except Exception as e:
            print(f"Error starting server: {e}")
            import traceback
            traceback.print_exc()
            raise
    
    def _kill_processes_on_port(self, port):
        """Attempt to kill any existing processes on the specified port"""
        try:
            if sys.platform == 'win32':
                # Windows
                cmd = f"FOR /F \"tokens=5\" %P IN ('netstat -ano ^| findstr :{port}') DO TaskKill /PID %P /F"
                subprocess.run(cmd, shell=True, capture_output=True)
            else:
                # Unix-like (macOS/Linux)
                cmd = f"lsof -i:{port} | grep LISTEN | awk '{{print $2}}' | xargs -r kill -9"
                subprocess.run(cmd, shell=True, capture_output=True)
            print(f"Attempted to kill processes on port {port}")
        except Exception as e:
            print(f"Warning: Failed to kill processes on port {port}: {e}")

    def _monitor_server_output(self):
        """Monitor server output in a separate thread"""
        print("Starting server output monitoring...")
        
        # Buffer for collecting stderr output
        error_buffer = []
        
        try:
            while True:
                # Check if process has terminated
                if self.server_process.poll() is not None:
                    exit_code = self.server_process.poll()
                    print(f"Django server process terminated with exit code: {exit_code}")
                    if error_buffer:
                        print("Server errors encountered during startup:")
                        for err in error_buffer:
                            print(f"  - {err}")
                    break
                
                # Read stdout
                output = self.server_process.stdout.readline()
                if output:
                    print(f"Django server: {output.strip()}")
                
                # Read stderr
                error = self.server_process.stderr.readline()
                if error:
                    error_str = error.strip()
                    print(f"Django server ERROR: {error_str}")
                    error_buffer.append(error_str)
                
                # Small delay to prevent CPU spinning
                time.sleep(0.1)
        
        except Exception as e:
            print(f"Error monitoring server output: {e}")
            import traceback
            traceback.print_exc()

    def _wait_for_server(self):
        """Wait for server to start"""
        max_attempts = 60  # Increased from 30 to 60 attempts
        attempt = 0
        last_error = None
        
        print(f"Waiting for Django server to start at http://127.0.0.1:8000 (max {max_attempts} attempts)")
        
        while attempt < max_attempts:
            try:
                print(f"Attempt {attempt+1}/{max_attempts} to connect to Django server...")
                response = requests.get("http://127.0.0.1:8000", timeout=2)  # Increased timeout to 2 seconds
                if response.status_code == 200:
                    print("Server detected as running!")
                    return
                else:
                    print(f"Server responded with status code: {response.status_code}")
            except requests.exceptions.RequestException as e:
                last_error = str(e)
                print(f"Connection attempt {attempt+1} failed: {last_error}")
            
            attempt += 1
            time.sleep(1)
        
        error_msg = f"Server failed to start within timeout period. Last error: {last_error}"
        print(error_msg)
        raise TimeoutError(error_msg)

    def enter_offline_mode(self):
        """Enter offline mode with proper server setup"""
        try:
            # Start local server
            self.start_local_server()
            
            def setup_and_configure():
                """Configure application for offline mode"""
                try:
                    print("Setting up offline mode...")
                    # Update application state in the main window
                    self.main_window.logged_in = True
                    self.main_window.offline_mode = True
                    self.main_window.current_user = "offlinemode"
                    
                    if hasattr(self.main_window, 'django_db'):
                        print("\nUpdating database connection...")
                        # First set offline mode
                        self.main_window.django_db.set_offline_mode(True)
                        
                        # Then update connection details
                        self.main_window.django_db.base_url = "http://127.0.0.1:8000"
                        self.main_window.django_db.username = "offlinemode"
                        self.main_window.django_db.password = "offline123"
                        self.main_window.django_db.token = None  # Reset token
                        
                        # Force re-authentication with offline credentials
                        print("Authenticating with offline credentials...")
                        success = self.main_window.django_db.verify_credentials()
                        if not success:
                            print("WARNING: Failed to authenticate with offline server")
                    
                    print("Updating UI...")
                    
                    # Ask user about browser
                    dialog = CTkMessagebox(
                        title="Local Server Ready",
                        message="Would you like to view the dashboard in your browser?\n"
                               "The local server is running at http://127.0.0.1:8000/dashboard/",
                        icon="question",
                        option_1="Yes",
                        option_2="No"
                    )

                    if dialog.get() == "Yes":
                        webbrowser.open("http://127.0.0.1:8000/dashboard/")
                    
                    # Update UI state and show home screen
                    self.main_window.update_login_status()
                    self.main_window.show_home()
                    
                    print("Offline mode setup complete!")
                    
                except Exception as e:
                    print(f"Error setting up offline mode: {e}")
                    self.show_error(str(e))
            
            setup_and_configure()
            
        except Exception as e:
            print(f"Error entering offline mode: {e}")
            self.show_error(str(e))

    def toggle_password_visibility(self):
        if self.show_password_var.get():
            self.password_entry.configure(show="")
        else:
            self.password_entry.configure(show="*")

    def login(self):
        self.login_button.configure(state="disabled", text="Logging in...")
        self.status_label.configure(text="Attempting to log in...")
        threading.Thread(target=self._login_thread, daemon=True).start()

    def _login_thread(self):
        try:
            start_time = time.time()
            
            # Instead of creating a new instance, update the existing one
            # with the current server and port values
            self.db_interface.host = self.server_var.get()
            self.db_interface.port = self.port_var.get()
            self.db_interface.update_base_url()
            
            # Set credentials before login attempt
            self.db_interface.username = self.username_var.get()
            self.db_interface.password = self.password_var.get()
            
            success = self.db_interface.login(
                self.username_var.get(),
                self.password_var.get(),
                self.server_var.get(),
                self.port_var.get(),
                self.remember_var.get()
            )
            end_time = time.time()
            print(f"Login attempt took {end_time - start_time:.2f} seconds")
            self.after(0, self._login_callback, success)
        except Exception as e:
            print(f"Exception during login: {e}")
            self.after(0, self._login_callback, False)

    def _login_callback(self, success):
        """Handle login callback"""
        try:
            # Get the main window instance
            if isinstance(self.parent, ctk.CTk):  # Main window
                main_window = self.parent
            else:  # Content frame
                main_window = self.parent.winfo_toplevel()
            
            if success:
                # Update main window state
                if hasattr(main_window, 'logged_in'):
                    main_window.logged_in = True
                    main_window.offline_mode = False
                    main_window.current_user = self.username_var.get()
                    
                    # Update login status if method exists
                    if hasattr(main_window, 'update_login_status'):
                        main_window.update_login_status()
                    
                    self.status_label.configure(text="Login successful.")
                    
                    # Show home if method exists
                    if hasattr(main_window, 'show_home'):
                        main_window.show_home()
                    
            else:
                self.status_label.configure(text="Login failed. Please try again.")
                CTkMessagebox(
                    title="Login Failed", 
                    message="Invalid credentials or server unreachable.",
                    icon="cancel"
                )
                self.login_button.configure(state="normal", text="Login")

            print(f"Login {'succeeded' if success else 'failed'}")
            
        except Exception as e:
            print(f"Error in login callback: {e}")
            self.status_label.configure(text="Login error occurred.")
            self.login_button.configure(state="normal", text="Login")

    def offline_button_click(self):
        """Handle offline mode button click"""
        print("Switching to offline mode...")
        if self.parent:
            self.enter_offline_mode()  # Call the method we just implemented
        self.destroy()
