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

from build_mmonitor_pyinstaller import ROOT

# Define colors here as well for consistency
COLORS = {
    "primary": "#2B6CB0",        # Deep blue
    "primary_hover": "#2C5282",  # Darker blue
    "background": "#ffffff",     # White
    "text": "#1A202C",          # Dark text
    "button_text": "#FFFFFF"     # White text for buttons
}

class LoginWindow(ctk.CTkFrame):
    def __init__(self, parent, db_path):
        super().__init__(parent)
        # Use db_config.json for login settings
        self.db_path = os.path.join(ROOT, "src", "resources", "db_config.json")
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
                       
        ctk.CTkCheckBox(checkbox_frame, text="Remember Me",
                       variable=self.remember_var,
                       font=("Helvetica", 11),
                       height=20).pack(side="right", padx=5)

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
            # Get server path from main window
            server_path = os.path.join("/Users", "timo", "Downloads", "home", "minion-computer", "mmonitor_production", "MMonitor", "server")
        
            if not os.path.exists(server_path):
                raise FileNotFoundError(f"Server directory not found at {server_path}")
            
            print(f"Starting Django server at {server_path}")
            
            # Change to server directory
            os.chdir(server_path)
            
            # Start server process
            cmd = [sys.executable, "manage.py", "runserver", "127.0.0.1:8000"]
            self.server_process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1
            )
            
            # Start thread to monitor server output
            threading.Thread(target=self._monitor_server_output, daemon=True).start()
            
            # Wait for server to start
            self._wait_for_server()
            
        except Exception as e:
            print(f"Error starting server: {e}")
            raise

    def _monitor_server_output(self):
        """Monitor server output in a separate thread"""
        while True:
            output = self.server_process.stdout.readline()
            if output:
                print("Server output:", output.strip())
            
            error = self.server_process.stderr.readline()
            if error:
                print("Server error:", error.strip())
            
            if self.server_process.poll() is not None:
                break

    def _wait_for_server(self):
        """Wait for server to start"""
        max_attempts = 30
        attempt = 0
        while attempt < max_attempts:
            try:
                response = requests.get("http://127.0.0.1:8000", timeout=1)
                if response.status_code == 200:
                    print("Server detected as running!")
                    return
            except requests.exceptions.RequestException:
                attempt += 1
                time.sleep(1)
        
        raise TimeoutError("Server failed to start within timeout period")

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
                        message="Would you like to view the results in your browser?\n"
                               "The local server is running at http://127.0.0.1:8000",
                        icon="question",
                        option_1="Yes",
                        option_2="No"
                    )
                    
                    if dialog.get() == "Yes":
                        webbrowser.open("http://127.0.0.1:8000")
                    
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
            self.db_interface = DjangoDBInterface(self.db_path)
            
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
