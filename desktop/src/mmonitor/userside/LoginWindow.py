import customtkinter as ctk
from mmonitor.database.django_db_interface import DjangoDBInterface
from CTkMessagebox import CTkMessagebox
import threading
import time

# Define colors here as well for consistency
COLORS = {
    "primary": "#2B6CB0",        # Deep blue
    "primary_hover": "#2C5282",  # Darker blue
    "secondary": "#718096",      # Slate gray
    "background": "#FFFFFF",     # White
    "surface": "#F7FAFC",        # Light gray background
    "border": "#E2E8F0",        # Light border color
    "text": "#2D3748",          # Dark gray text
    "text_secondary": "#4A5568", # Secondary text color
    "button_text": "#FFFFFF"    # White text for buttons
}

class LoginWindow(ctk.CTkFrame):
    def __init__(self, parent, db_path):
        super().__init__(parent)
        self.parent = parent
        self.db_path = db_path
        
        # Get the main window dimensions from parent
        main_window = self.winfo_toplevel()
        
        # Modern styling
        self.configure(fg_color=COLORS["background"])
        
        # Create centered container
        container = ctk.CTkFrame(
            self,
            fg_color=COLORS["surface"],
            corner_radius=10,
            width=600,  # Fixed width
            height=500  # Fixed height
        )
        container.place(relx=0.5, rely=0.5, anchor="center")
        container.pack_propagate(False)
        
        # Variables
        self.username_var = ctk.StringVar()
        self.password_var = ctk.StringVar()
        self.server_var = ctk.StringVar(value="mmonitor.org")
        self.port_var = ctk.StringVar(value="443")
        self.remember_var = ctk.BooleanVar(value=False)
        self.show_password_var = ctk.BooleanVar(value=False)

        self.logged_in = False
        self.db_interface = DjangoDBInterface(self.db_path)

        self.create_widgets()

    def create_widgets(self):
        # Main container with padding
        container = ctk.CTkFrame(self)
        container.pack(expand=True, padx=20, pady=20)
        
        # Center the content
        content_frame = ctk.CTkFrame(container)
        content_frame.pack(expand=True, padx=40, pady=20)

        # Title
        title_label = ctk.CTkLabel(content_frame, 
                                 text="MMonitor Login",
                                 font=("Helvetica", 24, "bold"))
        title_label.pack(pady=(0, 30))

        # Server details
        server_frame = ctk.CTkFrame(content_frame)
        server_frame.pack(fill="x", pady=(0, 15))
        
        ctk.CTkLabel(server_frame, text="Server:", 
                    font=("Helvetica", 12)).pack(anchor="w")
        server_entry = ctk.CTkEntry(server_frame, textvariable=self.server_var,
                                  height=32, placeholder_text="mmonitor.org")
        server_entry.pack(fill="x", pady=(2, 8))

        ctk.CTkLabel(server_frame, text="Port:", 
                    font=("Helvetica", 12)).pack(anchor="w")
        port_entry = ctk.CTkEntry(server_frame, textvariable=self.port_var,
                                height=32, placeholder_text="443")
        port_entry.pack(fill="x", pady=(2, 0))

        # Login details
        login_frame = ctk.CTkFrame(content_frame)
        login_frame.pack(fill="x", pady=15)
        
        ctk.CTkLabel(login_frame, text="Username:", 
                    font=("Helvetica", 12)).pack(anchor="w")
        username_entry = ctk.CTkEntry(login_frame, textvariable=self.username_var,
                                    height=32, placeholder_text="Enter username")
        username_entry.pack(fill="x", pady=(2, 8))

        ctk.CTkLabel(login_frame, text="Password:", 
                    font=("Helvetica", 12)).pack(anchor="w")
        self.password_entry = ctk.CTkEntry(login_frame, textvariable=self.password_var,
                                         show="*", height=32, 
                                         placeholder_text="Enter password")
        self.password_entry.pack(fill="x", pady=(2, 0))

        # Options
        options_frame = ctk.CTkFrame(content_frame)
        options_frame.pack(fill="x", pady=15)
        
        # Show password and Remember me in same row
        checkbox_frame = ctk.CTkFrame(options_frame)
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
        button_frame = ctk.CTkFrame(content_frame)
        button_frame.pack(fill="x", pady=(20, 0))
        
        # Button container for equal spacing
        btn_container = ctk.CTkFrame(button_frame)
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
        self.status_label = ctk.CTkLabel(content_frame, text="",
                                       font=("Helvetica", 11))
        self.status_label.pack(pady=(15, 0))

    def enter_offline_mode(self):
        """Handle offline mode login"""
        if isinstance(self.parent, ctk.CTk):  # Main window
            main_window = self.parent
        else:  # Content frame
            main_window = self.parent.winfo_toplevel()
        
        main_window.logged_in = True
        main_window.offline_mode = True
        main_window.current_user = "Offline"
        main_window.update_login_status()
        main_window.show_home()
        print("Entered offline mode")

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
        if isinstance(self.parent, ctk.CTk):  # Main window
            main_window = self.parent
        else:  # Content frame
            main_window = self.parent.winfo_toplevel()
        
        if success:
            main_window.logged_in = True
            main_window.offline_mode = False
            main_window.current_user = self.username_var.get()
            main_window.update_login_status()
            self.status_label.configure(text="Login successful.")
            main_window.show_home()
        else:
            self.status_label.configure(text="Login failed. Please try again.")
            CTkMessagebox(master=main_window, title="Login Failed", 
                         message="Invalid credentials or server unreachable.")
            self.login_button.configure(state="normal", text="Login")

        print(f"Login {'succeeded' if success else 'failed'}")
