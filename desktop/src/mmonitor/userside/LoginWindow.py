import customtkinter as ctk
from mmonitor.database.django_db_interface import DjangoDBInterface
from CTkMessagebox import CTkMessagebox
import threading

class LoginWindow(ctk.CTkToplevel):
    def __init__(self, parent, db_path):
        super().__init__(parent)
        self.parent = parent
        self.db_path = db_path
        self.title("Login")
        self.geometry("300x450")

        self.username_var = ctk.StringVar()
        self.password_var = ctk.StringVar()
        self.server_var = ctk.StringVar(value="mmonitor.org")
        self.port_var = ctk.StringVar(value="443")
        self.remember_var = ctk.BooleanVar(value=False)
        self.show_password_var = ctk.BooleanVar(value=False)

        self.logged_in = False

        self.create_widgets()

    def create_widgets(self):
        ctk.CTkLabel(self, text="Username:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.username_var).pack(pady=5)

        ctk.CTkLabel(self, text="Password:").pack(pady=5)
        self.password_entry = ctk.CTkEntry(self, textvariable=self.password_var, show="*")
        self.password_entry.pack(pady=5)

        ctk.CTkCheckBox(self, text="Show Password", variable=self.show_password_var, command=self.toggle_password_visibility).pack(pady=5)

        ctk.CTkLabel(self, text="Server:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.server_var).pack(pady=5)

        ctk.CTkLabel(self, text="Port:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.port_var).pack(pady=5)

        ctk.CTkCheckBox(self, text="Remember Me", variable=self.remember_var).pack(pady=10)

        self.login_button = ctk.CTkButton(self, text="Login", command=self.login)
        self.login_button.pack(pady=10)

    def toggle_password_visibility(self):
        if self.show_password_var.get():
            self.password_entry.configure(show="")
        else:
            self.password_entry.configure(show="*")

    def login(self):
        self.login_button.configure(state="disabled", text="Logging in...")
        threading.Thread(target=self._login_thread, daemon=True).start()

    def _login_thread(self):
        db_interface = DjangoDBInterface(self.db_path)
        success = db_interface.login(
            self.username_var.get(),
            self.password_var.get(),
            self.server_var.get(),
            self.port_var.get(),
            self.remember_var.get()
        )
        
        self.after(0, self._login_callback, success)

    def _login_callback(self, success):
        if success:
            self.logged_in = True
            self.parent.logged_in = True
            self.parent.current_user = self.username_var.get()
            self.parent.update_login_status()
            self.destroy()
        else:
            CTkMessagebox(title="Login Failed", message="Invalid credentials or server unreachable.")
        
        # Remove this line as the window might already be destroyed
        # self.login_button.configure(state="normal", text="Login")