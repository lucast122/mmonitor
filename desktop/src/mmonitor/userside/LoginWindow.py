import customtkinter as ctk

class LoginWindow(ctk.CTkToplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title("Login")
        self.geometry("300x200")

        self.username = ctk.StringVar()
        self.password = ctk.StringVar()
        self.server_ip = ctk.StringVar(value="134.2.78.150")
        self.port = ctk.StringVar(value="443")
        self.logged_in = False

        ctk.CTkLabel(self, text="Username:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.username).pack(pady=5)

        ctk.CTkLabel(self, text="Password:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.password, show="*").pack(pady=5)

        ctk.CTkLabel(self, text="Server IP:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.server_ip).pack(pady=5)

        ctk.CTkLabel(self, text="Port:").pack(pady=5)
        ctk.CTkEntry(self, textvariable=self.port).pack(pady=5)

        ctk.CTkButton(self, text="Login", command=self.login).pack(pady=10)

    def login(self):
        # Here you would typically verify the credentials with your server
        # For now, we'll just set logged_in to True
        self.logged_in = True
        self.destroy()