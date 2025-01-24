import os
import json
import customtkinter as ctk
from tkinter import messagebox
import tkinter as tk

# Get the src directory path
SRC_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

"""
Class DatabaseConfigForm

Window that lets user chose Host, User and Password to authenticate with the mmonitor webserver

Host: Public IP adress of the mmonitor webserver instance (server part of software)
User: Username on mmonitor webserver (used for registration)
Password: Password on mmonitor webserver (used for registration)

save_config(): saves the config as a json file under f"{SRC_DIR}/resources/db_config.json" 

NOTE: After registration an admin has to unlock account first
"""


class DataBaseConfigForm(ctk.CTkToplevel):
    def __init__(self, master=None):
        super().__init__(master=master)
        self.title("User Authentication")
        self.geometry("300x300")
        self.master = master

        self.db_config = {
            "host": tk.StringVar(),
            "user": tk.StringVar(),
            "password": tk.StringVar()
        }

        # Load last config if it exists
        self.last_config = {}
        if os.path.exists(f"{SRC_DIR}/resources/db_config.json"):
            try:
                with open(f"{SRC_DIR}/resources/db_config.json", "r") as f:
                    self.last_config = json.load(f)
            except json.JSONDecodeError:
                messagebox.showerror("Error", "Couldn't load the previous DB configuration")

        # Fill the form with the last config if it exists
        for key in self.db_config.keys():
            if key in self.last_config:
                self.db_config[key].set(self.last_config[key])

        self.password_censored = tk.BooleanVar(value=False)


        self.create_widgets()

    def create_widgets(self):
        frame = ctk.CTkFrame(self, corner_radius=10)
        frame.pack(pady=5, padx=10, fill="both", expand=True)
        for i, key in enumerate(self.db_config.keys()):
            ctk.CTkLabel(frame, text=key).pack(pady=2)
            ctk.CTkEntry(frame, textvariable=self.db_config[key]).pack(pady=5)
        ctk.CTkCheckBox(frame, text="Censor password", variable=self.password_censored,
                        command=self.toggle_password_censor).pack(pady=5)
        ctk.CTkButton(frame, text="Save Config", command=self.save_config).pack(pady=5)

    def toggle_password_censor(self):
        for widget in self.grid_slaves():
            if isinstance(widget, tk.Entry) and widget.cget('textvariable') == str(self.db_config['password']):
                widget.config(show="*" if self.password_censored.get() else None)

    def password_is_hash(self, password):
        """
        Check if the password is a hash.
        This is a simple check to see if the password is a 64 character string (which is the length of a SHA-256 hash).
        It is not a foolproof method, but will work for this demonstration.
        """
        return len(password) == 64

    def save_config(self):
        config = {key: self.db_config[key].get() for key in self.db_config.keys()}

        # Only hash the password if it's not already a hash
        # if not self.password_is_hash(config["password"]):
        #     config["password"] = hashlib.sha256(config["password"].encode()).hexdigest()
        # try without hashing for now

        with open(f"{SRC_DIR}/resources/db_config.json", "w") as f:
            json.dump(config, f)

        messagebox.showinfo("Success", "DB Config saved successfully")
        self.master.update_db_config_path()
        print("Updated db_config path")

        self.destroy()
