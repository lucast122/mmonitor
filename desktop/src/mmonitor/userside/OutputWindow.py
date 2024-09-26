import sys
import customtkinter as ctk
import tkinter as tk


class OutputWindow(ctk.CTkToplevel):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)

        # Set title
        self.title("MMonitor Console")
        app = self

        # Create a scrollable textbox
        self.text = tk.Text(app, wrap=tk.NONE, highlightthickness=0)
        self.text.grid(row=0, column=0, sticky="nsew")

        # Create CTk scrollbar
        ctk_textbox_scrollbar = ctk.CTkScrollbar(app, command=self.text.yview)
        ctk_textbox_scrollbar.grid(row=0, column=1, sticky="ns")

        # Connect textbox scroll event to CTk scrollbar
        self.text.configure(yscrollcommand=ctk_textbox_scrollbar.set)

        # Configure grid layout weights
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        # Redirect stdout to this window
        sys.stdout = self

    def write(self, txt):
        # Append text to the text widget
        self.text.insert(tk.END, txt)
        # Auto-scroll to the end
        self.text.see(tk.END)
        # Update the GUI (needed if writing in a loop or another thread)
        self.text.update_idletasks()

    def close(self):
        # Restore stdout to its original state
        sys.stdout = sys.__stdout__
        self.destroy()
