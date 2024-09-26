import csv
import os
import tkinter as tk
from tkinter import filedialog

import customtkinter as ctk
from CTkMessagebox import CTkMessagebox
from tkcalendar import Calendar


def get_files_from_folder(folder_path):
    """
    Gets a path to a folder, checks if path contains sequencing files with specified endings and returns list
    containing paths to sequencing files.
    """
    files = []

    # We will use os.walk for recursive search
    for file in os.listdir(folder_path):
        full_path = os.path.join(folder_path, file)
        print(full_path)
        if os.path.isfile(full_path) and file.endswith(
                (".fastq", ".fq", ".fasta", ".fastq.gz")) and "concatenated" not in file:
            files.append(full_path)

    # If no files were found, log an error
    if not files:
        print(f"No sequencing files (.fastq, .fq found at {folder_path}")

    return files


class InputWindow(ctk.CTkToplevel):
    def __init__(self, parent, emu_runner):
        super().__init__()
        self.selected_date = None
        self.sample_name = None
        self.project_name = None
        self.subproject_name = None
        self.parent = parent
        self.emu_runner = emu_runner
        self.file_paths_single_sample = []  # Store full paths of selected files
        self.multi_sample_input = {}  # dictionary that will contain file path lists and all other relevant information
        self.do_quit = False
        # Toplevel window

        font = ctk.CTkFont(family="Helvetica", size=12)

        self.title("Sample Data Input")
        self.geometry("500x800")
        self.minsize(500, 800)
        self.process_multiple_samples = False

        padding_y = 2
        label_width = 25
        ctk.CTkLabel(self, text="Sample Name", font=font).pack(pady=padding_y)
        self.sample_name_entry = ctk.CTkEntry(self)
        self.sample_name_entry.pack(pady=2)

        ctk.CTkLabel(self, text="Project Name", font=font).pack(pady=padding_y)
        self.project_name_entry = ctk.CTkEntry(self)
        self.project_name_entry.pack(pady=padding_y)
        # Trace change in project_name_entry to mirror its content into subproject_name_entry
        self.project_name_entry.bind("<KeyRelease>", self.mirror_project_name)

        ctk.CTkLabel(self, text="Subproject Name", font=font).pack(pady=padding_y)
        self.subproject_name_entry = ctk.CTkEntry(self)
        self.subproject_name_entry.pack(pady=padding_y)

        ctk.CTkLabel(self, text="Sample Date", font=font).pack(pady=padding_y)
        self.date_btn = ctk.CTkButton(self, text="Select Date", command=self.open_calendar)
        self.date_btn.pack(pady=padding_y)

        ctk.CTkLabel(self, text="Files", font=font).pack(pady=padding_y)
        self.file_display = tk.Text(self, width=400, height=15, wrap="word",
                                    state=tk.DISABLED)  # Adjusted size and made read-only
        self.file_display.pack(pady=padding_y)

        self.use_multiplexing = tk.BooleanVar(value=False)  # BooleanVar to store the checkbox value

        ctk.CTkButton(self, text="Add FASTQ from folder", command=self.add_data_single_sample).pack(pady=padding_y)
        ctk.CTkButton(self, text="Add multiple samples from CSV", command=self.load_from_csv).pack(pady=padding_y)
        self.multiplexing_checkbox = ctk.CTkCheckBox(self, text="Use Multiplexing", variable=self.use_multiplexing)
        self.multiplexing_checkbox.pack(pady=padding_y)

        ctk.CTkButton(self, text="Submit", command=self.submit).pack(pady=padding_y)
        ctk.CTkButton(self, text="Quit", command=self.quit).pack(pady=padding_y)

    def mirror_project_name(self, event):
        # Get the current value of project name entry and set it to subproject name entry
        project_name = self.project_name_entry.get()
        self.subproject_name_entry.delete(0, tk.END)
        self.subproject_name_entry.insert(0, project_name)

    def open_calendar(self):
        def on_close():
            selected_date = cal.selection_get()
            self.date_btn.configure(text=selected_date.strftime("%Y-%m-%d"))
            self.selected_date = selected_date
            date_win.destroy()

        date_win = tk.Toplevel(self)
        date_win.title("Select a Date")

        cal = Calendar(date_win, selectmode='day')
        cal.pack(pady=20, padx=20)

        ctk.CTkButton(date_win, text="OK", command=on_close).pack(pady=20)

    def add_data_single_sample(self):
        self.process_multiple_samples = False

        folder = filedialog.askdirectory(initialdir='/', title="Choose directory containing sequencing data")
        files = self.fetch_files_from_folder(folder)
        self.update_file_display(files)

    def submit(self):
        self.sample_name = self.sample_name_entry.get()
        self.project_name = self.project_name_entry.get()
        self.subproject_name = self.subproject_name_entry.get()

        self.destroy()

    def load_from_csv(self):
        self.process_multiple_samples = True
        file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")])

        if not file_path:
            print("No file selected.")
            return

        # Check if CSV is empty
        if os.path.getsize(file_path) == 0:
            print("Selected CSV file is empty.")
            return

        error_messages = []  # list to accumulate error messages

        # Prepare dictionary to hold multiple sample data
        self.multi_sample_input = {
            "file_paths_lists": [],
            "sample_names": [],
            "dates": [],
            "project_names": [],
            "subproject_names": []
        }

        with open(file_path, 'r') as file:
            reader = csv.DictReader(file)

            for row in reader:
                # Check if provided path exists
                if not os.path.exists(row["sample folder"].strip()):
                    error_message = f"Invalid path from CSV: {row['sample folder'].strip()}"
                    print(error_message)
                    error_messages.append(error_message)
                    continue

                # Look for the fastq_pass folder in the provided path and its child directories
                folder_path = None
                for root, dirs, files in os.walk(row["sample folder"].strip()):
                    if "fastq_pass" in dirs:
                        folder_path = os.path.join(root, "fastq_pass")
                        break

                if not folder_path:
                    error_message = f"'fastq_pass' directory not found for path: {row['sample folder'].strip()}"
                    print(error_message)
                    error_messages.append(error_message)
                    continue

                # If multiplexing is selected, navigate further to the barcode_x folder
                if self.use_multiplexing.get() == 1:
                    barcode_id_string = str(row['Barcode ID'])

                    if len(barcode_id_string) == 1:
                        barcode_id_string = f"0{barcode_id_string}"

                    barcode_folder = f"barcode{barcode_id_string}"
                    folder_path = os.path.join(folder_path, barcode_folder)

                    if not os.path.exists(folder_path):
                        error_message = f"Barcode folder '{barcode_folder}' not found."
                        print(error_message)
                        error_messages.append(error_message)
                        continue

                files = get_files_from_folder(folder_path)
                print(files)

                # Extract attributes from CSV and check for errors
                required_columns = ["sample_name", "date", "project_name", "subproject_name"]
                for col in required_columns:
                    if not row[col]:
                        error_message = f"Missing {col} in CSV for path: {row['sample folder']}"
                        print(error_message)
                        error_messages.append(error_message)

                # Store extracted CSV information to multi_sample_input
                self.multi_sample_input["file_paths_lists"].append(files)
                self.multi_sample_input["sample_names"].append(row["sample_name"])
                self.multi_sample_input["dates"].append(row["date"])
                self.multi_sample_input["project_names"].append(row["project_name"])
                self.multi_sample_input["subproject_names"].append(row["subproject_name"])

            # Provide feedback on the number of samples to be processed
            num_samples_to_process = len(self.multi_sample_input["file_paths_lists"])
            if num_samples_to_process <= 5:
                self.open_popup(f"Processing {num_samples_to_process} samples.", "Processing multiple samples.",
                                icon="check")
            else:
                self.open_popup(f"Processing {num_samples_to_process} samples. This may take a while.",
                                "Processing multiple samples.", icon="check")

        # Handle error messages from loading the CSV
        if error_messages:
            if len(error_messages) > 2:  # Limit the number of error popups
                self.open_popup("\n".join(error_messages[:2]) + "\n...and more", "Errors with CSV", icon="cancel")
            else:
                for error in error_messages:
                    self.open_popup(error, "Error with CSV", icon="cancel")

            self.destroy()

    def fetch_files_from_folder(self, sample_folder):
        files = [f for f in get_files_from_folder(sample_folder) if
                 f.endswith(("fastq", "fasta", "fq", "fa", 'fastq.gz', "fasta.gz", "fa.gz", "fq.gz"))]

        return files

    def update_file_display(self, files):
        self.file_display.config(state=tk.NORMAL)
        for file in files:
            self.file_paths_single_sample.append(file)
            self.file_display.insert(tk.END, os.path.abspath(file) + "\n")
        self.file_display.config(state=tk.DISABLED)

    def open_popup(self, text, title, icon):
        CTkMessagebox(message=text, title=title, icon=icon, option_1="Okay")

    def quit(self):
        self.do_quit = True
        self.destroy()
