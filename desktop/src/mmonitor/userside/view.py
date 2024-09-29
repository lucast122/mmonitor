import gzip
import hashlib
import json
import os
import sys
import tarfile
import time
import tkinter as tk
import numpy as np
import urllib.request
from datetime import datetime
from threading import Thread
from time import sleep
from tkinter import *
from tkinter import simpledialog
from tkinter import ttk
from webbrowser import open_new
from tkinter import messagebox, scrolledtext
import numpy as npfrom
from PIL import Image
from customtkinter import CTkImage
import CTkMessagebox
from customtkinter import filedialog
from requests import post
from tkcalendar import Calendar
# from mmonitor.Tooltip import ToolTip
from build_mmonitor_pyinstaller import ROOT, IMAGES_PATH
from mmonitor.dashapp.index import Index
from mmonitor.database.DBConfigForm import DataBaseConfigForm
from mmonitor.database.django_db_interface import DjangoDBInterface
from mmonitor.database.mmonitor_db import MMonitorDBInterface
from mmonitor.userside.CentrifugeRunner import CentrifugeRunner
from mmonitor.userside.EmuRunner import EmuRunner
from mmonitor.userside.FastqStatistics import FastqStatistics
from mmonitor.userside.InputWindow import InputWindow
from mmonitor.userside.PipelineWindow import PipelinePopup
from mmonitor.userside.FunctionalRunner import FunctionalRunner
from mmonitor.userside.FolderMonitor import FolderMonitor
from mmonitor.userside.MMonitorCMD import MMonitorCMD
import tkinter as tk
from tkinter import messagebox, scrolledtext
import customtkinter as ctk
import sys
import traceback
import subprocess
import csv

# Global constants for version and dimensions
VERSION = "v1.0.0"
MAIN_WINDOW_X: int = 260
MAIN_WINDOW_Y: int = 600

# Module description

"""
This file represents the basic gui for the desktop app. It is the entry point for the program and the only way 
for the user to create projects, select files and run MMonitor's computational engine (centrifuge at this moment)
"""

def compute_sha256(file_path):
    """Compute the sha256 checksum of a file. Used to check if index files like emu_db index have been downloaded correctly"""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def require_project(func):
    """Decorator that ensures that a database was selected or created by the user."""

    def func_wrapper(*args):
        obj: GUI = args[0]
        if obj.db_path is not None and len(obj.db_path) > 0:
            return func(*args)
        else:
            obj.open_popup("Please first create or choose a local DB.", "No data base chosen", icon="cancel")

    return func_wrapper

class RedirectText:
    def __init__(self, widget):
        self.widget = widget

    def write(self, text):
        self.widget.configure(state='normal')  # Enable the widget to insert text
        self.widget.insert(tk.END, text)
        self.widget.configure(state='disabled')  # Disable the widget to make it read-only
        self.widget.see(tk.END)

    # def flush(self):
    #     passdef



class GUI(ctk.CTk):
    """
    Main GUI class for the MMonitor desktop application. It initializes the GUI layout, provides methods for
    handling user interactions, and interfaces with other components of the application for various tasks.
    """

    """
    Initialize the GUI with default settings, prepare the database interfaces, and set up other essential attributes.
    """

    def __init__(self):
        super().__init__()
        self.centrifuge_index_path = f"{ROOT}/src/resources/dec22"
        # Move this line to the beginning of __init__
        self.analysis_var = tk.StringVar(value="taxonomy-wgs")  # Default to taxonomy-wgs
        
        SHOW_CONSOLE = False
        self.sample_files = {}  # Dictionary to hold sample names and their associated files
                
        self.folder_monitor = None
        self.create_file_list_window()
        # Add a button to select folders to monitor
        # hide MMonitor folder button until it works properly
        # self.monitor_folders_button = ctk.CTkButton(self, text="Monitor Folders", command=self.select_and_monitor_folders)
        # self.monitor_folders_button.pack(pady=10)
        self.create_file_list_window()
        # Create a console for logging
        self.console_frame = ctk.CTkFrame(self)
        self.console_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        
        self.console_text = scrolledtext.ScrolledText(self.console_frame, wrap=tk.WORD, height=10)
        self.console_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.console_scrollbar = ctk.CTkScrollbar(self.console_frame, command=self.console_text.yview)
        self.console_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.console_text.config(yscrollcommand=self.console_scrollbar.set, state="disabled")

        
        # Redirect stdout and stderr to the console
        if SHOW_CONSOLE:
            sys.stdout = RedirectText(self.console_text)
            sys.stderr = RedirectText(self.console_text)


        self.pipeline_popup = None
        self.django_db = DjangoDBInterface(f"{ROOT}/src/resources/db_config.json")
        self.progress_bar_exists = False
        self.input_window = None
        self.download_progress = 0
        self.cmd_runner = MMonitorCMD()

        # self.db_mysql.create_db()

        # declare data base class variable, to be chosen by user with choose_project()
        self.db: MMonitorDBInterface = None
        # self.db_mysql = None
        self.db_path = None
        self.centrifuge_runner = CentrifugeRunner()
        self.emu_runner = EmuRunner()
        self.functional_analysis_runner = FunctionalRunner()
        self.dashapp = None
        self.monitor_thread = None

        # self = ctk.CTk()
        mmonitor_logo = tk.PhotoImage(file=f"{IMAGES_PATH}/mmonitor_logo.png")
        self.iconphoto(True, mmonitor_logo)
        self.init_layout()
        self.taxonomy_nanopore_wgs_bool = tk.BooleanVar()
        self.taxonomy_nanopore_16s_bool = tk.BooleanVar()
        self.assembly = tk.BooleanVar()
        self.correction = tk.BooleanVar()
        self.binning = tk.BooleanVar()
        self.annotation = tk.BooleanVar()
        self.kegg = tk.BooleanVar()
        self.sample_date = None

        # Add this line to create the analysis_var
        self.analysis_var = tk.StringVar(value="taxonomy-wgs")  # Default to taxonomy-wgs

        
        
        

    def select_and_monitor_folders(self):
        folder = filedialog.askdirectory(title="Select Folder to Monitor")
        if folder:
            self.start_folder_monitoring([folder])  # Pass as a list

    def start_folder_monitoring(self, folders):
        if not isinstance(folders, list):
            folders = [folders]
        self.folder_monitor = FolderMonitor(folders, self)
        self.folder_monitor.start()
        self.show_file_list_window()


    def init_layout(self):
        ctk.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"

        # PRIMARY_COLOR = "#2C3E50"  # Dark blue for primary elements and backgrounds
        # SECONDARY_COLOR = "#EAECEE"  # Light gray for secondary backgrounds
        # TEXT_COLOR = "#FFFFFF"  # White text for better contrast on dark backgrounds
        # BUTTON_COLOR = "#ff8906"  # Bright blue for buttons
        # BUTTON_TEXT_COLOR = "#FFFFFF"

        def resize_icon_to_xpx(image_path, x):
            """Utility function to resize the image to 25px while maintaining aspect ratio."""
            icon = Image.open(image_path)
            base_width = x
            w_percent = base_width / float(icon.width)
            h_size = int(float(icon.height) * float(w_percent))

            icon = icon.resize((base_width, h_size), Image.LANCZOS)
            # Convert the PIL image to a CTkImage
            ctk_image = CTkImage(icon, size=(x, x))

            return ctk_image

        self.geometry(f"{MAIN_WINDOW_X}x{MAIN_WINDOW_Y}")
        self.title(f"MMonitor {VERSION}")
        self.minsize(MAIN_WINDOW_X, MAIN_WINDOW_Y)
        self.configure()  # Set primary color as window background

        # Style and theme

        # style.map('TButton',
        #           foreground=[('pressed', 'white'), ('active', 'black')],
        #           background=[('pressed', '!disabled', '#B22222'), ('active', '#B0E0E6')]  # Light blue when active
        #           )

        # Placeholder icon
        # placeholder_icon = tk.PhotoImage(width=16, height=)
        btn_icon_size = 35
        create_local_db_icon = resize_icon_to_xpx(f'{IMAGES_PATH}/mmonitor_button_add_db.png', btn_icon_size)
        add_db_icon = resize_icon_to_xpx(f'{IMAGES_PATH}/mmonitor_button_import_db.png', btn_icon_size)
        start_offline_monitoring_icon = resize_icon_to_xpx(f'{IMAGES_PATH}/offline_monitoring.png', btn_icon_size)
        user_authentication_icon = resize_icon_to_xpx(f"{IMAGES_PATH}/mmonitor_button4_authenticate.png", btn_icon_size)
        add_metadata_icon = resize_icon_to_xpx(f"{IMAGES_PATH}/mmonitor_button_importcsv.png", btn_icon_size)
        run_analysis_pipeline_icon = resize_icon_to_xpx(f"{IMAGES_PATH}/button_add_data2.png", btn_icon_size)
        quit_icon = resize_icon_to_xpx(f"{IMAGES_PATH}/mmonitor_button_quit.png", btn_icon_size)

        # placeholder_icon.put(("gray",), to=(0, 0, 15, 15))

        # Header with app title

        header_label = ctk.CTkLabel(self, text=f"Metagenome Monitor {VERSION}", font=("Helvetica", 18))
        # header_label.bg_color = PRIMARY_COLOR
        # header_label.fg_color = TEXT_COLOR
        header_label.pack(pady=10)

        # Categories and buttons

        # Define the 'categories' list
        categories = [
            # ("Local", [
                # ("Create Local DB", self.create_project, create_local_db_icon),
                # ("Choose Local DB", self.choose_project, add_db_icon),
                # ("Start Local dashboard", self.start_monitoring, start_offline_monitoring_icon)

            # ]),
            ("Webserver", [
                ("Select User", self.open_db_config_form, user_authentication_icon)

            ]),
            ("Add Data", [
                ("Run Analysis Pipeline", self.checkbox_popup, run_analysis_pipeline_icon)
            ])
        ]
        # Tooltips for the categories
        category_tooltips = {
            "Local": (
                "Local Database Options:\n"
                "- Create a new local database.\n"
                "- Choose an existing local database.\n"
                "- Begin offline monitoring.\n"
            ),
            "Webserver": (
                "Webserver Authentication:\n"
                "- Provide your username and password for the webserver.\n"
                "- Ensure you use the same credentials as on the MMonitor webpage.\n"
                "- Do not modify 'host' field (default 134.2.78.150).\n"
            ),
            "Add Data": (
                "Data Addition Options:\n"
                "- Append metadata using a CSV (e.g., 'meta.csv').\n"
                "- Process sequencing files for analysis.\n"
                "- Ensure either a local database is chosen or you're authenticated for data addition.\n"
            )
        }

        # Calculate the maximum button width based on text length
        button_texts = [btn[0] for cat in categories for btn in cat[1]]
        button_texts.append("Quit")  # Adding Quit button text
        max_text_length = max(map(len, button_texts))
        btn_width = max_text_length  # Adding an offset to account for padding and icon
        style = ttk.Style()

        # Styling improvements for buttons

        # style.theme_use("clam")

        # Configuring Button Style

        for category, btns in categories:
            cat_label = ctk.CTkLabel(self, text=category, font=("Helvetica", 20), anchor="center")
            # cat_label.bg_color = PRIMARY_COLOR
            # cat_label.fg_color = TEXT_COLOR
            cat_label.pack(pady=10)

            for text, cmd, img in btns:
                btn = ctk.CTkButton(self, text=text, command=cmd, image=img, width=210, height=40
                                    )

                # btn.bg_color = BUTTON_COLOR
                # btn.fg_color = BUTTON_TEXT_COLOR
                # btn.hover_bg_color = BUTTON_COLOR  # Set the hover background color
                # btn.hover_fg_color = TEXT_COLOR  # Set the hover foreground color
                # btn.active_bg_color = BUTTON_COLOR  # Set the active background color
                # btn.active_fg_color = TEXT_COLOR  # Set the active foreground color
                btn.pack(pady=2)

        # Quit button
        quit_btn = ctk.CTkButton(self, width=210, height=40, text="Quit", command=self.stop_app, image=quit_icon)
        # quit_btn.bg_color = BUTTON_COLOR
        # quit_btn.fg_color = BUTTON_TEXT_COLOR
        # quit_btn.hover_bg_color = BUTTON_COLOR  # Set the hover background color
        # quit_btn.hover_fg_color = TEXT_COLOR  # Set the hover foreground color
        # quit_btn.active_bg_color = BUTTON_COLOR  # Set the active background color
        # quit_btn.active_fg_color = TEXT_COLOR  # Set the active foreground color
        quit_btn.pack(pady=15)

    # create_tooltip(local_label, "This is the tooltip text for the Local category.")
    # create_tooltip(webserver_label, "This is the tooltip text for the Webserver category.")

    def create_tooltip(self, widget, text):
        tooltip = ToolTip(widget, text)
        widget.bind("<Enter>", tooltip.show_tip)
        widget.bind("<Leave>", tooltip.hide_tip)

    def open_calendar(self):
        calendar_window = tk.Toplevel()
        calendar_window.title("Choose a Date")

        # Use a StringVar to hold our date
        selected_date_var = tk.StringVar()

        calendar = Calendar(calendar_window, selectmode="day")
        calendar.pack(pady=10)

        def on_date_select():
            selected_date_var.set(calendar.get_date())
            calendar_window.destroy()
            calendar_window.quit()

        select_button = ttk.Button(calendar_window, text="Select Date", command=on_date_select, style="TButton")
        select_button.pack(pady=20)

        # calendar_window.mainloop()

        return selected_date_var.get()

    def open_popup(self, text, title, icon):

        CTkMessagebox(message=text, title=title, icon=icon, option_1="Okay")

    def create_project(self):
        filename = filedialog.asksaveasfilename(
            initialdir='projects/',
            title="Choose place to safe the project data"
        )
        filename += ".sqlite3"
        self.db_path = filename
        self.db = MMonitorDBInterface(filename)
        self.db.create_db(filename)

    def choose_project(self):

        self.db_path = filedialog.askopenfilename(
            initialdir='projects/',
            title="Choose project data base to use",
            filetypes=(("sqlite", "*.sqlite3"), ("all files", "*.*"))
        )
        self.db = MMonitorDBInterface(self.db_path)
        # self.db = "mmonitor.sqlite3"

    # def choose_index(self):
    #     self.centrifuge_index = filedialog.askopenfilename(
    #         initialdir='projects/',
    #         title="Choose project data base to use",
    #         filetypes=(("sqlite", "*.sqlite3"), ("all files", "*.*"))
    #     )

    # @require_project
    def append_metadata(self):
        csv_file = filedialog.askopenfilename(
            initialdir='projects/',
            title="Choose csv file containing metadata to append",
            filetypes=(("csv", "*.csv"), ("all files", "*.*"))
        )
        if csv_file is not None and len(csv_file) > 0:
            self.django_db.append_metadata_from_csv(csv_file)

        # ceck if a file exists and if not asks the user to download it. gets used to check if db are all present
        # TODO: also add checksum check to make sure the index is completely downloaded, if not remove file and download again

    def check_emu_db_exists(self):
        # unzip if tar exists, but not taxonomy.tsv
        if os.path.exists(f"{ROOT}/src/resources/emu_db/emu.py.tar") and not os.path.exists(
                f"{ROOT}/src/resources/emu_db/taxonomy.tsv"):
            self.open_popup("Unzipping tar", "Unzipping emu.py.tar", "check")
            self.unzip_tar(f"{ROOT}/src/resources/emu_db/emu.py.tar", f"{ROOT}/src/resources/emu_db/")

        if not os.path.exists(f"{ROOT}/src/resources/emu_db/taxonomy.tsv"):
            response = CTkMessagebox("Emu DB not found. Do you want to download it?", option_1="Yes", option_2="No")
            if response.selected_option == "Yes":
                try:
                    # download emu.py db from web if it doesn't exist
                    emu_db_path = f"{ROOT}/src/resources/emu_db/emu.tar"
                    self.download_file_from_web(emu_db_path,
                                                "https://software-ab.cs.uni-tuebingen.de/download/MMonitor/emu.tar")
                    computed_checksum = compute_sha256(emu_db_path)
                    print(f"Checksum of the downloaded file: {computed_checksum}")

                    self.unzip_tar(f"{ROOT}/src/resources/emu_db/emu.tar", f"{ROOT}/src/resources/emu_db/")

                except FileNotFoundError as e:
                    self.open_popup("Could not download the EMU DB. Please contact the MMonitor developer for help",
                                    "Could not find emu.py db", icon="cancel")

    def download_file_from_web(self, filepath, url):

        with urllib.request.urlopen(url) as response:
            # get file size from content-length header
            file_size = int(response.info().get("Content-Length"))
            # create progress bar widget
            progress = ttk.Progressbar(self, orient="horizontal", length=250, mode="determinate")
            progress.pack()
            self.progress_bar_exists = True
            progress["maximum"] = file_size
            progress["value"] = 0

            # Start the download in a separate thread
            download_thread = Thread(target=self._download_file, args=(filepath, url, progress))
            download_thread.start()

            # Periodically update the GUI
            self.after(100, self._check_download_progress, download_thread, progress)
            download_thread.join()

    def _download_file(self, filepath, url, progress):
        def update_progress(count, block_size, total_size):
            self.download_progress = count * block_size

        # Reset progress for a new download
        self.download_progress = 0

        # Download file and update progress bar
        urllib.request.urlretrieve(url, filepath, reporthook=update_progress)

        # Once download is complete
        self.progress_bar_exists = False
        progress.destroy()
        CTkMessagebox(message="Download complete. Unpacking files...", icon="check")

    def _check_download_progress(self, download_thread, progress):
        # Update the progress bar with the latest progress value

        if self.progress_bar_exists:
            progress["value"] = self.download_progress
            self.update_idletasks()

        # If the download thread is still running, keep checking
        if download_thread.is_alive():
            self.after(100, self._check_download_progress, download_thread, progress)

    #
    # def handle_kaiju_output(self):
    #     kaiju_runner = KaijuRunner()
    #
    #     # Using SampleInput window to get user input
    #     sample_input_window = SampleInputWindow()  # Assuming you have a class named SampleInputWindow
    #     sequence_list, sample_name = sample_input_window.get_user_input()  # Assuming this method returns a list of sequence files and a sample name
    #
    #     kaiju_runner.run_kaiju(sequence_list, sample_name)
    #
    #     # Now, we'll add the output to the DjangoDB
    #     db = DjangoDBInterface()
    #     db.add_kaiju_output_to_db(sample_name)
    #
    #
    #     db = DjangoDBInterface()  # Assuming this class has methods to handle Kaiju output
    #     db.add_kaiju_output_to_db(sample_name)

    def check_file_exists(self, filepath, url):
        if os.path.exists(f"{ROOT}/src/resources/dec22.tar"):
            if not os.path.exists(f"{ROOT}/src/resources/dec22.1.cf"):
                response = CTkMessagebox(
                    "Index is compressed. Do you want to decompress it?", option_1="Yes", option_2="No")
                if response.selected_option == "Yes":
                    try:
                        self.unzip_tar(f"{ROOT}/src/resources/dec22.tar", f"{ROOT}/src/resources/")
                    except FileNotFoundError as e:
                        print(f"Requested file not found")
            if not os.path.exists(f"{ROOT}/src/resources/dec22.1.cf") and os.path.exists(
                    f"{ROOT}/src/resources/dec22.1.cf.gz"):
                try:
                    self.unzip_gz(f"{ROOT}/src/resources/dec22.1.cf.gz")
                    self.unzip_gz(f"{ROOT}/src/resources/dec22.2.cf.gz")
                    self.unzip_gz(f"{ROOT}/src/resources/dec22.3.cf.gz")
                    self.unzip_gz(f"{ROOT}/src/resources/dec22.4.cf.gz")
                except FileNotFoundError as e:
                    print(f"Requested files not found.")





        else:
            response = CTkMessagebox("Centrifuge index not found. Do you want to download it?"
                                     " Might take some time, the tool is unusable while download.", option_1="Yes",
                                     option_2="No")
            if response.selected_option == "Yes":
                self.download_file_from_web(filepath, url)

    def unzip_tar(self, file, out_folder):
        my_tar = tarfile.open(file, mode='r')
        my_tar.extractall(out_folder)  # specify which folder to extract to
        my_tar.close()

    def unzip_gz(self, file):
        with gzip.open(file, 'rb') as f:
            file_content = f.read()
            gzip.decompress(file_content)

    # choose folder containing sequencing data
    # TODO: check if there is white space in path that causes problem
    # @require_project
    # @require_centrifuge
    def taxonomy_nanopore_wgs(self):
        def add_sample_to_databases(sample_name, project_name, subproject_name, sample_date):
            kraken_out = f"{ROOT}/src/resources/pipeline_out/{sample_name}_kraken_out"

            if self.db is not None:
                self.db.update_table_with_kraken_out(kraken_out, "species", sample_name, "project", self.sample_date)
            self.django_db.send_nanopore_record_centrifuge(kraken_out, sample_name, project_name, subproject_name,
                                                           sample_date, True)

        # check if centrifuge index exists, if not download it using check_file_exists method
        centrifuge_index_path = f"{ROOT}/src/resources/dec_22"
        self.centrifuge_index_path = centrifuge_index_path
        # download_thread = Thread(target=self.check_file_exists(f"{ROOT}/src/resources/dec22.tar",
        #                                                        "https://software-ab.cs.uni-tuebingen.de/download/MMonitor/dec22.tar"))
        # download_thread.start()
        # download_thread.join()

        # self.unzip_tar(f"{ROOT}/src/resources/dec22.tar", f"{ROOT}/src/resources/")
        # self.unzip_gz(f"{ROOT}/src/resources/dec22.1.cf.gz")
        # self.unzip_gz(f"{ROOT}/src/resources/dec22.2.cf.gz")
        # self.unzip_gz(f"{ROOT}/src/resources/dec22.3.cf.gz")
        # self.unzip_gz(f"{ROOT}/src/resources/dec22.4.cf.gz")

        # folder = filedialog.askdirectory(
        #     initialdir='/',
        #     title="Choose directory containing sequencing data"
        # )
        # files = self.centrifuge_runner.get_files_from_folder(folder)
        #
        # sample_name = simpledialog.askstring(
        #     "Input sample name",
        #     "What should the sample be called?",
        #     parent=self
        # )
        #
        #
        # sample_date = self.open_calendar()
        self.open_input_window_and_wait()
        if self.input_window.do_quit:
            return

        print("Created input window")
        # self.checkbox_popup()
        # get entries from input window
        # sample_name = str(self.input_window.sample_name_entry)  # Get the content of the entry and convert to string
        # project_name = str(self.input_window.project_name_entry)
        # subproject_name = str(self.input_window.subproject_name_entry)
        # sample_date = self.input_window.selected_date.strftime('%Y-%m-%d')  # Convert date to string format
        # files = self.input_window.file_paths_single_sample

        if self.input_window.process_multiple_samples:
            for index, file_path_list in enumerate(self.input_window.multi_sample_input["file_paths_lists"]):
                files = file_path_list
                sample_name = self.input_window.multi_sample_input["sample_names"][index]
                project_name = self.input_window.multi_sample_input["project_names"][index]
                subproject_name = self.input_window.multi_sample_input["subproject_names"][index]
                sample_date = self.input_window.multi_sample_input["dates"][index]
                self.centrifuge_runner.run_centrifuge(files, sample_name, centrifuge_index_path)
                self.centrifuge_runner.make_kraken_report(centrifuge_index_path)
                

                self.add_statistics(self.centrifuge_runner.concat_file_name, sample_name, project_name, subproject_name,
                                    sample_date)

                add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
            self.show_info("Analysis complete. You can start monitoring now.")

        # self.db.update_table_with_kraken_out(f"{ROOT}/src/resources/pipeline_out/{sample_name}_kraken_out", "species",
        #                                      sample_name, project_name, sample_date)

    def add_statistics(self, fastq_file, sample_name, project_name, subproject_name, sample_date):
        fastq_stats = FastqStatistics(fastq_file)

        # Calculate statistics
        fastq_stats.quality_statistics()
        fastq_stats.read_lengths_statistics()
        quality_vs_lengths_data = fastq_stats.qualities_vs_lengths()
        gc_contents = fastq_stats.gc_content_per_sequence()

        data = {
            'sample_name': sample_name,
            'project_id': project_name,
            'subproject_id': subproject_name,
            'date': sample_date,
            'mean_gc_content': fastq_stats.gc_content(),
            'mean_read_length': np.mean(fastq_stats.lengths),
            'median_read_length': np.median(fastq_stats.lengths),
            'mean_quality_score': np.mean([np.mean(q) for q in fastq_stats.qualities]),
            'read_lengths': json.dumps(quality_vs_lengths_data['read_lengths']),
            'avg_qualities': json.dumps(quality_vs_lengths_data['avg_qualities']),
            'number_of_reads': fastq_stats.number_of_reads(),
            'total_bases_sequenced': fastq_stats.total_bases_sequenced(),
            'q20_score': fastq_stats.q20_q30_scores()[0],
            'q30_score': fastq_stats.q20_q30_scores()[1],
            # 'avg_quality_per_read': fastq_stats.quality_score_distribution()[0],
            # 'base_quality_avg': fastq_stats.quality_score_distribution()[1],
            'gc_contents_per_sequence': json.dumps(gc_contents)

        }

        self.django_db.send_sequencing_statistics(data)

    def get_metadata_from_input_window(self):
        if not self.input_window.process_multiple_samples:
            sample_name = str(self.input_window.sample_name)  # Get the content of the entry and convert to string
            project_name = str(self.input_window.project_name)
            subproject_name = str(self.input_window.subproject_name)
            try:
                sample_date = self.input_window.selected_date.strftime('%Y-%m-%d')  # Convert date to string format
            except AttributeError as e:
                sample_date = datetime.today()
                self.open_popup(f"AttributeError. Please fill out all input fields or use CSV for sample input.",
                                f"AttributeError", icon="cancel")
                print(e)
                return
            files = self.input_window.file_paths_single_sample
            return sample_name, project_name, subproject_name, sample_date, files
        else:
            files = self.input_window.multi_sample_input["file_paths_lists"]
            sample_names = self.input_window.multi_sample_input["sample_names"]
            project_names = self.input_window.multi_sample_input["project_names"]
            subproject_names = self.input_window.multi_sample_input["subproject_names"]
            sample_dates = self.input_window.multi_sample_input["dates"]
            return sample_names, project_names, subproject_names, sample_dates, files




    def taxonomy_nanopore_16s(self):
        global sample_name

        def add_sample_to_databases(sample_name, project_name, subproject_name, sample_date):
            emu_out_path = f"{ROOT}/src/resources/pipeline_out/{sample_name}/"
            if self.db is not None:
                self.db.update_table_with_emu_out(emu_out_path, "species", sample_name, "project", self.sample_date)

            self.django_db.update_django_with_emu_out(emu_out_path, "species", sample_name, project_name, sample_date,
                                                      subproject_name, True)

        self.check_emu_db_exists()
        # create input window to input all relevant sample information and sequencing files

        self.open_input_window_and_wait()

        # quit the method when quit button is pressed instead of running the pipeline
        if self.input_window.do_quit:
            return
        # get entries from input window
        if not self.input_window.process_multiple_samples:
            sample_name = str(self.input_window.sample_name)  # Get the content of the entry and convert to string
            project_name = str(self.input_window.project_name)
            subproject_name = str(self.input_window.subproject_name)
            try:
                sample_date = self.input_window.selected_date.strftime('%Y-%m-%d')  # Convert date to string format
            except AttributeError as e:
                sample_date = datetime.today()
                self.open_popup(f"AttributeError. Please fill out all input fields or use CSV for sample input.",
                                f"AttributeError", icon="cancel")
                print(e)
                return
            files = self.input_window.file_paths_single_sample
            self.emu_runner.run_emu(files, sample_name, 0.01)
            print("add statistics")
            self.add_statistics(self.emu_runner.concat_file_name, sample_name, project_name, subproject_name,
                                sample_date)

            add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)
        else:
            print("Processing multiple samples")
            for index, file_path_list in enumerate(self.input_window.multi_sample_input["file_paths_lists"]):
                files = file_path_list
                sample_name = self.input_window.multi_sample_input["sample_names"][index]
                project_name = self.input_window.multi_sample_input["project_names"][index]
                subproject_name = self.input_window.multi_sample_input["subproject_names"][index]
                sample_date = self.input_window.multi_sample_input["dates"][index]
                self.emu_runner.run_emu(files, sample_name, 0.01)
                self.add_statistics(self.emu_runner.concat_file_name, sample_name, project_name, subproject_name,
                                    sample_date)

                add_sample_to_databases(sample_name, project_name, subproject_name, sample_date)

        # emu_out_path = f"{ROOT}/src/resources/pipeline_out/subset/"

        # if self.db is not None:
        #     self.db.update_table_with_emu_out(emu_out_path, "species", sample_name, "project", self.sample_date)
        #
        # self.db_mysql.update_django_with_emu_out(emu_out_path, "species", sample_name, project_name, sample_date,
        #                                          subproject_name)

        self.show_info("Analysis complete. You can start monitoring now.")

    def functional_pipeline(self):
        self.open_input_window_and_wait()
        if self.input_window.do_quit:
            return
        if not self.input_window.process_multiple_samples:
            sample_name, project_name, subproject_name, sample_date, files = self.get_metadata_from_input_window()
            self.cmd_runner.assembly_pipeline(sample_name, project_name, subproject_name, sample_date, files)
        else:
            sample_names, project_names, subproject_names, sample_dates, files = self.get_metadata_from_input_window()
            for idx, sample in enumerate(sample_names):
                self.cmd_runner.assembly_pipeline(sample_names[idx], project_names[idx], subproject_names[idx],
                                                  sample_dates[idx], files[idx])



    def open_input_window_and_wait(self):
        self.input_window = InputWindow(self, self.emu_runner)
        time.sleep(1)
        print("Before wait_window")
        self.wait_window(self.input_window)
        print("After wait_window")

    def checkbox_popup(self):
        self.pipeline_popup = PipelinePopup(self,
                                            self)  # Replace run_analysis_pipeline_function with your actual function

    def ask_sample_name(self):
        sample_name = simpledialog.askstring(
            "Input sample name",
            "What should the sample be called?",
            parent=self
        )
        return sample_name

    # @require_project

    def show_info(self, message):
        CTkMessagebox(message=message, icon="check", title="Info")

    @require_project
    def start_monitoring(self):
        try:
            self.dashapp = Index(self.db)
            self.monitor_thread = Thread(target=self.dashapp.run_server, args=(False,))
            self.monitor_thread.start()
        except IndexError:
            self.show_info(
                "No data found in database. Please first run analysis pipeline to fill DB with data.")
            return

        sleep(1)
        open_new('http://localhost:8050')

    def start_app(self):
        self.mainloop()

    def stop_app(self):
        if self.monitor_thread is not None and self.monitor_thread.is_alive():
            post('http://localhost:8050/shutdown')
            self.monitor_thread.join()
        if self.folder_monitor:
            self.folder_monitor.stop()
        self.destroy()


    def open_db_config_form(self):
        db_config_form = DataBaseConfigForm(master=self)
        print(db_config_form.last_config)

    def check_and_start_analysis(self, sample_name):
        # Example condition: if the sample has accumulated 10 files
        if len(self.sample_files[sample_name]) >= 10:
            print(f"Starting analysis for sample {sample_name}")
            self.run_analysis_for_sample(sample_name)

    def run_analysis_for_sample(self, sample_name):
        files = self.sample_files[sample_name]
        analysis_type = self.get_selected_analysis_type()
        if analysis_type == "taxonomy-16s":
            self.emu_runner.run_emu(files, sample_name, min_abundance=0.01)
        elif analysis_type == "taxonomy-wgs":
            self.centrifuge_runner.run_centrifuge(files, sample_name, self.centrifuge_index_path)
        self.sample_files[sample_name] = []
        self.refresh_file_tree()

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

    def run_pipeline(self, analysis_type):
        input_window = InputWindow(self, self.emu_runner)
        self.wait_window(input_window)

        if input_window.do_quit:
            return

        if input_window.process_multiple_samples:
            self.run_multi_sample_pipeline(analysis_type, input_window.multi_sample_input)
        else:
            self.run_single_sample_pipeline(analysis_type, input_window)

    def run_single_sample_pipeline(self, analysis_type, input_window):
        sample_name = input_window.sample_name
        project_name = input_window.project_name
        subproject_name = input_window.subproject_name
        sample_date = input_window.selected_date.strftime('%Y-%m-%d')
        files = input_window.file_paths_single_sample

        # Get the directory containing the mmonitor package
        mmonitor_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        
        # Use a default config file if self.db_path is None
        config_path = self.db_path if self.db_path else os.path.join(ROOT, "src", "resources", "db_config.json")
        
        # Construct the command
        cmd = [
            sys.executable,  # Use the current Python interpreter
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", analysis_type,
            "-c", config_path,
            "-i", *files,
            "-s", sample_name,
            "-d", sample_date,
            "-p", project_name,
            "-u", subproject_name
        ]
        
        # Set the PYTHONPATH to include the mmonitor directory
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{mmonitor_dir}:{env.get('PYTHONPATH', '')}"

        try:
            subprocess.Popen(cmd, env=env)
        except Exception as e:
            print(f"Error running pipeline: {e}")
            # You might want to show an error message to the user here

    def run_multi_sample_pipeline(self, analysis_type, multi_sample_input):
        # Create a temporary CSV file with the multi-sample input data
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

        # Get the directory containing the mmonitor package
        mmonitor_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

        # Use a default config file if self.db_path is None
        config_path = self.db_path if self.db_path else os.path.join(ROOT, "src", "resources", "db_config.json")

        # Construct the command
        cmd = [
            sys.executable,  # Use the current Python interpreter
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", analysis_type,
            "-c", config_path,
            "-m", temp_csv_path
        ]

        # Set the PYTHONPATH to include the mmonitor directory
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{mmonitor_dir}:{env.get('PYTHONPATH', '')}"

        try:
            subprocess.Popen(cmd, env=env)
        except Exception as e:
            print(f"Error running pipeline: {e}")
            # You might want to show an error message to the user here

    
class ToolTip:
    def __init__(self, widget, tip_text):
        self.widget = widget
        self.tip_text = tip_text
        self.tip_window = None

    def show_tip(self, event=None):
        x, y, _, _ = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 25

        self.tip_window = tk.Toplevel(self.widget)
        self.tip_window.wm_overrideredirect(True)
        self.tip_window.wm_geometry(f"+{x}+{y}")

        label = tk.Label(self.tip_window, text=self.tip_text, foreground="black", background="white", relief="solid", borderwidth=1,
                         font=("Helvetica", "12", "normal"))
        label.pack(ipadx=1)

    def hide_tip(self, event=None):
        if self.tip_window:
            self.tip_window.destroy()
            self.tip_window = None

