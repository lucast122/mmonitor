import tkinter as tk
from threading import Thread, Lock
import customtkinter as ctk
import os
import json
from tkinter import filedialog
from mmonitor.userside.MMonitorCMD import MMonitorCMD

from build_mmonitor_pyinstaller import ROOT

class PipelinePopup(ctk.CTkFrame):
    def __init__(self, parent, gui_ref):
        super().__init__(parent)
        ctk.set_default_color_theme("blue")
        self.parent = parent
        self.gui = gui_ref
        self.taxonomy_nanopore_wgs = tk.BooleanVar()
        self.taxonomy_nanopore_16s_bool = tk.BooleanVar()
        self.assembly = tk.BooleanVar()
        self.correction = tk.BooleanVar()
        self.binning = tk.BooleanVar()
        self.annotation = tk.BooleanVar()
        self.kegg = tk.BooleanVar()
        self.test_mag_upload = tk.BooleanVar()
        # default databases paths to be used if user doesn't supply own path
        self.emu_db_path = os.path.join(ROOT, "src", "resources", "emu_db/")
        self.centrifuge_db_path = os.path.join(ROOT, "src", "resources", "centrifuge_db/")
        
        self.create_widgets()
        
    def create_widgets(self):
        main_frame = ctk.CTkFrame(self)
        main_frame.pack(padx=20, pady=20, fill="both", expand=True)
        
        # Title and explanation
        title_label = ctk.CTkLabel(main_frame, text="Analysis Pipeline", font=("Helvetica", 24, "bold"))
        title_label.pack(pady=(0, 10))
        
        explanation = ("Run one or multiple analysis pipelines on existing files.\n"
                       "You can select which analysis to run and which databases to use."
                       "For more database options, go to Manage Databases."
                       "Please select either 16S-rRNA or whole genome nanopore analysis for taxonomy and"
                       "a functional analysis pipeline if whole genome sequencing was performed.")
        ctk.CTkLabel(main_frame, text=explanation, wraplength=500).pack(pady=(0, 20))
        
        # Taxonomy Analysis Section
        tax_frame = self.create_section(main_frame, "Taxonomic Analysis")
        ctk.CTkCheckBox(tax_frame, text='Taxonomy whole genome nanopore', variable=self.taxonomy_nanopore_wgs).pack(pady=5, anchor="w")
        ctk.CTkCheckBox(tax_frame, text='Taxonomy 16S-rRNA reads nanopore', variable=self.taxonomy_nanopore_16s_bool).pack(pady=5, anchor="w")
        
        # Functional Analysis Section
        func_frame = self.create_section(main_frame, "Functional Analysis")
        ctk.CTkCheckBox(func_frame, text='Assembly Pipeline', variable=self.assembly).pack(pady=5, anchor="w")
        ctk.CTkCheckBox(func_frame, text='KEGG', variable=self.kegg).pack(pady=5, anchor="w")
        
        # Database Selection Section
        db_frame = self.create_section(main_frame, "Database Selection")
        ctk.CTkButton(db_frame, text="Select Emu Database", command=self.select_emu_db).pack(pady=5, fill="x")
        ctk.CTkButton(db_frame, text="Select Centrifuge Database", command=self.select_centrifuge_db).pack(pady=5, fill="x")
        
        # Action Buttons
        button_frame = ctk.CTkFrame(main_frame)
        button_frame.pack(pady=20, fill="x")
        ctk.CTkButton(button_frame, text="Continue", command=self.run_analysis_pipeline).pack(side="left", padx=5, expand=True, fill="x")
        ctk.CTkButton(button_frame, text="Quit", command=self.destroy).pack(side="right", padx=5, expand=True, fill="x")

    def create_section(self, parent, title):
        frame = ctk.CTkFrame(parent)
        frame.pack(pady=10, fill="x")
        ctk.CTkLabel(frame, text=title, font=("Helvetica", 16, "bold")).pack(pady=5)
        return frame

    def select_emu_db(self):
        db_path = filedialog.askdirectory(title="Select Emu Database Directory")
        if db_path:
            self.emu_db_path = db_path

    def select_centrifuge_db(self):
        db_path = filedialog.askdirectory(title="Select Centrifuge Database Directory")
        if db_path:
            self.centrifuge_db_path = db_path

    def run_analysis_pipeline(self):
        if self.assembly.get():
            analysis_type = "assembly"
            cmd_runner = MMonitorCMD()
            args = cmd_runner.parse_arguments([
                "-a", analysis_type,
                "-c", self.gui.db_path,  # Use gui_ref instead of parent
                "-i"] + self.gui.sample_files.get(self.gui.sample_name, []) + [
                "-s", self.gui.sample_name,
                "-d", self.gui.selected_date.strftime("%Y-%m-%d") if hasattr(self.gui, 'selected_date') else datetime.now().strftime("%Y-%m-%d"),
                "-p", self.gui.project_name if hasattr(self.gui, 'project_name') else "",
                "-u", self.gui.subproject_name if hasattr(self.gui, 'subproject_name') else "",
                "--overwrite" if hasattr(self.gui, 'overwrite') and self.gui.overwrite else ""
            ])
            cmd_runner.initialize_from_args(args)
            cmd_runner.run()
        else:
            # Handle other analysis types as before
            pass