import tkinter as tk
import customtkinter as ctk
import os
import json
from tkinter import filedialog, messagebox
import multiprocessing
from mmonitor.userside.utils import create_tooltip
from mmonitor.userside.MMonitorCMD import MMonitorCMD
from build_mmonitor_pyinstaller import ROOT

class PipelinePopup(ctk.CTkFrame):
    def __init__(self, parent, gui_ref):
        super().__init__(parent)
        ctk.set_default_color_theme("blue")
        self.parent = parent
        self.gui = gui_ref
        
        # Initialize variables
        self.analysis_type = tk.StringVar(value="taxonomy-wgs")
        self.threads = tk.StringVar(value=str(multiprocessing.cpu_count()))
        self.min_length = tk.StringVar(value="1000")
        self.min_quality = tk.StringVar(value="10")
        self.emu_db = tk.StringVar(value="")
        self.centrifuge_db = tk.StringVar(value="")
        self.min_abundance = tk.StringVar(value="0.01")
        self.assembly_mode = tk.StringVar(value="nano-raw")
        self.medaka_model = tk.StringVar(value="r1041_e82_400bps_sup_v5.0.0")
        self.is_isolate = tk.BooleanVar(value=False)
        self.min_contig_length = tk.StringVar(value="1000")
        
        # Default database paths
        self.emu_db_path = os.path.join(ROOT, "src", "resources", "emu_db")
        self.centrifuge_db_path = os.path.join(ROOT, "src", "resources", "centrifuge_db")
        
        self.config_file = os.path.join(os.path.expanduser("~"), ".mmonitor_config.json")
        self.load_config()
        self.create_widgets()
        self.search_available_databases()

    def create_widgets(self):
        main_frame = ctk.CTkFrame(self)
        main_frame.pack(padx=10, pady=10, fill="both", expand=True)

        # Title and explanation
        title_label = ctk.CTkLabel(main_frame, text="Analysis Pipeline Configuration", font=("Helvetica", 20, "bold"))
        title_label.pack(pady=(0, 5))
        
        explanation = ("Configure the parameters for your analysis pipeline.\n"
                       "Set general parameters, taxonomy-specific settings, and assembly options.")
        ctk.CTkLabel(main_frame, text=explanation, wraplength=400).pack(pady=(0, 10))

        # Scrollable frame for content
        content_frame = ctk.CTkScrollableFrame(main_frame)
        content_frame.pack(fill="both", expand=True, padx=5, pady=5)

        # General parameters
        self.create_section(content_frame, "General Parameters", [
            ("Analysis Type:", self.analysis_type, ["taxonomy-wgs", "taxonomy-16s", "assembly", "functional"], "option"),
            ("CPU cores:", self.threads),
            ("Min read length:", self.min_length),
            ("Min read quality:", self.min_quality)
        ])

        # Taxonomy parameters
        self.create_section(content_frame, "Taxonomy Parameters", [
            ("Emu Database:", self.emu_db),
            ("Centrifuge Database:", self.centrifuge_db),
            ("Min abundance:", self.min_abundance)
        ])

        # Assembly parameters
        self.create_section(content_frame, "Assembly Parameters", [
            ("Assembly mode:", self.assembly_mode, ["nano-raw", "nano-corrected"], "option"),
            ("Medaka model:", self.medaka_model, self.get_sorted_medaka_models(), "option"),
            ("Is isolate:", self.is_isolate, None, "checkbox"),
            ("Min contig length:", self.min_contig_length)
        ])

        # Save Parameters button
        save_button = ctk.CTkButton(main_frame, text="Save Parameters", command=self.save_parameters)
        save_button.pack(pady=10)
        create_tooltip(save_button, "Save the current parameters to a configuration file")

    def create_section(self, parent, title, fields):
        section_frame = ctk.CTkFrame(parent)
        section_frame.pack(pady=5, fill="x")
        ctk.CTkLabel(section_frame, text=title, font=("Helvetica", 14, "bold")).pack(pady=2)

        for field in fields:
            field_frame = ctk.CTkFrame(section_frame)
            field_frame.pack(fill="x", pady=1)
            label = ctk.CTkLabel(field_frame, text=field[0], width=100, anchor="w")
            label.pack(side="left", padx=2)
            
            if len(field) > 3 and field[3] == "option":
                widget = ctk.CTkOptionMenu(field_frame, variable=field[1], values=field[2])
            elif len(field) > 3 and field[3] == "checkbox":
                widget = ctk.CTkCheckBox(field_frame, text="", variable=field[1])
            else:
                widget = ctk.CTkEntry(field_frame, textvariable=field[1])
            
            widget.pack(side="left", expand=True, fill="x", padx=2)
            create_tooltip(widget, f"Set the {field[0].lower()[:-1]}")

            # Store references to the database entry widgets
            if field[0] == "Emu Database:":
                self.emu_db_entry = widget
            elif field[0] == "Centrifuge Database:":
                self.centrifuge_db_entry = widget

    def load_config(self):
        if os.path.exists(self.config_file):
            with open(self.config_file, 'r') as f:
                self.config = json.load(f)
            self.update_widgets_from_config()
        else:
            self.config = {}

    def update_widgets_from_config(self):
        for key, value in self.config.items():
            if hasattr(self, key):
                getattr(self, key).set(value)

    def save_parameters(self):
        config = {
            'analysis_type': self.analysis_type.get(),
            'threads': self.threads.get(),
            'min_length': self.min_length.get(),
            'min_quality': self.min_quality.get(),
            'emu_db': self.emu_db.get(),
            'centrifuge_db': self.centrifuge_db.get(),
            'min_abundance': self.min_abundance.get(),
            'assembly_mode': self.assembly_mode.get(),
            'medaka_model': self.medaka_model.get(),
            'is_isolate': self.is_isolate.get(),
            'min_contig_length': self.min_contig_length.get()
        }
        
        missing_params = self.check_missing_parameters(config)
        if missing_params:
            self.fill_missing_parameters(missing_params, config)
        
        with open(self.config_file, 'w') as f:
            json.dump(config, f, indent=4)
        
        messagebox.showinfo("Success", "Parameters saved successfully!")

    def check_missing_parameters(self, config):
        missing = []
        essential_params = [
            'analysis_type', 'threads', 'min_length', 'min_quality',
            'emu_db', 'centrifuge_db', 'min_abundance', 'assembly_mode',
            'medaka_model', 'min_contig_length'
        ]
        for param in essential_params:
            if param not in config or not config[param]:
                missing.append(param)
        return missing

    def fill_missing_parameters(self, missing_params, config):
        for param in missing_params:
            if param == 'emu_db':
                if os.path.exists(self.emu_db_path):
                    config[param] = self.emu_db_path
                else:
                    messagebox.showerror("Error", f"No default Emu database found. Please specify manually.")
                    return
            elif param == 'centrifuge_db':
                if os.path.exists(self.centrifuge_db_path):
                    config[param] = self.centrifuge_db_path
                else:
                    messagebox.showerror("Error", f"No default Centrifuge database found. Please specify manually.")
                    return
            elif param == 'threads':
                config[param] = str(multiprocessing.cpu_count())
            else:
                messagebox.showerror("Error", f"Missing parameter: {param}. Please specify manually.")
                return

    def search_available_databases(self):
        if os.path.exists(self.emu_db_path):
            self.emu_db.set(self.emu_db_path)
            self.emu_db_entry.configure(state="normal")
            self.emu_db_entry.delete(0, tk.END)
            self.emu_db_entry.insert(0, self.emu_db_path)
            self.emu_db_entry.configure(state="readonly")
        
        if os.path.exists(self.centrifuge_db_path):
            self.centrifuge_db.set(self.centrifuge_db_path)
            self.centrifuge_db_entry.configure(state="normal")
            self.centrifuge_db_entry.delete(0, tk.END)
            self.centrifuge_db_entry.insert(0, self.centrifuge_db_path)
            self.centrifuge_db_entry.configure(state="readonly")

    def get_sorted_medaka_models(self):
        models = [
            "r1041_e82_400bps_sup_v5.0.0", "r1041_e82_400bps_sup_variant_v5.0.0",
            "r1041_e82_400bps_sup_v4.3.0", "r1041_e82_400bps_sup_variant_v4.3.0",
            "r1041_e82_400bps_sup_v4.2.0", "r1041_e82_400bps_sup_variant_v4.2.0",
            "r1041_e82_400bps_sup_v4.1.0", "r1041_e82_400bps_sup_variant_v4.1.0",
            "r1041_e82_400bps_sup_v4.0.0", "r1041_e82_400bps_sup_variant_v4.0.0",
            "r1041_e82_400bps_hac_v5.0.0", "r1041_e82_400bps_hac_variant_v5.0.0",
            "r1041_e82_400bps_hac_v4.3.0", "r1041_e82_400bps_hac_variant_v4.3.0",
            "r1041_e82_400bps_hac_v4.2.0", "r1041_e82_400bps_hac_variant_v4.2.0",
            "r1041_e82_400bps_hac_v4.1.0", "r1041_e82_400bps_hac_variant_v4.1.0",
            "r1041_e82_400bps_hac_v4.0.0",
            "r1041_e82_400bps_fast_g632", "r1041_e82_400bps_fast_variant_g632",
            "r1041_e82_400bps_hac_g632", "r1041_e82_400bps_hac_variant_g632",
            "r1041_e82_400bps_sup_g615", "r1041_e82_400bps_sup_variant_g615",
            "r1041_e82_400bps_hac_g615", "r1041_e82_400bps_hac_variant_g615",
            "r1041_e82_400bps_fast_g615", "r1041_e82_400bps_fast_variant_g615",
            "r1041_e82_260bps_sup_v4.1.0", "r1041_e82_260bps_sup_variant_v4.1.0",
            "r1041_e82_260bps_sup_v4.0.0",
            "r1041_e82_260bps_hac_v4.1.0", "r1041_e82_260bps_hac_variant_v4.1.0",
            "r1041_e82_260bps_hac_v4.0.0",
            "r1041_e82_260bps_fast_g632", "r1041_e82_260bps_fast_variant_g632",
            "r1041_e82_260bps_hac_g632", "r1041_e82_260bps_hac_variant_g632",
            "r1041_e82_260bps_sup_g632", "r1041_e82_260bps_sup_variant_g632",
            "r941_min_sup_g507", "r941_min_sup_variant_g507",
            "r941_min_hac_g507", "r941_min_hac_variant_g507",
            "r941_min_fast_g507", "r941_min_fast_variant_g507",
            "r941_prom_sup_g507", "r941_prom_sup_variant_g507",
            "r941_prom_hac_g507", "r941_prom_hac_variant_g507",
            "r941_prom_fast_g507", "r941_prom_fast_variant_g507",
            "r941_e81_sup_g514", "r941_e81_sup_variant_g514",
            "r941_e81_hac_g514", "r941_e81_hac_variant_g514",
            "r941_e81_fast_g514", "r941_e81_fast_variant_g514",
            "r103_sup_g507", "r103_sup_variant_g507",
            "r103_hac_g507", "r103_hac_variant_g507",
            "r103_fast_g507", "r103_fast_variant_g507",
            "r104_e81_sup_g610", "r104_e81_sup_variant_g610",
            "r104_e81_sup_g5015",
            "r104_e81_hac_g5015", "r104_e81_hac_variant_g5015",
            "r104_e81_fast_g5015", "r104_e81_fast_variant_g5015",
        ]
        return models

    def run_pipeline(self):
        analysis_type = self.analysis_type.get()
        params = {
            "threads": int(self.threads.get()),
            "min_length": int(self.min_length.get()),
            "min_quality": float(self.min_quality.get()),
            "emu_db": self.emu_db.get(),
            "centrifuge_db": self.centrifuge_db.get(),
            "min_abundance": float(self.min_abundance.get()),
            "assembly_mode": self.assembly_mode.get(),
            "is_isolate": self.is_isolate.get(),
            "min_contig_length": int(self.min_contig_length.get()),
            "medaka_model": self.medaka_model.get(),
        }
        self.gui.run_pipeline(analysis_type, params)