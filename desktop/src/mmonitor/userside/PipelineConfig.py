import os
import json
import customtkinter as ctk
from tkinter import messagebox
from ..paths import SRC_DIR, RESOURCES_DIR
import tkinter as tk
import filedialog
import multiprocessing
from mmonitor.userside.utils import create_tooltip
from mmonitor.userside.MMonitorCMD import MMonitorCMD

class PipelineConfig(ctk.CTkFrame):
    def __init__(self, parent, main_window=None):
        super().__init__(parent)
        self.parent = parent
        self.main_window = main_window
        
        # Initialize configuration dictionary
        self.pipeline_config = {}
        
        # Load configuration
        self.config_file = os.path.join(SRC_DIR, "resources", "pipeline_config.json")
        
        # Initialize variables first
        self.analysis_type = ctk.StringVar(value='taxonomy-wgs')  # Default value
        self.threads = ctk.StringVar(value='4')
        self.min_length = ctk.StringVar(value='1000')
        self.min_quality = ctk.StringVar(value='10')
        self.emu_db = ctk.StringVar(value='')
        self.centrifuger_db = ctk.StringVar(value='')
        self.min_abundance = ctk.StringVar(value='0.01')
        
        # Initialize assembly variables
        self.assembly_mode = ctk.StringVar(value='nano-raw')
        self.medaka_model = ctk.StringVar(value='r1041_e82_400bps_sup_v5.0.0')
        self.is_isolate = ctk.BooleanVar(value=False)
        self.min_contig_length = ctk.StringVar(value='1000')
        
        # Now load config and update variables
        self.load_config()
        
        # Create widgets
        self.create_widgets()

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
            ("Centrifuger Database:", self.centrifuger_db),
            ("Min abundance:", self.min_abundance)
        ])

        # Assembly parameters
        self.create_section(content_frame, "Assembly Parameters", [
            ("Assembly mode:", self.assembly_mode, ["nano-raw", "nano-corrected","nano-hq"], "option"),
            ("Medaka model:", self.medaka_model, self.get_sorted_medaka_models(), "option"),
            ("Isolate Mode:", self.is_isolate, None, "checkbox"),
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

        # Group fields that should appear on the same line
        if title == "General Parameters":
            # Create a frame for the first line (Analysis Type)
            line_frame = ctk.CTkFrame(section_frame, fg_color="transparent")
            line_frame.pack(fill="x", pady=1)
            
            # Analysis Type (full width)
            label = ctk.CTkLabel(line_frame, text="Analysis Type:", width=100, anchor="w")
            label.pack(side="left", padx=2)
            widget = ctk.CTkOptionMenu(line_frame, variable=fields[0][1], values=fields[0][2])
            widget.pack(side="left", padx=2)
            
            # Create a frame for the second line (CPU, Length, Quality)
            line_frame = ctk.CTkFrame(section_frame, fg_color="transparent")
            line_frame.pack(fill="x", pady=1)
            
            # Pack CPU cores, Min length, and Min quality on one line
            for field in fields[1:4]:  # CPU cores, min length, min quality
                field_container = ctk.CTkFrame(line_frame, fg_color="transparent")
                field_container.pack(side="left", expand=True, fill="x", padx=5)
                
                label = ctk.CTkLabel(field_container, text=field[0], anchor="w")
                label.pack(side="left")
                widget = ctk.CTkEntry(field_container, textvariable=field[1], width=80)
                widget.pack(side="left", padx=2)
                create_tooltip(widget, f"Set the {field[0].lower()[:-1]}")

        elif title == "Assembly Parameters":
            # First line: Assembly mode, Medaka model
            line_frame = ctk.CTkFrame(section_frame, fg_color="transparent")
            line_frame.pack(fill="x", pady=1)
            
            # Assembly mode
            label = ctk.CTkLabel(line_frame, text=fields[0][0], width=100, anchor="w")
            label.pack(side="left", padx=2)
            widget = ctk.CTkOptionMenu(line_frame, variable=fields[0][1], values=fields[0][2])
            widget.pack(side="left", padx=(2, 15))
            
            # Medaka model
            label = ctk.CTkLabel(line_frame, text=fields[1][0], width=100, anchor="w")
            label.pack(side="left", padx=2)
            widget = ctk.CTkOptionMenu(line_frame, variable=fields[1][1], values=fields[1][2])
            widget.pack(side="left", padx=2)
            
            # Second line: Isolate Mode, Min contig length
            line_frame = ctk.CTkFrame(section_frame, fg_color="transparent")
            line_frame.pack(fill="x", pady=1)
            
            # Isolate Mode
            label = ctk.CTkLabel(line_frame, text=fields[2][0], width=100, anchor="w")
            label.pack(side="left", padx=2)
            widget = ctk.CTkCheckBox(line_frame, text="", variable=fields[2][1])
            widget.pack(side="left", padx=(2, 15))
            
            # Min contig length
            label = ctk.CTkLabel(line_frame, text=fields[3][0], width=100, anchor="w")
            label.pack(side="left", padx=2)
            widget = ctk.CTkEntry(line_frame, textvariable=fields[3][1], width=80)
            widget.pack(side="left", padx=2)

        else:  # Taxonomy Parameters
            for field in fields:
                field_frame = ctk.CTkFrame(section_frame)
                field_frame.pack(fill="x", pady=1)
                
                label_text = "Centrifuger Database:" if "Centrifuger" in field[0] else "Emu Database:" if "Emu" in field[0] else "Min Abundance:" if "abundance" in field[0] else field[0]
            
                label = ctk.CTkLabel(field_frame, text=label_text, width=100, anchor="w")
                label.pack(side="left", padx=2)
                
                if field[0] in ["Emu Database:", "Centrifuger Database:"]:
                    db_frame = ctk.CTkFrame(field_frame)
                    db_frame.pack(side="left", expand=True, fill="x", padx=2)
                    
                    widget = ctk.CTkEntry(db_frame, textvariable=field[1])
                    widget.pack(side="left", expand=True, fill="x")
                    
                    browse_btn = ctk.CTkButton(db_frame, text="Browse", 
                                           command=lambda f=field[0]: self.browse_database(f))
                    browse_btn.pack(side="right", padx=2)
                    
                    if field[0] == "Emu Database:":
                        self.emu_db_entry = widget
                    else:
                        self.centrifuger_db_entry = widget
                else:
                    widget = ctk.CTkEntry(field_frame, textvariable=field[1], width=80)
                    widget.pack(side="left", padx=2)
                
                create_tooltip(widget, f"Set the {field[0].lower()[:-1]}")

    def browse_database(self, db_type):
        """Browse for database directory and check for valid database files"""
        print("\n=== Database Selection Debug ===")
        print(f"Opening file dialog for {db_type} database selection...")
        
        db_path = filedialog.askdirectory(title=f"Select {db_type} Database Directory")
        print(f"Selected database path: {db_path}")
        
        if db_path:
            if db_type == "Centrifuger Database:":  # Match the exact string from create_section
                print(f"Looking for Centrifuger database files in: {db_path}")
                
                # Find all .cfr files recursively
                cfr_files = []
                for root, _, files in os.walk(db_path):
                    for file in files:
                        if file.endswith('.cfr'):
                            cfr_files.append(os.path.join(root, file))
                
                print(f"Found {len(cfr_files)} .cfr files: {cfr_files}")
                
                if cfr_files:
                    # Get the common prefix from the first file
                    first_file = cfr_files[0]
                    prefix = '.'.join(first_file.split('.')[:-2])  # Remove .N.cfr
                    print(f"Extracted prefix: {prefix}")
                    
                    # Update both the StringVar and the Entry widget
                    self.centrifuger_db.set(prefix)
                    if hasattr(self, 'centrifuger_db_entry'):
                        print("Updating entry widget with prefix...")
                        self.centrifuger_db_entry.delete(0, tk.END)
                        self.centrifuger_db_entry.insert(0, prefix)
                        print(f"Entry widget updated with: {prefix}")
                        
                        # Update config and save
                        self.pipeline_config['centrifuger_db'] = prefix  # Changed from self.config
                        self.save_parameters()
                        print("Configuration saved")
                    else:
                        print("Error: centrifuger_db_entry widget not found!")
                else:
                    print("No .cfr files found!")
                    messagebox.showerror("Error", 
                        "No valid Centrifuger database found in selected directory.\n"
                        "Directory should contain .cfr files with a common prefix.")
            
            elif db_type == "Emu Database:":
                self.emu_db.set(db_path)
                if hasattr(self, 'emu_db_entry'):
                    print("Updating EMU entry widget...")
                    self.emu_db_entry.delete(0, tk.END)
                    self.emu_db_entry.insert(0, db_path)
                    print(f"Entry widget updated with: {db_path}")
                    
                    # Update config and save
                    self.pipeline_config['emu_db'] = db_path  # Changed from self.config
                    self.save_parameters()
                    print("Configuration saved")
                else:
                    print("Error: emu_db_entry widget not found!")
        else:
            print("No directory selected")
        
        print("=== Database Selection Complete ===\n")

    def find_database_prefix(self, db_dir):
        """Find Centrifuger database prefix from .cfr files"""
        print(f"Scanning directory: {db_dir}")
        try:
            # List all files in directory
            files = os.listdir(db_dir)
            print(f"Files found: {files}")
            
            # Find all .cfr files
            cfr_files = [f for f in files if f.endswith('.cfr')]
            if not cfr_files:
                print("No .cfr files found in", db_dir)
                return None
            
            # Extract common prefix from .cfr files
            # Example: if files are ['centrifuger.1.cfr', 'centrifuger.2.cfr'],
            # prefix should be 'centrifuger'
            prefixes = set()
            for cfr_file in cfr_files:
                # Split by dots and take all parts except the last two (number and extension)
                parts = cfr_file.split('.')
                if len(parts) >= 3:  # Ensure we have at least prefix.number.cfr
                    prefix = '.'.join(parts[:-2])  # Join all parts except last two
                    if prefix:
                        prefixes.add(prefix)
            
            if len(prefixes) == 1:
                prefix = prefixes.pop()
                full_path = os.path.join(db_dir, prefix)
                print(f"Found database prefix: {full_path}")
                return full_path
            elif len(prefixes) > 1:
                print(f"Multiple database prefixes found: {prefixes}")
                return None
            else:
                print("No valid database prefix found")
                return None
            
        except Exception as e:
            print(f"Error scanning directory: {e}")
            return None

    def show_database_window(self):
        # Import DatabaseWindow here instead
        from mmonitor.userside.DatabaseWindow import DatabaseWindow
        if not hasattr(self.parent, 'database_window'):
            self.parent.database_window = DatabaseWindow(self.parent)
            self.parent.database_window.pack(fill="both", expand=True)
        self.parent.database_window.lift()

    def load_config(self):
        """Load existing configuration"""
        if os.path.exists(self.config_file):
            try:
                with open(self.config_file, 'r') as f:
                    config = json.load(f)
                    print(f"Loading pipeline config: {config}")
                    
                    # Update variables from config
                    if 'analysis_type' in config:
                        self.analysis_type.set(config['analysis_type'])
                    if 'threads' in config:
                        self.threads.set(config['threads'])
                    if 'min_length' in config:
                        self.min_length.set(config['min_length'])
                    if 'min_quality' in config:
                        self.min_quality.set(config['min_quality'])
                    if 'emu_db' in config:
                        self.emu_db.set(config['emu_db'])
                    if 'centrifuger_db' in config:
                        self.centrifuger_db.set(config['centrifuger_db'])
                    if 'min_abundance' in config:
                        self.min_abundance.set(config['min_abundance'])
                    if 'assembly_mode' in config:
                        self.assembly_mode.set(config['assembly_mode'])
                    if 'medaka_model' in config:
                        self.medaka_model.set(config['medaka_model'])
                    if 'is_isolate' in config:
                        self.is_isolate.set(config['is_isolate'])
                    if 'min_contig_length' in config:
                        self.min_contig_length.set(config['min_contig_length'])
                    
            except Exception as e:
                print(f"Error loading pipeline config: {e}")
                # Initialize with defaults if loading fails
                self.save_parameters()
        else:
            self.save_parameters()

    def update_widgets_from_config(self):
        for key, value in self.config.items():
            if hasattr(self, key):
                getattr(self, key).set(value)

    def save_parameters(self):
        """Save parameters to config file"""
        # Update pipeline_config with current values
        self.pipeline_config.update({  # Changed from self.config
            'analysis_type': self.analysis_type.get(),
            'threads': self.threads.get(),
            'min_length': self.min_length.get(),
            'min_quality': self.min_quality.get(),
            'emu_db': self.emu_db.get(),
            'centrifuger_db': self.centrifuger_db.get(),
            'min_abundance': self.min_abundance.get(),
            'assembly_mode': self.assembly_mode.get(),
            'medaka_model': self.medaka_model.get(),
            'is_isolate': self.is_isolate.get(),
            'min_contig_length': self.min_contig_length.get()
        })
        
        # Save to file
        with open(self.config_file, 'w') as f:
            json.dump(self.pipeline_config, f, indent=4)  # Changed from self.config

    def search_available_databases(self):
        """Search for available databases in default locations"""
        if os.path.exists(self.emu_db_path):
            self.emu_db.set(self.emu_db_path)
        
        if os.path.exists(self.centrifuger_db_path):
            prefix = self.find_database_prefix(self.centrifuger_db_path)
            if prefix:
                full_path = os.path.join(self.centrifuger_db_path, prefix)
                self.centrifuger_db.set(full_path)
                print(f"Found Centrifuger database with prefix: {prefix}")
                print(f"Full path: {full_path}")

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
            "centrifuger_db": self.centrifuger_db.get(),
            "min_abundance": float(self.min_abundance.get()),
            "assembly_mode": self.assembly_mode.get(),
            "is_isolate": self.is_isolate.get(),
            "min_contig_length": int(self.min_contig_length.get()),
            "medaka_model": self.medaka_model.get(),
        }
        self.gui.run_pipeline(analysis_type, params)
