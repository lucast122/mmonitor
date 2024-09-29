import tkinter as tk
from threading import Thread, Lock
import customtkinter as ctk

class PipelinePopup(ctk.CTkToplevel):
    def __init__(self, parent, gui_ref):
        super().__init__()
        ctk.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"
        self.parent = parent
        self.gui = gui_ref
        self.taxonomy_nanopore_wgs = tk.BooleanVar()
        self.taxonomy_nanopore_16s_bool = tk.BooleanVar()
        self.assembly = tk.BooleanVar()
        self.correction = tk.BooleanVar()
        self.binning = tk.BooleanVar()
        self.annotation = tk.BooleanVar()
        self.kegg = tk.BooleanVar()
        self.test_mag_upload = tk.BooleanVar()  # New variable for test MAG upload
        self.geometry("440x420")  # Increased height to accommodate new checkbox
        self.minsize(440, 420)
        self.title("Select analysis steps to perform.")
        frame_taxonomy = ctk.CTkFrame(self, corner_radius=10)
        frame_functional = ctk.CTkFrame(self, corner_radius=10)
        frame_taxonomy.pack(pady=5, padx=10, fill="both", expand=True)
        frame_functional.pack(pady=5, padx=10, fill="both", expand=True)
        label_taxonomy = ctk.CTkLabel(frame_taxonomy, text="Taxonomic analysis")
        label_functional = ctk.CTkLabel(frame_functional, text="Functional analysis")        
        label_taxonomy.pack(pady=10)
        label_functional.pack(pady=10)
        # Taxonomy checkboxes
        taxonomy_wgs_checkbox = ctk.CTkCheckBox(frame_taxonomy, text='Taxonomy whole genome nanopore', variable=self.taxonomy_nanopore_wgs)
        taxonomy_wgs_checkbox.pack(pady=2)
        taxonomy_wgs_checkbox.bind("<Enter>", lambda e: self.show_tooltip(e, "Select this if your data is whole genome nanopore sequencing data (not targeted to 16S rRNA genes)."))
        taxonomy_wgs_checkbox.bind("<Leave>", lambda e: self.hide_tooltip(e))
        taxonomy_16s_checkbox = ctk.CTkCheckBox(frame_taxonomy, text='Taxonomy 16S-rRNA reads nanopore',
                                                variable=self.taxonomy_nanopore_16s_bool)
        taxonomy_16s_checkbox.pack(pady=2)
        taxonomy_16s_checkbox.bind("<Enter>", lambda e: self.show_tooltip(e, "Select this if your data is 16S-rRNA reads from nanopore sequencing data (targeted to 16S rRNA genes only)."))
        taxonomy_16s_checkbox.bind("<Leave>", lambda e: self.hide_tooltip(e))
        # Functional analysis checkboxes
        assembly_checkbox = ctk.CTkCheckBox(frame_functional, text='Assembly Pipeline', variable=self.assembly)
        assembly_checkbox.pack(pady=2)
        assembly_checkbox.bind("<Enter>", lambda e: self.show_tooltip(e, "Select this to assemble a metagenome. This will run the full assembly pipeline, results will be MAGs."))
        assembly_checkbox.bind("<Leave>", lambda e: self.hide_tooltip(e))
        ctk.CTkCheckBox(frame_functional, text='KEGG', variable=self.kegg).pack(pady=2)
        # Test MAG Upload checkbox
        ctk.CTkCheckBox(self, text='Test MAG Upload', variable=self.test_mag_upload).pack(pady=5)
        # Continue and Quit buttons
        continue_btn = ctk.CTkButton(self, text="Continue", command=self.run_analysis_pipeline, corner_radius=10)
        continue_btn.pack(pady=5)
        quit_btn = ctk.CTkButton(self, text="Quit", command=self.destroy, corner_radius=10)
        quit_btn.pack(pady=5)
        self.lock = Lock()  # Add a threading lock

    def show_tooltip(self, event, text):
        # Create a tooltip window
        self.tooltip = tk.Toplevel(self)
        self.tooltip.wm_overrideredirect(True)  # Remove window decorations
        self.tooltip.wm_geometry(f"+{event.x_root + 10}+{event.y_root + 10}")  # Position the tooltip
        # Create a label inside the tooltip window
        label = tk.Label(self.tooltip, text=text, background="black", relief="solid", borderwidth=1)
        label.pack()

    def hide_tooltip(self, event):
        if hasattr(self, 'tooltip'):
            self.tooltip.destroy()
            del self.tooltip

    def run_analysis_pipeline(self):
        if self.taxonomy_nanopore_wgs.get():
            self.gui.run_pipeline("taxonomy-wgs")
        if self.taxonomy_nanopore_16s_bool.get():
            self.gui.run_pipeline("taxonomy-16s")
        if self.assembly.get():
            self.gui.run_pipeline("assembly")
        if self.kegg.get():
            self.gui.run_pipeline("kegg")
        self.destroy()