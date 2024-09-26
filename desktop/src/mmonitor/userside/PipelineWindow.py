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
        assembly_checkbox.bind("<Enter>", lambda e: self.show_tooltip(e, "Select this to assemble a metagenome. This will run the full assembly pipeline including read filtering, de novo assembly with flye, assembly correction with medaka, contig binning with metabat2, and taxonomic assignment of bins."))
        assembly_checkbox.bind("<Leave>", lambda e: self.hide_tooltip(e))
        # ctk.CTkCheckBox(frame_functional, text='Correction', variable=self.correction).pack(pady=2)
        # ctk.CTkCheckBox(frame_functional, text='Binning', variable=self.binning).pack(pady=2)
        # ctk.CTkCheckBox(frame_functional, text='Annotation', variable=self.annotation).pack(pady=2)
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
        # if self.assembly.get() or self.correction.get():
        #     seq_file = filedialog.askopenfilename(title="Please select a sequencing file")

        # if (self.assembly.get() or self.correction.get() or
        #         self.annotation.get() or self.binning.get()):
        #     sample_name = self.gui.ask_sample_name()
        #     self.gui.functional_analysis_runner.check_software_avail()

        if self.taxonomy_nanopore_wgs.get():
            thread_wgs = Thread(target=self.gui.taxonomy_nanopore_wgs)
            thread_wgs.start()

        if self.taxonomy_nanopore_16s_bool.get():
            thread_16s = Thread(target=self.gui.taxonomy_nanopore_16s)
            thread_16s.start()

        if self.assembly.get():
            thread_assembly = Thread(target=self.gui.functional_pipeline)
            thread_assembly.start()
            pass

        if self.correction.get():
            self.gui.django_db.process_mags_and_upload("/Users/timo/Documents/spirito/mags", overwrite=True)

        if self.test_mag_upload.get():
            thread_mag_upload = Thread(target=self.upload_mag_thread)
            thread_mag_upload.start()

        # if self.binning.get():
        #     self.gui.functional_analysis_runner.run_binning(sample_name)
        # if self.annotation.get():
        #     bins_path = f"{ROOT}/src/resources/{sample_name}/bins/"
        #     self.gui.functional_analysis_runner.run_prokka(bins_path)

        # if only kegg analysis is selected then the user needs to chose the path to the annotations
        # if self.kegg.get() and not self.assembly.get() and not self.correction.get() and not self.binning.get() and not self.annotation.get():
        #     sample_name = self.gui.ask_sample_name()
        #     pipeline_out = filedialog.askdirectory(
        #         title="Please select the path to the prokka output (folder with tsv files with annotations).")
        #     pipeline_out = f"{pipeline_out}/"
        #     pipeline_out = f"{ROOT}/src/resources/pipeline_out/{sample_name}/"
        #     self.kegg_thread1 = Thread(target=self.functional_analysis_runner.create_keggcharter_input(pipeline_out))
        #     self.kegg_thread1.start()
        #     self.kegg_thread2 = Thread(
        #         self.functional_analysis_runner.run_keggcharter(pipeline_out, f"{pipeline_out}keggcharter.tsv"))
        #     self.kegg_thread2.start()
        #
        # self.destroy()

        # if kegg and annotation is chosen then the user only needs to select the sample name, then the tsv files from the results
        # of the annotations will be used as input for creating keggcharter input and creating kegg maps
    # if self.kegg.get() and self.annotation.get():
    #     sample_name = self.ask_sample_name()
    #     # pipeline_out = filedialog.askdirectory(title="Please select the path to the prokka output (folder with tsv files with annotations).")
    #     pipeline_out = f"{ROOT}/src/resources/pipeline_out/{sample_name}/"
    #     self.kegg_thread1 = Thread(target=self.functional_analysis_runner.create_keggcharter_input(pipeline_out))
    #     self.kegg_thread1.start()
    #     self.kegg_thread2 = Thread(self.functional_analysis_runner.run_keggcharter(pipeline_out, f"{pipeline_out}/keggcharter.tsv"))
    #     self.kegg_thread2.start()

    def upload_mag_thread(self):
        with self.lock:
            self.gui.django_db.upload_mag(
                name="p.caeni_test",
                taxonomy="Pseudoclavibacter caeni",  # Added a taxonomy field
                sample_name="test",
                fasta_file_path="/Users/timo/Downloads/GCA_008831125.1_ASM883112v1_genomic.fna",
                gff_file_path="/Users/timo/Downloads/genomic.gff"
            )