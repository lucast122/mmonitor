import os
import subprocess
import logging
import glob
import json
from typing import List, Optional, Dict, Any
from .PipelineStateTracker import PipelineStateTracker
from .ToolInstaller import ToolInstaller
from .MetaBatRunner import MetaBatRunner
from .FunctionalRunner import FunctionalRunner
from ..paths import SRC_DIR, LIB_DIR, RESOURCES_DIR

logger = logging.getLogger(__name__)

class AssemblyPipeline:
    def __init__(self, sample_dir: str, sample_name: str, db_interface=None):
        """Initialize the assembly pipeline
        
        Args:
            sample_dir (str): Directory for sample data and results
            sample_name (str): Name of the sample
            db_interface: Interface for database operations
        """
        self.sample_dir = sample_dir
        self.sample_name = sample_name
        self.db_interface = db_interface
        self.default_threads = os.cpu_count() or 1
        self.tool_installer = ToolInstaller()
        self.state_tracker = PipelineStateTracker(os.path.join(sample_dir, "pipeline_state.json"))
        self.functional_runner = FunctionalRunner()

    def run_assembly(self, input_files: List[str], threads: Optional[int] = None) -> bool:
        """Run the complete assembly pipeline
        
        Args:
            input_files (List[str]): List of input FASTQ files
            threads (Optional[int]): Number of threads to use
            
        Returns:
            bool: True if successful
        """
        if not threads:
            threads = self.default_threads

        try:
            # Create output directories
            bins_dir = os.path.join(self.sample_dir, "bins")
            os.makedirs(bins_dir, exist_ok=True)

            # Run assembly steps
            if not self.state_tracker.is_step_complete("concatenate"):
                print("Concatenating input files...")
                if len(input_files) > 1:
                    concat_file = os.path.join(self.sample_dir, "concatenated.fastq")
                    self._concatenate_files(input_files, concat_file)
                    input_file = concat_file
                else:
                    input_file = input_files[0]
                self.state_tracker.mark_step_complete("concatenate", input_file)
            else:
                input_file = self.state_tracker.get_step_output("concatenate")

            # Run Flye assembly
            if not self.state_tracker.is_step_complete("flye_assembly"):
                print("Running Flye assembly...")
                assembly_file = os.path.join(self.sample_dir, "assembly.fasta")
                self._run_flye(input_file, assembly_file, threads)
                self.state_tracker.mark_step_complete("flye_assembly", assembly_file)
            else:
                assembly_file = self.state_tracker.get_step_output("flye_assembly")

            # Run Medaka polishing
            if not self.state_tracker.is_step_complete("medaka_polish"):
                print("Running Medaka polishing...")
                polished_file = os.path.join(self.sample_dir, "polished.fasta")
                self._run_medaka(input_file, assembly_file, polished_file, threads)
                self.state_tracker.mark_step_complete("medaka_polish", polished_file)
            else:
                polished_file = self.state_tracker.get_step_output("medaka_polish")

            # Run MetaBAT2 binning
            if not self.state_tracker.is_step_complete("metabat_binning"):
                print("Running MetaBAT2 binning...")
                metabat_runner = MetaBatRunner()
                metabat_runner.run_metabat2(polished_file, input_file, bins_dir, threads)
                self.state_tracker.mark_step_complete("metabat_binning", bins_dir)

            # Run CheckM2 quality assessment
            if not self.state_tracker.is_step_complete("checkm2_quality"):
                print("Running CheckM2 quality assessment...")
                checkm2_dir = os.path.join(self.sample_dir, "checkm2")
                os.makedirs(checkm2_dir, exist_ok=True)
                self._run_checkm2(bins_dir, checkm2_dir, threads)
                self.state_tracker.mark_step_complete("checkm2_quality", checkm2_dir)

            # Run GTDB-TK taxonomy classification
            if not self.state_tracker.is_step_complete("gtdbtk_taxonomy"):
                print("Running GTDB-TK classification...")
                gtdbtk_dir = os.path.join(self.sample_dir, "gtdbtk")
                os.makedirs(gtdbtk_dir, exist_ok=True)
                self._run_gtdbtk(bins_dir, gtdbtk_dir, threads)
                self.state_tracker.mark_step_complete("gtdbtk_taxonomy", gtdbtk_dir)

            # Run Bakta annotation
            if not self.state_tracker.is_step_complete("bakta_annotation"):
                print("Running Bakta annotation...")
                bakta_dir = os.path.join(self.sample_dir, "bakta")
                os.makedirs(bakta_dir, exist_ok=True)
                self._run_bakta(bins_dir, bakta_dir, threads)
                self.state_tracker.mark_step_complete("bakta_annotation", bakta_dir)

            # Run KEGGCharter
            if not self.state_tracker.is_step_complete("keggcharter"):
                print("Running KEGGCharter...")
                kegg_dir = os.path.join(self.sample_dir, "kegg")
                os.makedirs(kegg_dir, exist_ok=True)
                
                # Process each bin's Bakta output with KEGGCharter
                bakta_dir = self.state_tracker.get_step_output("bakta_annotation")
                for bin_dir in glob.glob(os.path.join(bakta_dir, "*")):
                    bin_name = os.path.basename(bin_dir)
                    kegg_out = os.path.join(kegg_dir, bin_name)
                    os.makedirs(kegg_out, exist_ok=True)
                    
                    # Create KEGGCharter input from Bakta output
                    keggcharter_input = os.path.join(bin_dir, "keggcharter.tsv")
                    self.functional_runner.create_keggcharter_input(bin_dir)
                    
                    # Run KEGGCharter
                    self.functional_runner.run_keggcharter(kegg_out, keggcharter_input)
                
                self.state_tracker.mark_step_complete("keggcharter", kegg_dir)

            # Upload results to server
            if not self.state_tracker.is_step_complete("upload_to_server"):
                print("Uploading results to server...")
                self._upload_results()
                self.state_tracker.mark_step_complete("upload_to_server", True)

            return True

        except Exception as e:
            print(f"Error in assembly pipeline: {str(e)}")
            return False

    def _concatenate_files(self, input_files: List[str], output_file: str) -> None:
        """Concatenate input FASTQ files using Rust implementation, falls back to Python if Rust fails
        
        Args:
            input_files: List of input FASTQ files
            output_file: Output concatenated file
        """
        try:
            # Try Rust implementation first
            from .MMonitorCMD import MMonitorCMD
            cmd = MMonitorCMD()
            print("Using Rust implementation for FASTQ concatenation...")
            cmd.concatenate_fastq_files(input_files, output_file)
            
        except Exception as e:
            print(f"Rust concatenation failed: {str(e)}")
            print("Falling back to Python implementation...")
            
            # Python fallback
            with open(output_file, 'wb') as outfile:
                for file in input_files:
                    with open(file, 'rb') as infile:
                        outfile.write(infile.read())

    def _run_flye(self, input_file: str, output_file: str, threads: int) -> None:
        """Run Flye assembly"""
        cmd = f"flye --nano-raw {input_file} --out-dir {os.path.dirname(output_file)} --threads {threads}"
        subprocess.run(cmd, shell=True, check=True)

    def _run_medaka(self, input_file: str, assembly_file: str, output_file: str, threads: int) -> None:
        """Run Medaka polishing"""
        cmd = f"medaka_consensus -i {input_file} -d {assembly_file} -o {os.path.dirname(output_file)} -t {threads}"
        subprocess.run(cmd, shell=True, check=True)

    def _run_checkm2(self, bins_dir: str, output_dir: str, threads: int) -> None:
        """Run CheckM2 quality assessment"""
        cmd = f"checkm2 predict -i {bins_dir} -o {output_dir}/quality_report.tsv -t {threads} --force"
        subprocess.run(cmd, shell=True, check=True)

    def _run_gtdbtk(self, bins_dir: str, output_dir: str, threads: int) -> None:
        """Run GTDB-TK classification"""
        if 'GTDBTK_DATA_PATH' not in os.environ:
            db_path = self._check_gtdbtk_db()
            if not db_path:
                raise RuntimeError("GTDB-TK database not found")
            os.environ['GTDBTK_DATA_PATH'] = db_path
        
        cmd = f"gtdbtk classify_wf --genome_dir {bins_dir} --out_dir {output_dir} --cpus {threads} --extension fa"
        subprocess.run(cmd, shell=True, check=True)

    def _run_bakta(self, bins_dir: str, output_dir: str, threads: int) -> None:
        """Run Bakta annotation"""
        db_path = self._check_bakta_db()
        if not db_path:
            raise RuntimeError("Bakta database not found")
        
        for bin_file in glob.glob(os.path.join(bins_dir, "*.fa")):
            bin_name = os.path.basename(bin_file).rsplit('.', 1)[0]
            bin_output = os.path.join(output_dir, bin_name)
            os.makedirs(bin_output, exist_ok=True)
            
            cmd = f"bakta --db {db_path} --output {bin_output} --prefix {bin_name} --threads {threads} {bin_file}"
            subprocess.run(cmd, shell=True, check=True)

    def _check_gtdbtk_db(self) -> Optional[str]:
        """Check GTDB-TK database location"""
        db_paths = [
            os.path.expanduser("~/.gtdbtk/db"),
            "/opt/conda/db/gtdbtk",
            os.path.join(LIB_DIR, "gtdbtk_db")
        ]
        return next((path for path in db_paths if os.path.exists(os.path.join(path, "metadata"))), None)

    def _check_bakta_db(self) -> Optional[str]:
        """Check Bakta database location"""
        db_paths = [
            os.path.expanduser("~/.bakta/db"),
            "/opt/conda/db/bakta",
            os.path.join(LIB_DIR, "bakta_db")
        ]
        return next((path for path in db_paths if os.path.exists(os.path.join(path, "version.json"))), None)

    def _upload_results(self) -> None:
        """Upload results to Django server"""
        if not self.db_interface:
            return

        try:
            # Get paths from state tracker
            bins_dir = self.state_tracker.get_step_output("metabat_binning")
            checkm2_dir = self.state_tracker.get_step_output("checkm2_quality")
            gtdbtk_dir = self.state_tracker.get_step_output("gtdbtk_taxonomy")
            bakta_dir = self.state_tracker.get_step_output("bakta_annotation")
            kegg_dir = self.state_tracker.get_step_output("keggcharter")

            # Read CheckM2 quality results
            quality_results = {}
            with open(os.path.join(checkm2_dir, "quality_report.tsv"), 'r') as f:
                next(f)  # Skip header
                for line in f:
                    bin_name, completeness, contamination, strain_het = line.strip().split('\t')
                    quality_results[bin_name] = {
                        'completeness': float(completeness),
                        'contamination': float(contamination),
                        'strain_heterogeneity': float(strain_het)
                    }

            # Read GTDB-TK taxonomy results
            taxonomy_results = {}
            gtdbtk_summary = os.path.join(gtdbtk_dir, "classify", "gtdbtk.summary.tsv")
            if os.path.exists(gtdbtk_summary):
                with open(gtdbtk_summary, 'r') as f:
                    next(f)  # Skip header
                    for line in f:
                        fields = line.strip().split('\t')
                        bin_name = os.path.basename(fields[0])
                        taxonomy = fields[1] if len(fields) > 1 else "Unknown"
                        taxonomy_results[bin_name] = taxonomy

            # Process each bin
            for bin_path in glob.glob(os.path.join(bins_dir, "*.fa")):
                bin_name = os.path.basename(bin_path)
                bin_prefix = bin_name.rsplit('.', 1)[0]
                
                # Get quality metrics and taxonomy
                quality = quality_results.get(bin_name, {})
                taxonomy = taxonomy_results.get(bin_name, "Unknown")
                
                # Get Bakta annotation stats
                bakta_stats = {}
                bakta_json = os.path.join(bakta_dir, bin_prefix, f"{bin_prefix}.json")
                if os.path.exists(bakta_json):
                    with open(bakta_json) as f:
                        stats = json.load(f)
                        bakta_stats = {
                            'gene_count': stats.get('genes', {}).get('total', 0),
                            'trna_count': stats.get('genes', {}).get('trna', 0),
                            'rrna_count': stats.get('genes', {}).get('rrna', 0),
                            'cds_count': stats.get('genes', {}).get('cds', 0)
                        }
                
                # Get KEGG maps
                kegg_maps = []
                kegg_bin_dir = os.path.join(kegg_dir, bin_prefix)
                if os.path.exists(kegg_bin_dir):
                    for map_file in glob.glob(os.path.join(kegg_bin_dir, "*.png")):
                        map_id = os.path.basename(map_file).split('.')[0]
                        kegg_maps.append({
                            'map_id': map_id,
                            'description': f"KEGG map {map_id}",
                            'image_path': map_file
                        })
                
                # Upload to Django server
                self.db_interface.upload_mag(
                    name=f"{self.sample_name}_{bin_prefix}",
                    taxonomy=taxonomy,
                    sample_name=self.sample_name,
                    gff_file_path=os.path.join(bakta_dir, bin_prefix, f"{bin_prefix}.gff3"),
                    fasta_file_path=bin_path,
                    completeness=quality.get('completeness'),
                    contamination=quality.get('contamination'),
                    strain_heterogeneity=quality.get('strain_heterogeneity'),
                    **bakta_stats
                )
                
                # Upload KEGG maps
                for kegg_map in kegg_maps:
                    self.db_interface.upload_kegg_map(
                        mag_name=f"{self.sample_name}_{bin_prefix}",
                        **kegg_map
                    )

        except Exception as e:
            print(f"Error uploading results: {str(e)}")
            raise