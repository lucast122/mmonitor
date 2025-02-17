import os
import subprocess
import logging
import glob
import gzip
import shutil
import tempfile
import argparse
from typing import List, Optional, Dict, Any
from .PipelineStateTracker import PipelineStateTracker
from .ToolInstaller import ToolInstaller
from .MetaBatRunner import MetaBatRunner
from .FunctionalRunner import FunctionalRunner
from .FileUtils import concatenate_fastq_files
from ..paths import SRC_DIR, LIB_DIR, RESOURCES_DIR

logger = logging.getLogger(__name__)

class AssemblyPipeline:
    def __init__(self, output_dir: str, sample_name: str, django_db=None):
        """Initialize assembly pipeline
        
        Args:
            output_dir: Base output directory
            sample_name: Name of the sample
            django_db: Optional Django database interface
        """
        self.output_dir = output_dir
        self.sample_name = sample_name
        self.django_db = django_db
        
        # Create sample directory
        self.sample_dir = os.path.join(output_dir, sample_name)
        os.makedirs(self.sample_dir, exist_ok=True)
        
        # Initialize pipeline state tracker
        self.state_tracker = PipelineStateTracker(output_dir, sample_name)
        self.state_tracker.initialize_pipeline_state()
        
        # Initialize runners
        self.functional_runner = FunctionalRunner()
        self.tool_installer = ToolInstaller()
        
        # Create temp directory for intermediate files
        self.temp_dir = tempfile.mkdtemp(prefix=f"{self.sample_name}_")
        
        # Initialize default threads
        self.default_threads = os.cpu_count() or 1

    def run_assembly(self, input_files: List[str], threads: Optional[int] = None) -> bool:
        """Run the assembly pipeline
        
        Args:
            input_files: List of input fastq files
            threads: Number of threads to use, defaults to config value
            
        Returns:
            bool: True if assembly was successful, False otherwise
        """
        try:
            # Initialize pipeline state
            self.state_tracker.initialize_pipeline_state()
            
            # Create temp directory for intermediate files
            with tempfile.TemporaryDirectory() as temp_dir:
                # Concatenate input files
                input_file = os.path.join(temp_dir, f"{self.sample_name}_concat.fastq.gz")
                
                try:
                    # Run concatenation using FileUtils
                    print(f"Concatenating {len(input_files)} files...")
                    concatenate_fastq_files(input_files, input_file, threads or self.default_threads)
                    
                    if not os.path.exists(input_file):
                        raise FileNotFoundError(f"Concatenated file not found: {input_file}")
                    
                    # Run Flye assembly
                    assembly_file = os.path.join(self.sample_dir, "assembly.fasta")
                    if not self.state_tracker.is_step_complete("flye_assembly"):
                        self._run_flye(input_file, assembly_file, threads)
                    
                    # Run Medaka polishing if assembly succeeded
                    if os.path.exists(assembly_file):
                        polished_file = os.path.join(self.sample_dir, "polished_assembly.fasta")
                        if not self.state_tracker.is_step_complete("medaka_polish"):
                            self._run_medaka(input_file, assembly_file, polished_file, threads)
                            
                    # Run MetaBAT2 binning
                    if not self.state_tracker.is_step_complete("metabat_binning"):
                        print("Running MetaBAT2 binning...")
                        metabat_runner = MetaBatRunner()
                        try:
                            metabat_runner.run_metabat2(polished_file, input_file, os.path.join(self.sample_dir, "bins"), threads)
                            self.state_tracker.mark_step_complete("metabat_binning", os.path.join(self.sample_dir, "bins"))
                        except Exception as e:
                            self.state_tracker.mark_step_failed("metabat_binning", str(e))
                            raise
                    else:
                        bins_dir = self.state_tracker.get_step_output("metabat_binning")

                    # Run CheckM2 quality assessment
                    if not self.state_tracker.is_step_complete("checkm2_quality"):
                        print("Running CheckM2 quality assessment...")
                        checkm2_dir = os.path.join(self.sample_dir, "checkm2")
                        os.makedirs(checkm2_dir, exist_ok=True)
                        try:
                            self._run_checkm2(bins_dir, checkm2_dir, threads)
                            self.state_tracker.mark_step_complete("checkm2_quality", checkm2_dir)
                        except Exception as e:
                            self.state_tracker.mark_step_failed("checkm2_quality", str(e))
                            raise
                    else:
                        checkm2_dir = self.state_tracker.get_step_output("checkm2_quality")

                    # Run GTDB-TK taxonomy classification
                    if not self.state_tracker.is_step_complete("gtdbtk_taxonomy"):
                        print("Running GTDB-TK classification...")
                        gtdbtk_dir = os.path.join(self.sample_dir, "gtdbtk")
                        os.makedirs(gtdbtk_dir, exist_ok=True)
                        try:
                            self._run_gtdbtk(bins_dir, gtdbtk_dir, threads)
                            self.state_tracker.mark_step_complete("gtdbtk_taxonomy", gtdbtk_dir)
                        except Exception as e:
                            self.state_tracker.mark_step_failed("gtdbtk_taxonomy", str(e))
                            raise
                    else:
                        gtdbtk_dir = self.state_tracker.get_step_output("gtdbtk_taxonomy")

                    # Run Bakta annotation
                    if not self.state_tracker.is_step_complete("bakta_annotation"):
                        print("Running Bakta annotation...")
                        bakta_dir = os.path.join(self.sample_dir, "bakta")
                        os.makedirs(bakta_dir, exist_ok=True)
                        try:
                            self._run_bakta(bins_dir, bakta_dir, threads)
                            self.state_tracker.mark_step_complete("bakta_annotation", bakta_dir)
                        except Exception as e:
                            self.state_tracker.mark_step_failed("bakta_annotation", str(e))
                            raise
                    else:
                        bakta_dir = self.state_tracker.get_step_output("bakta_annotation")

                    # Run KEGGCharter
                    if not self.state_tracker.is_step_complete("keggcharter"):
                        print("Running KEGGCharter...")
                        kegg_dir = os.path.join(self.sample_dir, "kegg")
                        os.makedirs(kegg_dir, exist_ok=True)
                        
                        # Process each bin's Bakta output with KEGGCharter
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
                        try:
                            self._upload_results()
                            self.state_tracker.mark_step_complete("upload_to_server", True)
                        except Exception as e:
                            self.state_tracker.mark_step_failed("upload_to_server", str(e))
                            raise

                    # Finalize pipeline state
                    self.state_tracker.finalize_pipeline_state()
                    return True
                    
                except Exception as e:
                    error_msg = str(e)
                    print(f"Error in assembly pipeline: {error_msg}")
                    self.state_tracker.mark_step_failed("assembly", error_msg)
                    return False
                    
        except Exception as e:
            error_msg = str(e)
            print(f"Error in assembly pipeline: {error_msg}")
            self.state_tracker.mark_step_failed("assembly", error_msg)
            return False

    def _concatenate_files(self, input_files: List[str], output_file: str) -> None:
        """Concatenate input FASTQ files and gzip the output
        
        Args:
            input_files: List of input FASTQ files
            output_file: Output concatenated file (should end with .gz)
        """
        try:
            with gzip.open(output_file, 'wb') as outf:
                for file in input_files:
                    if file.endswith('.gz'):
                        with gzip.open(file, 'rb') as inf:
                            shutil.copyfileobj(inf, outf)
                    else:
                        with open(file, 'rb') as inf:
                            shutil.copyfileobj(inf, outf)
        except Exception as e:
            print(f"Error concatenating files: {e}")
            raise

    def _run_flye(self, input_file: str, output_file: str, threads: int) -> None:
        """Run Flye assembly"""
        # Verify input file exists
        if not os.path.exists(input_file):
            raise RuntimeError(f"Input file not found: {input_file}")
            
        # Create output directory
        out_dir = os.path.dirname(output_file)
        os.makedirs(out_dir, exist_ok=True)
        
        # Run Flye with error capture
        cmd = f"flye --nano-raw {input_file} --out-dir {out_dir} --threads {threads}"
        try:
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            error_msg = f"Flye assembly failed: {e.stderr if e.stderr else str(e)}"
            raise RuntimeError(error_msg)

    def _run_medaka(self, input_file: str, assembly_file: str, output_file: str, threads: int) -> None:
        """Run Medaka polishing"""
        # Verify input file exists
        if not os.path.exists(input_file):
            raise RuntimeError(f"Input file not found: {input_file}")
            
        # Verify assembly file exists
        if not os.path.exists(assembly_file):
            raise RuntimeError(f"Assembly file not found: {assembly_file}")
            
        # Create output directory
        out_dir = os.path.dirname(output_file)
        os.makedirs(out_dir, exist_ok=True)
        
        # Run Medaka with error capture
        cmd = f"medaka_consensus -i {input_file} -d {assembly_file} -o {out_dir} -t {threads}"
        try:
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            error_msg = f"Medaka polishing failed: {e.stderr if e.stderr else str(e)}"
            raise RuntimeError(error_msg)

    def _run_checkm2(self, bins_dir: str, output_dir: str, threads: int) -> None:
        """Run CheckM2 quality assessment"""
        # Verify bins directory exists
        if not os.path.exists(bins_dir):
            raise RuntimeError(f"Bins directory not found: {bins_dir}")
            
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Run CheckM2 with error capture
        cmd = f"checkm2 predict -i {bins_dir} -o {output_dir}/quality_report.tsv -t {threads} --force"
        try:
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            error_msg = f"CheckM2 quality assessment failed: {e.stderr if e.stderr else str(e)}"
            raise RuntimeError(error_msg)

    def _run_gtdbtk(self, bins_dir: str, output_dir: str, threads: int) -> None:
        """Run GTDB-TK classification"""
        # Verify bins directory exists
        if not os.path.exists(bins_dir):
            raise RuntimeError(f"Bins directory not found: {bins_dir}")
            
        # Verify GTDB-TK database exists
        if 'GTDBTK_DATA_PATH' not in os.environ:
            db_path = self._check_gtdbtk_db()
            if not db_path:
                raise RuntimeError("GTDB-TK database not found")
            os.environ['GTDBTK_DATA_PATH'] = db_path
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Run GTDB-TK with error capture
        cmd = f"gtdbtk classify_wf --genome_dir {bins_dir} --out_dir {output_dir} --cpus {threads} --extension fa"
        try:
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            error_msg = f"GTDB-TK classification failed: {e.stderr if e.stderr else str(e)}"
            raise RuntimeError(error_msg)

    def _run_bakta(self, bins_dir: str, output_dir: str, threads: int) -> None:
        """Run Bakta annotation"""
        # Verify bins directory exists
        if not os.path.exists(bins_dir):
            raise RuntimeError(f"Bins directory not found: {bins_dir}")
            
        # Verify Bakta database exists
        db_path = self._check_bakta_db()
        if not db_path:
            raise RuntimeError("Bakta database not found")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Run Bakta with error capture
        for bin_file in glob.glob(os.path.join(bins_dir, "*.fa")):
            bin_name = os.path.basename(bin_file)
            bin_output = os.path.join(output_dir, bin_name)
            os.makedirs(bin_output, exist_ok=True)
            
            cmd = f"bakta --db {db_path} --output {bin_output} --prefix {bin_name} --threads {threads} {bin_file}"
            try:
                result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                error_msg = f"Bakta annotation failed: {e.stderr if e.stderr else str(e)}"
                raise RuntimeError(error_msg)

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
        if not self.django_db:
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
                self.django_db.upload_mag(
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
                    self.django_db.upload_kegg_map(
                        mag_name=f"{self.sample_name}_{bin_prefix}",
                        **kegg_map
                    )

        except Exception as e:
            print(f"Error uploading results: {str(e)}")
            raise

    def __del__(self):
        """Clean up temp directory when pipeline is done"""
        try:
            if hasattr(self, 'temp_dir') and os.path.exists(self.temp_dir):
                shutil.rmtree(self.temp_dir)
        except Exception as e:
            print(f"Warning: Could not remove temp directory {self.temp_dir}: {e}")
        finally:
            # Make sure to mark pipeline as failed if any step failed
            if hasattr(self, 'state_tracker'):
                self.state_tracker.finalize_pipeline_state()