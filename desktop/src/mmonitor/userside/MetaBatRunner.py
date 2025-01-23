import os
import platform
import subprocess
import shutil
from pathlib import Path

class MetaBatRunner:
    def __init__(self):
        self.metabat_path = shutil.which('metabat2')
        if not self.metabat_path:
            print("MetaBAT2 not found. Installing...")
            self.ensure_metabat_installed()
            self.metabat_path = shutil.which('metabat2')
            if not self.metabat_path:
                raise RuntimeError("Failed to find metabat2 after installation")

    def ensure_metabat_installed(self):
        """Ensure MetaBAT2 is installed via conda"""
        try:
            # Try installing with conda
            print("Installing MetaBAT2 via conda...")
            subprocess.run(
                ["conda", "install", "-c", "bioconda", "metabat2", "-y"],
                check=True,
                capture_output=True,
                text=True
            )
            print("MetaBAT2 installed successfully via conda")
        except subprocess.CalledProcessError as e:
            print(f"Error installing MetaBAT2 via conda: {e}")
            print("Output:", e.output)
            print("Error:", e.stderr)
            raise RuntimeError("Failed to install MetaBAT2. Please ensure conda is installed and try again.")

    def run_binning(self, assembly_file, output_dir, min_contig_length=1500, threads=1):
        """Run MetaBAT2 binning on assembly"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate depth file
        depth_file = os.path.join(output_dir, "depth.txt")
        cmd_depth = [
            self.metabat_path,
            "-i", assembly_file,
            "-o", os.path.join(output_dir, "bin"),
            "-m", str(min_contig_length),
            "-t", str(threads),
            "--unbinned"
        ]
        
        print(f"Running MetaBAT2: {' '.join(cmd_depth)}")
        try:
            subprocess.run(cmd_depth, check=True)
            return {
                "bins_dir": output_dir,
                "depth_file": depth_file
            }
        except subprocess.CalledProcessError as e:
            print(f"Error running MetaBAT2: {e}")
            raise
