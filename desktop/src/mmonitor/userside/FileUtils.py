import os
import gzip
import shutil
import tempfile
import subprocess
from typing import List
from pathlib import Path

def concatenate_fastq_files(input_files: List[str], output_file: str, threads: int = 1) -> None:
    """Concatenate multiple FASTQ files into a single file using Rust implementation
    with Python fallback if Rust fails. If only one input file is provided, creates
    a symlink instead of concatenating.
    
    Args:
        input_files: List of input FASTQ files
        output_file: Path to output file
        threads: Number of threads to use
    """
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # If only one input file, create a symlink instead of concatenating
    if len(input_files) == 1:
        input_file = input_files[0]
        # Remove output file if it exists
        if os.path.exists(output_file):
            os.remove(output_file)
        # Create symlink
        os.symlink(os.path.abspath(input_file), output_file)
        print(f"Created symlink from {input_file} to {output_file}")
        return
    
    # Try Rust implementation first
    try:
        # Get path to Rust binary
        current_dir = os.path.dirname(os.path.abspath(__file__))
        rust_binary = os.path.join(current_dir, "..", "..", "lib", "fastq_tools", "target", "release", "fastq_concat")
        if os.name == 'nt':  # Windows
            rust_binary += '.exe'
            
        # Ensure binary exists and is executable
        if not os.path.exists(rust_binary):
            raise FileNotFoundError(f"Rust binary not found: {rust_binary}")
        
        if os.name != 'nt':  # Not Windows
            os.chmod(rust_binary, 0o755)  # Make executable
            
        # Run Rust concatenation
        print(f"Running fast concatenation with {threads} threads...")
        cmd = [
            rust_binary,
            "--threads", str(threads),
            "--output", output_file,
            *input_files
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Rust concatenation failed: {result.stderr}")
            
        print("Successfully finished Rust concatenation")
        
    except Exception as e:
        print(f"Rust concatenation failed ({str(e)}), falling back to Python implementation...")
        
        # Python fallback implementation
        all_gzipped = all(f.endswith('.gz') for f in input_files)
        out_mode = 'wb' if output_file.endswith('.gz') else 'w'
        out_open = gzip.open if output_file.endswith('.gz') else open
        
        with out_open(output_file, out_mode) as outfile:
            for input_file in input_files:
                # Check if input file exists
                if not os.path.exists(input_file):
                    raise FileNotFoundError(f"Input file not found: {input_file}")
                
                # Open input file with appropriate mode
                in_mode = 'rb' if input_file.endswith('.gz') else 'r'
                in_open = gzip.open if input_file.endswith('.gz') else open
                
                with in_open(input_file, in_mode) as infile:
                    shutil.copyfileobj(infile, outfile)
        
        print("Successfully finished Python concatenation")
