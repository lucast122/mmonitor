import os
import subprocess
import logging
import tempfile
import gzip
from pathlib import Path

logger = logging.getLogger(__name__)

class FiltlongRunner:
    def __init__(self):
        """Initialize the Filtlong runner"""
        pass
        
    def run_filtlong(self, input_file, output_file, min_length=1000, min_mean_q=7, 
                    target_bases=None, keep_percent=90, trim=False, split=False, 
                    max_length=None, threads=1):
        """
        Run Filtlong on input FASTQ file
        
        Args:
            input_file: Path to input FASTQ file
            output_file: Path to output filtered FASTQ file
            min_length: Minimum read length to keep (default: 1000)
            min_mean_q: Minimum read quality to keep (default: 7)
            target_bases: Target number of bases to retain (default: None)
            keep_percent: Keep only this percentage of best reads (default: 90)
            trim: Trim adapters from reads (default: False)
            split: Split reads at adapters (default: False)
            max_length: Maximum read length to keep (default: None)
            threads: Number of threads to use (default: 1)
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Ensure output directory exists
            os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
            
            # Create a temporary file for uncompressed output
            with tempfile.NamedTemporaryFile(suffix='.fastq', delete=False) as temp_out:
                temp_output = temp_out.name
            
            # Check if filtlong is available in PATH
            filtlong_cmd = "filtlong"
            try:
                # Try to run filtlong to check if it's in PATH
                subprocess.run([filtlong_cmd, "--version"], capture_output=True, text=True, check=False)
            except FileNotFoundError:
                # If not in PATH, try to use the one from lib directory
                script_dir = os.path.dirname(os.path.abspath(__file__))
                base_dir = os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))
                filtlong_path = os.path.join(base_dir, "lib", "Filtlong", "bin", "filtlong")
                
                if os.path.exists(filtlong_path):
                    filtlong_cmd = filtlong_path
                else:
                    raise FileNotFoundError(f"Filtlong not found in PATH or at {filtlong_path}")
            
            # Build the command
            cmd = [filtlong_cmd]
            
            # Add parameters
            cmd.extend(["--min_length", str(min_length)])
            cmd.extend(["--min_mean_q", str(min_mean_q)])
            
            if target_bases:
                cmd.extend(["--target_bases", str(target_bases)])
                
            if keep_percent and keep_percent < 100:
                cmd.extend(["--keep_percent", str(keep_percent)])
                
            if trim:
                cmd.append("--trim")
                
            if split:
                cmd.append("--split")
                
            if max_length:
                cmd.extend(["--max_length", str(max_length)])
            
            # Handle input file - decompress if needed
            if input_file.endswith('.gz'):
                with tempfile.NamedTemporaryFile(suffix='.fastq', delete=False) as temp_in:
                    with gzip.open(input_file, 'rb') as gz_in:
                        temp_in.write(gz_in.read())
                    input_path = temp_in.name
            else:
                input_path = input_file
            
            # Add input file to command
            cmd.append(input_path)
            
            # Run the command and redirect output to temporary file
            logger.info(f"Running Filtlong: {' '.join(cmd)}")
            with open(temp_output, 'w') as out_file:
                result = subprocess.run(cmd, stdout=out_file, stderr=subprocess.PIPE, text=True)
            
            # Clean up temporary input file if we created one
            if input_file.endswith('.gz'):
                os.unlink(input_path)
            
            if result.returncode != 0:
                logger.error(f"Filtlong error: {result.stderr}")
                os.unlink(temp_output)
                return False
            
            # Compress the output
            logger.info("Compressing filtered output...")
            with open(temp_output, 'rb') as f_in:
                with gzip.open(output_file, 'wb', compresslevel=6) as f_out:
                    while True:
                        chunk = f_in.read(1024*1024)  # 1MB chunks
                        if not chunk:
                            break
                        f_out.write(chunk)
            
            # Clean up temporary output file
            os.unlink(temp_output)
            
            logger.info(f"Filtlong completed successfully. Output: {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error running Filtlong: {e}")
            # Clean up temporary files if they exist
            if 'temp_output' in locals() and os.path.exists(temp_output):
                os.unlink(temp_output)
            if 'input_path' in locals() and input_file.endswith('.gz') and os.path.exists(input_path):
                os.unlink(input_path)
            return False

    def filter_fastq_files(self, input_files, output_dir, sample_name, **kwargs):
        """
        Filter multiple FASTQ files
        
        Args:
            input_files: List of input FASTQ files
            output_dir: Directory to output filtered files
            sample_name: Sample name for output file naming
            **kwargs: Additional parameters to pass to run_filtlong
            
        Returns:
            List of filtered FASTQ file paths or None if any step fails
        """
        if not input_files:
            logger.error("No input files provided")
            return None
            
        filtered_files = []
        os.makedirs(output_dir, exist_ok=True)
        
        for idx, input_file in enumerate(input_files):
            # Generate output filename
            input_basename = os.path.basename(input_file)
            if input_basename.endswith('.gz'):
                input_basename = input_basename[:-3]
            filtered_filename = f"{sample_name}_filtered_{idx}_{input_basename}.gz"
            output_file = os.path.join(output_dir, filtered_filename)
            
            # Run filtlong
            success = self.run_filtlong(input_file, output_file, **kwargs)
            if not success:
                logger.error(f"Failed to filter {input_file}")
                return None
                
            filtered_files.append(output_file)
            
        return filtered_files 