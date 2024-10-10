import os
import time
import subprocess
import argparse
from datetime import datetime

def simulate_nanopore_reads(output_dir, interval=20, duration=300):
    """
    Simulate nanopore reads using NanoSim and save them to the specified directory.
    
    :param output_dir: Directory to save the simulated reads
    :param interval: Time interval between file generations (in seconds)
    :param duration: Total duration to run the simulation (in seconds)
    """
    start_time = time.time()
    file_count = 0
    
    while time.time() - start_time < duration:
        # Generate a unique filename based on the current timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"simulated_nanopore_reads_{timestamp}.fastq"
        output_file = os.path.join(output_dir, filename)
        
        # Run NanoSim to generate simulated reads
        command = [
            "nanosim",
            "-p", "ecoli_R9_2D",  # Pre-trained model for E. coli R9 2D reads
            "-n", "1000",  # Number of reads to generate
            "--perfect",  # Generate perfect reads without errors
            "-o", output_file
        ]
        
        try:
            subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"Generated: {filename}")
            file_count += 1
        except subprocess.CalledProcessError as e:
            print(f"Error generating reads: {e}")
        
        # Wait for the specified interval before generating the next file
        time.sleep(interval)
    
    print(f"Simulation complete. Generated {file_count} files.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate nanopore reads at regular intervals.")
    parser.add_argument("output_dir", help="Directory to save the simulated reads")
    parser.add_argument("--interval", type=int, default=20, help="Time interval between file generations (in seconds)")
    parser.add_argument("--duration", type=int, default=300, help="Total duration to run the simulation (in seconds)")
    
    args = parser.parse_args()
    
    simulate_nanopore_reads(args.output_dir, args.interval, args.duration)