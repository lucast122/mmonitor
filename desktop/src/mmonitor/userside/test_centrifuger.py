#!/usr/bin/env python3
import os
import sys
import argparse
from pathlib import Path
import logging

# Add the current directory to the path to make local imports work
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

# Import from the current directory
from centrifuger_handler import CentrifugerHandler

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('test_centrifuger.log')
    ]
)
logger = logging.getLogger('test_centrifuger')

def centrifuger_callback(success, message, details=None):
    """Callback function for CentrifugerHandler"""
    if success:
        logger.info(f"Centrifuger completed successfully: {message}")
    else:
        logger.error(f"Centrifuger failed: {message}")
        
    if details:
        logger.info(f"Return code: {details.get('returncode')}")
        
        # Log a sample of stdout (first 10 lines)
        stdout = details.get('stdout', '').split('\n')
        if stdout:
            logger.info("Sample of stdout (first 10 lines):")
            for i, line in enumerate(stdout[:10]):
                logger.info(f"  {line}")
            if len(stdout) > 10:
                logger.info(f"  ... and {len(stdout) - 10} more lines")
        
        # Log full stderr
        stderr = details.get('stderr', '').split('\n')
        if stderr:
            logger.info("Full stderr:")
            for line in stderr:
                if line.strip():
                    logger.info(f"  {line}")

def main():
    parser = argparse.ArgumentParser(description='Test the centrifuger wrapper')
    parser.add_argument('--index-prefix', required=True, help='Centrifuger index prefix')
    parser.add_argument('--watch-dir', required=True, help='Directory to watch for FASTQ files')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--min-files', type=int, default=5, help='Minimum files for analysis')
    parser.add_argument('--mode', choices=['auto', 'realtime', 'batch'], default='auto', 
                        help='Mode to run (auto, realtime, batch)')
    
    args = parser.parse_args()
    
    # Create handler
    handler = CentrifugerHandler()
    
    # Prepare parameters
    params = {
        'centrifuge_db': args.index_prefix,
        'threads': args.threads,
        'watch_dir': args.watch_dir,
        'output_file': args.output,
        'min_files_for_auto_analysis': args.min_files
    }
    
    # Determine if we want realtime mode
    realtime = args.mode != 'batch'
    
    # Log parameters
    logger.info(f"Running centrifuger with:")
    logger.info(f"  Index prefix: {params['centrifuge_db']}")
    logger.info(f"  Watch directory: {params['watch_dir']}")
    logger.info(f"  Output file: {params['output_file']}")
    logger.info(f"  Threads: {params['threads']}")
    logger.info(f"  Min files: {params['min_files_for_auto_analysis']}")
    logger.info(f"  Mode: {args.mode}")
    
    # Check paths
    if not os.path.exists(os.path.dirname(params['centrifuge_db'])):
        logger.error(f"Index directory does not exist: {os.path.dirname(params['centrifuge_db'])}")
        return 1
    
    # Ensure watch and output directories exist
    os.makedirs(args.watch_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Run centrifuger
    process_id = handler.run_centrifuger(params, realtime=realtime, callback=centrifuger_callback)
    
    if not process_id:
        logger.error("Failed to start centrifuger process")
        return 1
    
    logger.info(f"Centrifuger process started with ID: {process_id}")
    logger.info("Press Ctrl+C to terminate the process")
    
    try:
        # Wait for process to complete
        while handler.is_process_running(process_id):
            import time
            time.sleep(1)
        
        logger.info("Centrifuger process completed")
        return 0
        
    except KeyboardInterrupt:
        logger.info("Received keyboard interrupt, terminating process...")
        handler.terminate_process(process_id)
        return 130

if __name__ == "__main__":
    sys.exit(main()) 