import os
import sys
import subprocess
import signal
import time
import shutil
import argparse
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('centrifuger_wrapper.log')
    ]
)
logger = logging.getLogger('centrifuger_wrapper')

def run_centrifuger(args):
    """Run centrifuger-realtime with fallback to regular centrifuger mode if realtime fails"""
    centrifuger_path = args.centrifuger_path
    index_prefix = args.index_prefix
    threads = args.threads
    topk = args.topk
    watch_dir = args.watch_dir
    min_files = args.min_files
    output_file = args.output
    
    # Ensure the watch directory exists
    os.makedirs(watch_dir, exist_ok=True)
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)
    
    # First try realtime mode
    realtime_cmd = [
        centrifuger_path,
        '-x', index_prefix,
        '-t', str(threads),
        '-k', str(topk),
        '--realtime',
        '--watch-dir', watch_dir,
        '--min-files', str(min_files),
        '--output', output_file
    ]
    
    logger.info(f"Attempting to run in realtime mode: {' '.join(realtime_cmd)}")
    
    try:
        result = subprocess.run(
            realtime_cmd, 
            capture_output=True, 
            text=True, 
            check=False
        )
        
        # Check if command failed with segfault
        if result.returncode == -signal.SIGSEGV:
            logger.warning("Realtime mode crashed with segmentation fault. Falling back to batch mode.")
            run_batch_mode(args)
        elif result.returncode != 0:
            logger.error(f"Realtime mode failed with error code {result.returncode}")
            logger.error(f"STDOUT: {result.stdout}")
            logger.error(f"STDERR: {result.stderr}")
            run_batch_mode(args)
        else:
            logger.info("Realtime mode completed successfully.")
            logger.info(f"STDOUT: {result.stdout}")
            logger.info(f"STDERR: {result.stderr}")
    
    except Exception as e:
        logger.error(f"Exception running realtime mode: {str(e)}")
        run_batch_mode(args)

def run_batch_mode(args):
    """Run centrifuger in batch mode as a fallback"""
    logger.info("Starting batch mode as fallback")
    
    watch_dir = args.watch_dir
    centrifuger_path = args.centrifuger_path.replace('-realtime', '')  # Use regular centrifuger
    index_prefix = args.index_prefix
    threads = args.threads
    topk = args.topk
    output_file = args.output
    
    # Find all FASTQ files in the watch directory
    fastq_files = []
    for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
        fastq_files.extend(list(Path(watch_dir).glob(f'*{ext}')))
    
    if not fastq_files:
        logger.warning(f"No FASTQ files found in {watch_dir}")
        return
    
    logger.info(f"Found {len(fastq_files)} FASTQ files")
    
    # Run centrifuger on all files
    # Check if we need to use --out or -o based on the command availability
    output_cmd = '--out'
    
    # Build command
    batch_cmd = [
        centrifuger_path,
        '-x', index_prefix,
        '-t', str(threads),
        '-k', str(topk)
    ]
    
    # Add output file parameter
    batch_cmd.extend([output_cmd, output_file])
    
    # Add all fastq files
    for fastq in fastq_files:
        batch_cmd.append('-U')
        batch_cmd.append(str(fastq))
    
    logger.info(f"Running batch command: {' '.join(batch_cmd)}")
    
    try:
        # Run the command and capture output
        process = subprocess.Popen(
            batch_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1
        )
        
        # Process stdout
        for line in process.stdout:
            logger.info(line.strip())
            
        # Process stderr
        for line in process.stderr:
            logger.error(line.strip())
            
        # Wait for process to complete
        returncode = process.wait()
        
        if returncode != 0:
            logger.error(f"Batch processing failed with return code {returncode}")
            # Try again with -o instead of --out if that might be the issue
            if output_cmd == '--out':
                logger.info("Retrying with -o instead of --out...")
                # Update the command with -o
                batch_cmd = [cmd if cmd != '--out' else '-o' for cmd in batch_cmd]
                logger.info(f"Running modified batch command: {' '.join(batch_cmd)}")
                result = subprocess.run(batch_cmd, capture_output=True, text=True, check=True)
                logger.info("Batch processing completed successfully on second attempt")
                if result.stdout:
                    logger.info(result.stdout)
                if result.stderr:
                    logger.error(result.stderr)
        else:
            logger.info("Batch processing completed successfully")
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Batch processing failed with error code {e.returncode}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
    except Exception as e:
        logger.error(f"Exception during batch processing: {str(e)}")
        logger.error(f"Command was: {' '.join(batch_cmd)}")

def monitor_directory(args):
    """Watch directory and process new files periodically"""
    watch_dir = args.watch_dir
    check_interval = args.check_interval
    
    logger.info(f"Starting directory monitoring of {watch_dir} (interval: {check_interval}s)")
    
    processed_files = set()
    
    while True:
        # Find all FASTQ files
        fastq_files = []
        for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
            fastq_files.extend([str(p) for p in Path(watch_dir).glob(f'*{ext}')])
        
        # Find new files
        new_files = [f for f in fastq_files if f not in processed_files]
        
        if new_files:
            logger.info(f"Found {len(new_files)} new files to process")
            
            # Create a temporary directory to store new files for batch processing
            temp_dir = os.path.join(watch_dir, "temp_batch")
            os.makedirs(temp_dir, exist_ok=True)
            
            # Copy new files to temp directory
            for file in new_files:
                shutil.copy(file, temp_dir)
            
            # Process the batch
            batch_args = argparse.Namespace(
                centrifuger_path=args.centrifuger_path.replace('-realtime', ''),
                index_prefix=args.index_prefix,
                threads=args.threads,
                topk=args.topk,
                watch_dir=temp_dir,
                min_files=args.min_files,
                output=args.output,
                check_interval=args.check_interval
            )
            
            try:
                run_batch_mode(batch_args)
                
                # Add processed files to our set
                processed_files.update(new_files)
                
                # Clean up temp directory
                shutil.rmtree(temp_dir)
                
            except Exception as e:
                logger.error(f"Error processing batch: {str(e)}")
        
        # Wait for the next check
        time.sleep(check_interval)

def main():
    parser = argparse.ArgumentParser(description='Centrifuger wrapper with realtime/batch fallback')
    parser.add_argument('-x', '--index-prefix', required=True, help='Index prefix path')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('-k', '--topk', type=int, default=1, help='Report up to k hits')
    parser.add_argument('--watch-dir', required=True, help='Directory to watch for new files')
    parser.add_argument('--min-files', type=int, default=5, help='Minimum files to process')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--centrifuger-path', default='centrifuger-realtime', help='Path to centrifuger-realtime binary')
    parser.add_argument('--mode', choices=['auto', 'realtime', 'batch', 'monitor'], default='auto', 
                        help='Mode to run (auto: try realtime then batch, realtime: only realtime, batch: only batch, monitor: poll directory)')
    parser.add_argument('--check-interval', type=int, default=30, help='Interval (seconds) to check for new files in monitor mode')
    
    args = parser.parse_args()
    
    logger.info(f"Starting centrifuger wrapper in {args.mode} mode")
    
    if args.mode == 'auto':
        run_centrifuger(args)
    elif args.mode == 'realtime':
        # Run only in realtime mode
        realtime_cmd = [
            args.centrifuger_path,
            '-x', args.index_prefix,
            '-t', str(args.threads),
            '-k', str(args.topk),
            '--realtime',
            '--watch-dir', args.watch_dir,
            '--min-files', str(args.min_files),
            '--output', args.output
        ]
        subprocess.run(realtime_cmd, check=True)
    elif args.mode == 'batch':
        run_batch_mode(args)
    elif args.mode == 'monitor':
        monitor_directory(args)

if __name__ == "__main__":
    main() 