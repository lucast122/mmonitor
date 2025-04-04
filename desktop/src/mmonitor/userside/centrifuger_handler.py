import os
import sys
import subprocess
import logging
import threading
from pathlib import Path
import time

# Define ROOT_DIR based on the file's location
current_file_path = os.path.abspath(__file__)
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(current_file_path)))
ROOT_DIR = os.path.dirname(os.path.dirname(src_dir))

class CentrifugerHandler:
    def __init__(self, log_file=None):
        self.logger = self._setup_logging(log_file)
        self.running_processes = {}
        
    def _setup_logging(self, log_file=None):
        logger = logging.getLogger('CentrifugerHandler')
        logger.setLevel(logging.INFO)
        
        # Create console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        
        # Create file handler if specified
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(logging.INFO)
            file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
        
        return logger
    
    def run_centrifuger(self, params, realtime=True, callback=None):
        """Run centrifuger using the wrapper script"""
        wrapper_script = os.path.join(os.path.dirname(__file__), "centrifuger_wrapper.py")
        
        # Make sure the wrapper script exists
        if not os.path.exists(wrapper_script):
            self.logger.error(f"Wrapper script not found: {wrapper_script}")
            if callback:
                callback(False, f"Wrapper script not found: {wrapper_script}")
            return False
        
        # Get parameters
        index_prefix = params.get("centrifuge_db", "")
        threads = params.get("threads", 1)
        watch_dir = params.get("watch_dir", "")
        output_file = params.get("output_file", "")
        min_files = params.get("min_files_for_auto_analysis", 5)
        
        # Validate parameters
        if not index_prefix or not os.path.exists(os.path.dirname(index_prefix)):
            self.logger.error(f"Invalid index prefix: {index_prefix}")
            if callback:
                callback(False, f"Invalid index prefix: {index_prefix}")
            return False
        
        if not watch_dir:
            self.logger.error("Watch directory not specified")
            if callback:
                callback(False, "Watch directory not specified")
            return False
        
        if not output_file:
            self.logger.error("Output file not specified")
            if callback:
                callback(False, "Output file not specified")
            return False
        
        # Ensure directories exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        os.makedirs(watch_dir, exist_ok=True)
        
        # Try to find centrifuger executables in several locations
        possible_paths = [
            # Current working directory
            os.path.join(os.getcwd(), "centrifuger-realtime" if realtime else "centrifuger"),
            # ROOT_DIR
            os.path.join(ROOT_DIR, "centrifuger-realtime" if realtime else "centrifuger"),
            # In PATH
            "centrifuger-realtime" if realtime else "centrifuger"
        ]
        
        # Find first existing path
        centrifuger_path = None
        for path in possible_paths:
            if os.path.exists(path) or self._is_executable_in_path(path):
                centrifuger_path = path
                break
        
        if not centrifuger_path:
            self.logger.error(f"Centrifuger executable not found. Tried: {possible_paths}")
            if callback:
                callback(False, f"Centrifuger executable not found. Tried: {possible_paths}")
            return False
        
        # Build command
        cmd = [
            sys.executable,
            wrapper_script,
            '--index-prefix', index_prefix,
            '--threads', str(threads),
            '--watch-dir', watch_dir,
            '--output', output_file,
            '--centrifuger-path', centrifuger_path,
            '--min-files', str(min_files),
            '--mode', 'auto' if realtime else 'batch',
            '--check-interval', '30'
        ]
        
        # Generate a unique ID for this process
        process_id = f"centrifuger_{int(time.time())}"
        
        # Log the command
        self.logger.info(f"Running centrifuger command (ID: {process_id}): {' '.join(cmd)}")
        
        # Run in a separate thread
        thread = threading.Thread(
            target=self._run_process,
            args=(cmd, process_id, callback)
        )
        thread.daemon = True
        thread.start()
        
        return process_id
    
    def _run_process(self, cmd, process_id, callback=None):
        """Run a subprocess and handle its output"""
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1
        )
        
        # Store process for potential termination
        self.running_processes[process_id] = process
        
        # Collect output
        stdout_lines = []
        stderr_lines = []
        
        # Process stdout
        for line in process.stdout:
            line = line.strip()
            if line:
                self.logger.info(f"[{process_id}] {line}")
                stdout_lines.append(line)
        
        # Process stderr
        for line in process.stderr:
            line = line.strip()
            if line:
                self.logger.error(f"[{process_id}] {line}")
                stderr_lines.append(line)
        
        # Wait for process to complete
        returncode = process.wait()
        
        # Remove from running processes
        if process_id in self.running_processes:
            del self.running_processes[process_id]
        
        # Determine success
        success = returncode == 0
        
        # Prepare message
        if success:
            message = "Centrifuger completed successfully"
        else:
            message = f"Centrifuger failed with return code {returncode}"
            if stderr_lines:
                message += f"\nError: {stderr_lines[-1]}"
        
        self.logger.info(f"[{process_id}] {message}")
        
        # Call callback if provided
        if callback:
            callback(success, message, {
                'stdout': '\n'.join(stdout_lines),
                'stderr': '\n'.join(stderr_lines),
                'returncode': returncode
            })
    
    def terminate_process(self, process_id):
        """Terminate a running process by ID"""
        if process_id in self.running_processes:
            process = self.running_processes[process_id]
            self.logger.info(f"Terminating process {process_id}")
            
            try:
                process.terminate()
                # Wait briefly to see if it terminates
                for _ in range(5):
                    if process.poll() is not None:
                        del self.running_processes[process_id]
                        self.logger.info(f"Process {process_id} terminated successfully")
                        return True
                    time.sleep(0.1)
                
                # If not terminated, kill it
                process.kill()
                process.wait()
                del self.running_processes[process_id]
                self.logger.info(f"Process {process_id} killed")
                return True
                
            except Exception as e:
                self.logger.error(f"Error terminating process {process_id}: {str(e)}")
                return False
        else:
            self.logger.warning(f"Process {process_id} not found")
            return False
    
    def is_process_running(self, process_id):
        """Check if a process is still running"""
        if process_id in self.running_processes:
            process = self.running_processes[process_id]
            return process.poll() is None
        return False
    
    def _is_executable_in_path(self, executable):
        """Check if an executable is in the PATH"""
        try:
            # Use 'which' on Unix/Mac or 'where' on Windows
            if sys.platform == "win32":
                result = subprocess.run(['where', executable], 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        check=False)
            else:
                result = subprocess.run(['which', executable], 
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE, 
                                        check=False)
            return result.returncode == 0
        except Exception:
            return False 