
# Runtime hook to set up environment variables
import os
import sys
import signal
import atexit
import traceback
import subprocess
from pathlib import Path

# Global flag to track whether we're already running
_is_running = False

# Patch subprocess.Popen to handle centrifuger calls
original_popen = subprocess.Popen

def patched_popen(cmd, *args, **kwargs):
    # Check if this is a centrifuger call
    if isinstance(cmd, (list, tuple)) and len(cmd) > 0 and 'centrifuger' in str(cmd[0]):
        print(f"Intercepted centrifuger call: {cmd}")
        
        # Check if we're in a frozen app
        if getattr(sys, 'frozen', False):
            # Get the application resources directory
            if sys.platform == 'darwin':
                app_dir = os.path.abspath(os.path.join(os.path.dirname(sys.executable), '..', 'Resources'))
                wrapper_path = os.path.join(app_dir, 'lib', 'centrifuger-realtime', 'centrifuger-wrapper.sh')
                
                if os.path.exists(wrapper_path):
                    print(f"Using centrifuger wrapper: {wrapper_path}")
                    
                    # Replace the centrifuger command with the wrapper
                    if isinstance(cmd, list):
                        cmd[0] = wrapper_path
                    else:
                        cmd = [wrapper_path] + list(cmd[1:])
                    
                    print(f"Modified command: {cmd}")
    
    # Call the original Popen
    return original_popen(cmd, *args, **kwargs)

# Apply the patch
subprocess.Popen = patched_popen

# Set up proper signal handling
def handle_exit(signum, frame):
    print(f"Received signal {signum}, exiting cleanly")
    sys.exit(0)

# Register signal handlers for clean exit
signal.signal(signal.SIGINT, handle_exit)
signal.signal(signal.SIGTERM, handle_exit)

def cleanup():
    # Perform any cleanup needed
    print("Performing cleanup before exit")
    
    # Remove any temporary files
    try:
        import tempfile
        temp_dir = os.path.join(tempfile.gettempdir(), 'mmonitor')
        if os.path.exists(temp_dir):
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)
    except Exception as e:
        print(f"Error during cleanup: {e}")

atexit.register(cleanup)

# Exception hook to log unhandled exceptions
def exception_hook(exctype, value, tb):
    # Log the exception
    error_msg = ''.join(traceback.format_exception(exctype, value, tb))
    print(f"Unhandled exception: {error_msg}")
    
    # Call the original exception hook
    sys.__excepthook__(exctype, value, tb)

sys.excepthook = exception_hook

# Function to ensure centrifuger can be found
def check_centrifuger():
    # Check if centrifuger is in PATH
    import shutil
    centrifuger_path = shutil.which('centrifuger')
    
    if centrifuger_path:
        print(f"Found centrifuger in PATH: {centrifuger_path}")
        return True
    
    # Check if we're in a frozen app
    if getattr(sys, 'frozen', False):
        # Get the application resources directory
        if sys.platform == 'darwin':
            app_dir = os.path.abspath(os.path.join(os.path.dirname(sys.executable), '..', 'Resources'))
            bundled_path = os.path.join(app_dir, 'lib', 'centrifuger-realtime', 'centrifuger')
            
            if os.path.exists(bundled_path):
                print(f"Found bundled centrifuger: {bundled_path}")
                
                # Add the directory to PATH
                centrifuger_dir = os.path.dirname(bundled_path)
                os.environ['PATH'] = centrifuger_dir + os.pathsep + os.environ.get('PATH', '')
                print(f"Added {centrifuger_dir} to PATH")
                
                # Make it executable
                try:
                    import stat
                    os.chmod(bundled_path, os.stat(bundled_path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
                    print(f"Made centrifuger executable")
                except Exception as e:
                    print(f"Error setting executable permission: {e}")
                
                return True
    
    print("Warning: Centrifuger not found in PATH or bundle")
    return False

def setup_environment():
    global _is_running
    
    # Check if we're already running to prevent multiple instances
    if _is_running:
        print("Warning: setup_environment called multiple times")
        return
    
    _is_running = True
    
    # Print startup information
    print("="*50)
    print(f"Starting MMonitor - Python {sys.version}")
    print(f"Platform: {sys.platform}")
    print(f"Executable: {sys.executable}")
    print(f"Args: {sys.argv}")
    print("="*50)
    
    # Get the application path
    if getattr(sys, 'frozen', False):
        # Running as bundled app
        app_path = os.path.dirname(sys.executable)
        print(f"Running as frozen application from {app_path}")
        if sys.platform == 'darwin':
            # For macOS bundles, adjust the path
            app_path = os.path.abspath(os.path.join(app_path, '..', 'Resources'))
            print(f"Adjusted macOS resource path: {app_path}")
    else:
        # Running as script
        app_path = os.path.dirname(os.path.abspath(__file__))
        print(f"Running as script from {app_path}")
    
    # Set environment variables
    os.environ['MMONITOR_APP_PATH'] = app_path
    os.environ['PYTHONIOENCODING'] = 'utf-8'
    os.environ['PYTHONUNBUFFERED'] = '1'
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    
    # Prevent multiprocessing from spawning extra processes
    os.environ['PYTHONOPTIMIZE'] = '1'
    
    # Create logs directory if it doesn't exist
    logs_dir = os.path.expanduser('~/.mmonitor/logs')
    os.makedirs(logs_dir, exist_ok=True)
    
    # Create a log file with timestamp
    from datetime import datetime
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(logs_dir, f'mmonitor_{timestamp}.log')
    
    # Configure logging
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logging.info("Environment setup complete")
    
    # Log Python path
    logging.info(f"Python path: {sys.path}")
    
    # Check centrifuger availability
    check_centrifuger()
    
    # Check if common modules are importable
    try:
        import tkinter
        logging.info("tkinter is available")
    except ImportError as e:
        logging.error(f"tkinter import error: {e}")
    
    try:
        import customtkinter
        logging.info(f"customtkinter version: {customtkinter.__version__}")
    except ImportError as e:
        logging.error(f"customtkinter import error: {e}")
    except Exception as e:
        logging.error(f"customtkinter error: {e}")

setup_environment()
