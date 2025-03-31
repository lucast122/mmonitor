#!/usr/bin/env python3
import os
import json
import shutil
import subprocess
import sys
import glob
from create_icns import create_icns

def clean_build_dirs():
    """Clean up build directories before building."""
    dirs_to_clean = ['build', 'dist']
    for dir_name in dirs_to_clean:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
            print(f"Cleaned {dir_name} directory")

def build_mac():
    """Build macOS application."""
    # Clean build directories
    clean_build_dirs()
    
    # Setup paths
    project_root = os.path.dirname(os.path.abspath(__file__))
    main_path = os.path.join(project_root, 'desktop', 'src', 'mmonitor', '__main__.py')
    hooks_path = os.path.join(project_root, 'hooks')
    runtime_hooks_path = os.path.join(project_root, 'runtime_hooks.py')
    theme_path = os.path.join(project_root, 'desktop', 'src', 'mmonitor', 'resources', 'grey_theme.json')
    src_dir = os.path.join(project_root, 'desktop', 'src')
    emu_db_path = os.path.join(project_root, 'desktop', 'src', 'resources', 'emu_db')
    
    # Only create icon if it doesn't exist
    icon_path = os.path.join(project_root, 'desktop', 'src', 'resources', 'icons', 'MMonitor.png')
    if not os.path.exists(icon_path):
        print("Creating application icon...")
        icon_path = create_icns()
    
    # Server paths and essential files
    server_root = os.path.join(project_root, 'server')
    server_essential_dirs = [
        'mmonitor',
        'media',
        'static',
        'staticfiles',
        'templates'
    ]
    
    server_essential_files = [
        'manage.py',
        'requirements.txt',
    ]
    
    # Create temporary resources directory for writable files
    temp_resources_dir = os.path.join(project_root, 'temp_resources')
    os.makedirs(temp_resources_dir, exist_ok=True)
    
    # Create empty pipeline_config.json if it doesn't exist
    config_file = os.path.join(temp_resources_dir, 'pipeline_config.json')
    if not os.path.exists(config_file):
        with open(config_file, 'w') as f:
            json.dump({}, f)
    
    # Create temporary server directory
    temp_server_dir = os.path.join(project_root, 'temp_server')
    if os.path.exists(temp_server_dir):
        shutil.rmtree(temp_server_dir)
    os.makedirs(temp_server_dir)
    
    # Create all server directories in one pass
    for dir_name in server_essential_dirs:
        dst = os.path.join(temp_server_dir, dir_name)
        src = os.path.join(server_root, dir_name)
        os.makedirs(dst, exist_ok=True)
        if os.path.exists(src):
            shutil.copytree(src, dst, dirs_exist_ok=True)
            print(f"Copied {dir_name} directory")
    
    # Copy essential server files
    for file_name in server_essential_files:
        src = os.path.join(server_root, file_name)
        dst = os.path.join(temp_server_dir, file_name)
        if os.path.exists(src):
            shutil.copy2(src, dst)
        else:
            open(dst, 'a').close()
    
    # Create empty db.sqlite3 file
    open(os.path.join(temp_server_dir, 'db.sqlite3'), 'w').close()
    
    # Create error logging script
    error_log_script = os.path.join(project_root, 'error_logging_hook.py')
    with open(error_log_script, 'w') as f:
        f.write('''
import os
import sys
import traceback
import logging
import datetime
from pathlib import Path

def setup_error_logging():
    """
    Set up error logging for the application.
    This redirects stdout and stderr to a log file.
    """
    # Determine app location and create log directory
    try:
        # Get the app's location (works in both dev and bundled app)
        if getattr(sys, 'frozen', False):
            # We're in a PyInstaller bundle
            app_dir = os.path.dirname(sys.executable)
            if sys.platform == 'darwin':  # macOS
                # Navigate up from .app/Contents/MacOS/
                app_dir = os.path.dirname(os.path.dirname(os.path.dirname(app_dir)))
        else:
            # We're running in a normal Python process
            app_dir = os.path.dirname(os.path.abspath(__file__))
            
        # Create a logs directory next to the application
        log_dir = os.path.join(os.path.dirname(app_dir), 'logs')
        os.makedirs(log_dir, exist_ok=True)
            
        # Set up logging with timestamp
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = os.path.join(log_dir, f"mmonitor_error_{timestamp}.log")
            
        # Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            filename=log_file,
            filemode='w'
        )
            
        # Add console handler to also print to stdout if we're not in a bundled app
        if not getattr(sys, 'frozen', False):
            console = logging.StreamHandler()
            console.setLevel(logging.DEBUG)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            console.setFormatter(formatter)
            logging.getLogger('').addHandler(console)
            
        # Redirect stdout and stderr to the log file
        sys.stdout = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
        sys.stderr = StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)
            
        # Log system information
        logging.info(f"Starting MMonitor application")
        logging.info(f"Python version: {sys.version}")
        logging.info(f"Platform: {sys.platform}")
        logging.info(f"Application directory: {app_dir}")
        logging.info(f"Log file: {log_file}")
            
        # Create an alternative log file (more accessible location)
        home_dir = str(Path.home())
        home_log_dir = os.path.join(home_dir, '.mmonitor', 'logs')
        os.makedirs(home_log_dir, exist_ok=True)
        home_log_file = os.path.join(home_log_dir, f"mmonitor_error_{timestamp}.log")
            
        # Link to the home directory log file
        with open(home_log_file, 'w') as f:
            f.write(f"MMonitor started at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
            f.write(f"Main log file is at: {log_file}\\n\\n")
            f.write("Application startup log:\\n\\n")
            
        # Set up exception handler
        original_excepthook = sys.excepthook
            
        def exception_handler(exc_type, exc_value, exc_tb):
            # Log the exception
            logging.error("Uncaught exception:", exc_info=(exc_type, exc_value, exc_tb))
                
            # Also log to home directory for easier access
            with open(home_log_file, 'a') as f:
                f.write("\\n\\nUncaught exception:\\n")
                f.write(''.join(traceback.format_exception(exc_type, exc_value, exc_tb)))
                
            # Call the original exception handler
            original_excepthook(exc_type, exc_value, exc_tb)
            
        sys.excepthook = exception_handler
            
        return log_file
            
    except Exception as e:
        # If we can't set up logging, at least try to write to a file in the home directory
        try:
            home_dir = str(Path.home())
            fallback_log = os.path.join(home_dir, 'mmonitor_error.log')
            with open(fallback_log, 'w') as f:
                f.write(f"Error setting up logging: {str(e)}\\n")
                f.write(traceback.format_exc())
            return fallback_log
        except:
            # Last resort - can't do anything more
            return None

class StreamToLogger:
    """
    File-like object that logs everything written to it.
    """
    def __init__(self, logger, log_level):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        temp_linebuf = self.linebuf + buf
        self.linebuf = ''
        for line in temp_linebuf.splitlines(True):
            if line.endswith(('\\n', '\\r')):
                self.logger.log(self.log_level, line.rstrip())
                self.linebuf = ''
            else:
                self.linebuf = line

    def flush(self):
        if self.linebuf:
            self.logger.log(self.log_level, self.linebuf.rstrip())
            self.linebuf = ''

# Setup error logging when this module is imported
log_file = setup_error_logging()
print(f"Log file: {log_file}")
''')
    
    # Create base command list with optimized settings
    cmd = [
        'pyinstaller',
        '--name=MMonitor',
        '--windowed',  # Create a macOS .app bundle
        '--clean',
        '--noconfirm',
    ]
    
    # Add icon and app settings
    if os.path.exists(icon_path):
        cmd.extend([
            f'--icon={icon_path}',
            '--osx-bundle-identifier=com.mmonitor.app',
        ])
    
    # Target current architecture
    cmd.append('--target-arch=x86_64' if os.uname().machine == 'x86_64' else '--target-arch=arm64')
    
    # Add core imports
    cmd.extend([
        '--hidden-import=sqlite3',
        '--hidden-import=numpy',
        '--hidden-import=customtkinter',
        '--hidden-import=PIL',
        '--hidden-import=PIL._webp',
        '--hidden-import=PIL._imaging',
        '--hidden-import=CTkMessagebox',
        '--hidden-import=mmonitor.paths',
    ])
    
    # Create temporary resources directory with only necessary files
    temp_build_resources = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build', 'app_resources')
    if os.path.exists(temp_build_resources):
        shutil.rmtree(temp_build_resources)
    os.makedirs(temp_build_resources, exist_ok=True)
    
    # Define and create necessary resource directories
    resource_dirs = {
        'icons': ['*.png', '*.icns'],
        'images': ['*.png', '*.jpg', '*.jpeg', '*.gif'],
        'themes': ['*.json'],
        'config': ['*.json', '*.yml', '*.yaml'],
    }
    
    # Copy only necessary resources with specific file patterns
    src_resources = os.path.join(project_root, "desktop", "src", "resources")
    print(f"Copying resources from: {src_resources}")
    print(f"To: {temp_build_resources}")
    
    for dir_name, patterns in resource_dirs.items():
        src_dir = os.path.join(src_resources, dir_name)
        dst_dir = os.path.join(temp_build_resources, dir_name)
        
        if os.path.exists(src_dir):
            os.makedirs(dst_dir, exist_ok=True)
            print(f"Processing {dir_name}...")
            
            # Copy files matching patterns
            for pattern in patterns:
                for src_file in glob.glob(os.path.join(src_dir, "**", pattern), recursive=True):
                    rel_path = os.path.relpath(src_file, src_dir)
                    dst_file = os.path.join(dst_dir, rel_path)
                    os.makedirs(os.path.dirname(dst_file), exist_ok=True)
                    shutil.copy2(src_file, dst_file)
                    print(f"  Copied: {rel_path}")
    
    # Copy specific required files
    required_files = [
        (os.path.abspath(emu_db_path), 'resources/emu_db'),
        (os.path.abspath(theme_path), 'mmonitor/resources'),
    ]
    
    # Add data files with absolute paths and correct destinations
    for src, dst in required_files:
        if os.path.exists(src):
            cmd.extend([f'--add-data={src}:{dst}'])
            print(f"Adding data file: {src} -> {dst}")
    
    # Add the processed resources
    if os.path.exists(temp_build_resources):
        cmd.extend([f'--add-data={temp_build_resources}:resources'])
        print(f"Adding resources directory: {temp_build_resources}")
    
    # Exclude large database directories and unnecessary files
    cmd.extend([
        '--exclude-module=tensorflow',
        '--exclude-module=torch',
        '--exclude-module=custom_centrifuger_db',
        '--exclude-module=ncbi_build',
    ])
    
    # Add specific file/directory exclusions
    cmd.extend([
        '--exclude=*.gz',
        '--exclude=*.fna',
        '--exclude=*.fasta',
        '--exclude=*.fa',
        '--exclude=*.db',
        '--exclude=custom_centrifuger_db',
        '--exclude=ncbi_build*',
        '--exclude=__pycache__',
        '--exclude=*.pyc',
        '--exclude=*.pyo',
        '--exclude=*.pyd',
        '--exclude=*.so',
    ])
    
    # Add paths and hooks
    cmd.extend([
        f'--paths={src_dir}',
        f'--paths={server_root}',
        f'--additional-hooks-dir={hooks_path}',
        f'--runtime-hook={os.path.join(hooks_path, "exclude_files.py")}',
        f'--runtime-hook={runtime_hooks_path}',
        f'--runtime-hook={error_log_script}',  # Add our error logging hook
    ])
    
    # Add required package collections
    cmd.extend([
        '--collect-all=customtkinter',
        '--collect-all=CTkMessagebox',
        '--collect-all=PIL',
    ])
    
    # Add main script
    cmd.append(main_path)
    
    # Run PyInstaller with optimized environment
    env = os.environ.copy()
    env['PYTHONOPTIMIZE'] = '1'
    env['PYTHONUTF8'] = '1'
    
    # Clean previous build
    for dir_name in ['build', 'dist']:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
            print(f"Cleaned {dir_name} directory")
    
    print("Starting build process...")
    try:
        # First try to import all required modules to catch any import errors
        import_check = subprocess.run(
            [sys.executable, '-c', 'import customtkinter, PIL, numpy, sqlite3, CTkMessagebox'],
            env=env,
            capture_output=True,
            text=True
        )
        if import_check.returncode != 0:
            print("Error: Failed to import required modules:")
            print(import_check.stderr)
            sys.exit(1)
            
        # Run the build
        build_process = subprocess.run(cmd, env=env, capture_output=True, text=True)
        
        if build_process.returncode != 0:
            print("Error during build process:")
            print(build_process.stdout)
            print(build_process.stderr)
            sys.exit(1)
        
        app_path = os.path.join('dist', 'MMonitor.app')
        if os.path.exists(app_path):
            print(f"\nApplication bundle created at: {app_path}")
            print("You can now:")
            print("1. Drag MMonitor.app to your Applications folder")
            print("2. Double-click to run the application")
            print("\nNote: Large database files are excluded from the bundle.")
            print("They should be downloaded or copied separately to the appropriate location.")
            print("\nDEBUG INFO: Error logs will be created in:")
            print(f"1. ~/logs directory next to the app")
            print(f"2. ~/.mmonitor/logs/ for easier access")
            
        print("Build completed successfully!")
    except Exception as e:
        print(f"Error during build process: {str(e)}")
        sys.exit(1)
    finally:
        # Clean up temporary directories
        if os.path.exists(temp_build_resources):
            shutil.rmtree(temp_build_resources)
        for temp_dir in [temp_server_dir, temp_resources_dir]:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

if __name__ == '__main__':
    build_mac()
