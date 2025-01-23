import shutil
import subprocess
import queue
import threading
import os
from tkinter import messagebox
import customtkinter as ctk
import sys

class InstallDialog(ctk.CTkToplevel):
    def __init__(self, parent, title, message):
        super().__init__(parent)
        
        self.title(title)
        self.result = False
        
        # Make dialog modal
        self.transient(parent)
        self.grab_set()
        
        # Center the dialog
        self.geometry(f"400x150")
        self.resizable(False, False)
        
        # Add message
        message_label = ctk.CTkLabel(self, text=message, wraplength=350)
        message_label.pack(pady=20)
        
        # Add buttons
        button_frame = ctk.CTkFrame(self)
        button_frame.pack(pady=10)
        
        yes_button = ctk.CTkButton(button_frame, text="Yes", command=self._on_yes)
        yes_button.pack(side="left", padx=10)
        
        no_button = ctk.CTkButton(button_frame, text="No", command=self._on_no)
        no_button.pack(side="left", padx=10)
        
        # Center the dialog on the parent window
        if parent:
            self.geometry("+%d+%d" % (
                parent.winfo_rootx() + parent.winfo_width()/2 - 200,
                parent.winfo_rooty() + parent.winfo_height()/2 - 75
            ))

    def _on_yes(self):
        self.result = True
        self.destroy()

    def _on_no(self):
        self.result = False
        self.destroy()

class ToolInstaller:
    TOOL_SPECS = {
        'minimap2': {'conda_package': 'minimap2', 'conda_channel': 'bioconda'},
        'flye': {'conda_package': 'flye', 'conda_channel': 'bioconda'},
        'medaka': {'env_name': 'medaka_env', 'python_version': '3.9'},
        'bakta': {'conda_package': 'bakta', 'conda_channel': 'bioconda'},
        'gtdbtk': {'conda_package': 'gtdbtk', 'conda_channel': 'bioconda'},
        'checkm2': {'pip_package': 'checkm2'},
        'metabat2': {'conda_package': 'metabat2', 'conda_channel': 'bioconda'}
    }

    _main_window = None
    _install_queue = queue.Queue()
    _response_queue = queue.Queue()
    _installation_lock = threading.Lock()
    _current_dialog = None

    @classmethod
    def set_main_window(cls, window):
        cls._main_window = window

    @classmethod
    def _setup_conda_env(cls):
        """Setup conda environment for all tools"""
        try:
            env_name = "mmonitor_tools"
            print(f"Setting up conda environment: {env_name}")
            
            # Create environment if it doesn't exist
            result = subprocess.run(
                ["conda", "env", "list"],
                capture_output=True,
                text=True
            )
            if env_name not in result.stdout:
                subprocess.run(
                    ["conda", "create", "-n", env_name, "python=3.9", "-y"],
                    check=True,
                    capture_output=True,
                    text=True
                )

            # Get conda activation script path
            conda_path = os.path.dirname(shutil.which('conda'))
            activate_script = os.path.join(conda_path, 'activate')
            
            return env_name, activate_script
        except Exception as e:
            print(f"Error setting up conda environment: {e}")
            return None, None

    @classmethod
    def _create_medaka_env(cls):
        """Create a dedicated Python 3.9 environment for Medaka"""
        env_name = "medaka_env"
        try:
            # Create environment with Python 3.9
            subprocess.run(['conda', 'create', '-n', env_name, 'python=3.9', '-y'], check=True)
            
            # Install Medaka using pip in the new environment
            subprocess.run(['conda', 'run', '-n', env_name, 'pip', 'install', 'medaka'], check=True)
            
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error creating Medaka environment: {str(e)}")
            return False

    @classmethod
    def _install_in_env(cls, tool_name, spec):
        """Install a tool in the conda environment"""
        env_name = "mmonitor_tools"
        
        try:
            print(f"Installing {tool_name}...")
            
            # Special handling for specific tools
            if tool_name == 'flye':
                # Install Flye with specific version that works with Python 3.9
                cmd = [
                    "conda", "install",
                    "-c", "bioconda",
                    "-c", "conda-forge",
                    "flye=2.9.2",
                    "-y"
                ]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                
            elif tool_name == 'checkm2':
                # First ensure we have the right Python environment
                subprocess.run(
                    ["conda", "install", "-c", "conda-forge", 
                     "python=3.9", "numpy", "scipy", "scikit-learn", "h5py",
                     "-y"],
                    check=True, capture_output=True, text=True
                )
                
                # Then install CheckM2 via pip
                subprocess.run(
                    ["pip", "install", "--no-deps", "checkm2"],
                    check=True, capture_output=True, text=True
                )
                
                # Install dependencies separately
                subprocess.run(
                    ["pip", "install", 
                     "numpy==1.24.3",
                     "scipy==1.10.1",
                     "scikit-learn==1.2.2",
                     "h5py==3.8.0"],
                    check=True, capture_output=True, text=True
                )
                
            elif tool_name == 'medaka':
                return cls._create_medaka_env()
                
            elif spec.get('conda_package'):
                # Default conda installation
                cmd = [
                    "conda", "install",
                    "-c", spec['conda_channel'],
                    "-c", "conda-forge",
                    spec['conda_package'],
                    "-y"
                ]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                
            elif spec.get('pip_package'):
                # Default pip installation
                cmd = ["pip", "install", spec['pip_package']]
                subprocess.run(cmd, check=True, capture_output=True, text=True)
            else:
                raise RuntimeError(f"No installation method available for {tool_name}")

            # Verify installation
            which_cmd = ["which", tool_name]
            result = subprocess.run(which_cmd, check=True, capture_output=True, text=True)
            tool_path = result.stdout.strip()
            
            if not tool_path:
                raise RuntimeError(f"Tool {tool_name} not found after installation")
                
            print(f"Successfully installed {tool_name} at {tool_path}")
            return tool_path

        except subprocess.CalledProcessError as e:
            print(f"Error installing {tool_name}:")
            if e.stdout:
                print("stdout:", e.stdout)
            if e.stderr:
                print("stderr:", e.stderr)
            return None
        except Exception as e:
            print(f"Error installing {tool_name}: {e}")
            return None

    @classmethod
    def _show_install_dialog(cls):
        """Show installation dialog and handle response"""
        if cls._current_dialog is not None:
            return

        try:
            tool_name, spec = cls._install_queue.get_nowait()
            dialog = InstallDialog(
                cls._main_window,
                "Install Required Tool",
                f"The tool '{tool_name}' is required but not installed. Would you like to install it?"
            )
            
            # Store dialog reference
            cls._current_dialog = dialog
            
            # Wait for dialog to close
            cls._main_window.wait_window(dialog)
            
            # Get result and clean up
            result = dialog.result
            cls._current_dialog = None
            cls._response_queue.put(result)
            
        except queue.Empty:
            pass

    @classmethod
    def check_tool(cls, tool_name):
        """Check if a tool is installed and offer to install it if not"""
        with cls._installation_lock:
            if tool_name not in cls.TOOL_SPECS:
                raise ValueError(f"Unknown tool: {tool_name}")

            # First check system PATH
            tool_path = shutil.which(tool_name)
            if tool_path:
                print(f"Found {tool_name} in system PATH: {tool_path}")
                return tool_path

            # Get tool specs
            tool_spec = cls.TOOL_SPECS.get(tool_name, {})
            
            # Special handling for Medaka
            if tool_name == 'medaka':
                env_name = tool_spec.get('env_name', 'medaka_env')
                env_bin = os.path.join(os.path.dirname(shutil.which('conda')), 'envs', env_name, 'bin', tool_name)
                
                if os.path.exists(env_bin):
                    print(f"Found {tool_name} in {env_name} environment: {env_bin}")
                    return env_bin
                    
                print(f"{tool_name} not found, creating dedicated Python 3.9 environment...")
                if cls._create_medaka_env():
                    if os.path.exists(env_bin):
                        print(f"Successfully installed {tool_name} in {env_name} environment: {env_bin}")
                        return env_bin
                return None

            # Then check conda environment
            env_name = "mmonitor_tools"
            conda_path = os.path.dirname(shutil.which('conda'))
            activate_script = os.path.join(conda_path, 'activate')
            
            try:
                which_cmd = [
                    "bash", "-c",
                    f"source {activate_script} {env_name} && which {tool_name}"
                ]
                result = subprocess.run(which_cmd, capture_output=True, text=True)
                if result.returncode == 0 and result.stdout.strip():
                    tool_path = result.stdout.strip()
                    print(f"Found {tool_name} in conda environment: {tool_path}")
                    return tool_path
            except:
                pass

            # Tool not found anywhere, ask user if they want to install it
            if cls._main_window is None:
                print(f"Warning: {tool_name} not found and no main window set for installation dialog")
                return None

            print(f"{tool_name} not found in system PATH or conda environment")

            # Put request in queue and wait for response
            cls._install_queue.put((tool_name, cls.TOOL_SPECS[tool_name]))
            
            # Schedule the dialog in the main thread
            cls._main_window.after(0, cls._show_install_dialog)
            
            # Wait for response
            try:
                result = cls._response_queue.get(timeout=60)
                if not result:
                    return None

                # Install the tool
                return cls._install_in_env(tool_name, cls.TOOL_SPECS[tool_name])

            except queue.Empty:
                print(f"Timeout waiting for user response about installing {tool_name}")
                if cls._current_dialog:
                    cls._current_dialog.destroy()
                    cls._current_dialog = None
                return None
