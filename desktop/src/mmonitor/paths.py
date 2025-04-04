import os
import sys

# Get the application's path dynamically
def get_app_path():
    if getattr(sys, 'frozen', False):
        # Running as bundled executable
        app_path = os.path.dirname(sys.executable)
        
        if sys.platform == 'darwin':
            # On macOS, we need to handle both .app and _internal structures
            if 'Contents/MacOS' in app_path:
                # We're in the .app bundle
                # app_path is /Path/to/MMonitor.app/Contents/MacOS
                app_path = os.path.abspath(os.path.join(app_path, '..', '..', '..'))
                
                # Check if we have an _internal directory alongside the .app
                internal_dir = os.path.join(os.path.dirname(app_path), "MMonitor", "_internal")
                if os.path.exists(internal_dir):
                    print(f"Found _internal directory: {internal_dir}")
                    return os.path.dirname(internal_dir)  # Return the MMonitor directory containing _internal
                
                # If we don't have _internal alongside .app, use the .app directory
                print(f"Using .app directory as app path: {app_path}")
                return app_path
            else:
                # We might be in the _internal directory structure
                if '_internal' in app_path:
                    # Find the root directory (containing _internal)
                    root_dir = app_path
                    while os.path.basename(root_dir) != '_internal' and os.path.dirname(root_dir) != root_dir:
                        root_dir = os.path.dirname(root_dir)
                    
                    if os.path.basename(root_dir) == '_internal':
                        app_path = os.path.dirname(root_dir)
                        print(f"Using _internal parent as app path: {app_path}")
                        return app_path
        
        # For other platforms or fallback
        return app_path
    else:
        # Running in development mode
        # Get the directory containing mmonitor module
        module_dir = os.path.dirname(os.path.abspath(__file__))
        # Go up to src directory
        src_dir = os.path.abspath(os.path.join(module_dir, '..', '..'))
        # Project root is one level above src
        return os.path.abspath(os.path.join(src_dir, '..'))

# Get the project root directory path
PROJECT_ROOT = get_app_path()
print(f"PROJECT_ROOT set to: {PROJECT_ROOT}")

# Define important paths
SRC_DIR = os.path.join(PROJECT_ROOT, "desktop", "src")
RESOURCES_DIR = os.path.join(SRC_DIR, "resources")
IMAGES_DIR = os.path.join(RESOURCES_DIR, "images")
LIB_DIR = os.path.join(SRC_DIR, "mmonitor", "lib")

# Check if we're in a bundled app with _internal directory
if getattr(sys, 'frozen', False) and os.path.exists(os.path.join(PROJECT_ROOT, "_internal")):
    INTERNAL_DIR = os.path.join(PROJECT_ROOT, "_internal")
    print(f"Using INTERNAL_DIR: {INTERNAL_DIR}")
    
    # Update paths to use _internal structure if we're in the bundled app
    RESOURCES_DIR = os.path.join(INTERNAL_DIR, "mmonitor", "src", "resources")
    IMAGES_DIR = os.path.join(RESOURCES_DIR, "images")
    LIB_DIR = os.path.join(INTERNAL_DIR, "lib")
    
    # Define CENTRIFUGER_DIR to help locate the binary
    CENTRIFUGER_DIR = os.path.join(INTERNAL_DIR, "lib", "centrifuger-realtime", "centrifuger")
    print(f"Set CENTRIFUGER_DIR: {CENTRIFUGER_DIR}")
else:
    # In development mode or .app without _internal
    CENTRIFUGER_DIR = os.path.join(PROJECT_ROOT, "lib", "centrifuger-realtime", "centrifuger")

# Set ROOT the same as PROJECT_ROOT for compatibility
ROOT = PROJECT_ROOT
SERVER_DIR = os.path.join(ROOT, "server")

# Print detected paths to help with debugging
print(f"SRC_DIR: {SRC_DIR}")
print(f"RESOURCES_DIR: {RESOURCES_DIR}")
print(f"IMAGES_DIR: {IMAGES_DIR}")
print(f"LIB_DIR: {LIB_DIR}")
print(f"ROOT: {ROOT}")
print(f"SERVER_DIR: {SERVER_DIR}")
