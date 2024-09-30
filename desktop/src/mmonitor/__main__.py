import os
import sys
import traceback

# Determine the absolute path to the 'desktop' directory to append to sys.path
desktop_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if desktop_path not in sys.path:
    sys.path.insert(0, desktop_path)
print(f"Desktop path: {desktop_path}")

# add mmonitor to sys path
mmonitor_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if mmonitor_path not in sys.path:
    sys.path.insert(0, mmonitor_path)
print(f"MMonitor path: {mmonitor_path}")

try:
    from build_mmonitor_pyinstaller import ROOT
    from mmonitor.userside.view import GUI
    print(f"ROOT path: {ROOT}")
    print("Successfully imported GUI and ROOT")
except ImportError as e:
    print(f"Error importing GUI or ROOT: {e}")
    traceback.print_exc()

# set up PYTHONPATH correctly, when running __main__.py

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if base_path not in sys.path:
    sys.path.insert(0, base_path)



def main():
    print("Starting MMonitor application...")
    try:
        app = GUI()
        app.start_app()
        print("Application closed successfully.")
    except Exception as e:
        print(f"Error in main function: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    main()

import sys

def exception_handler(exception_type, exception, traceback):
    print(f"{exception_type.__name__}: {exception}")

sys.excepthook = exception_handler
