import os
import sys

# Determine the absolute path to the 'desktop' directory to append to sys.path
# This allows importing the ROOT variable from build_mmonitor_pyinstaller
desktop_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if desktop_path not in sys.path:
    sys.path.insert(0, desktop_path)

# add mmonitor to sys path
mmonitor_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if mmonitor_path not in sys.path:
    sys.path.insert(0, mmonitor_path)
    
from build_mmonitor_pyinstaller import ROOT
from mmonitor.userside.view import GUI

# set up PYTHONPATH correctly, when running __main__.py

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if base_path not in sys.path:
    sys.path.insert(0, base_path)


print(ROOT)

def main():
    GUI().start_app()
    # EnhancedView().start_app()


if __name__ == '__main__':
    main()
        