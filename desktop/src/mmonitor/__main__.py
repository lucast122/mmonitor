import os
import sys
import traceback
# fix module not found mmonitor
# Add the parent directory of 'src' to sys.path
src_parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, src_parent_dir)
src_parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, src_parent_dir)

from mmonitor.userside.view import GUI

# Determine the absolute path to the 'src' directory to append to sys.path
# Now we can import from mmonitor


import tkinter as tk
from mmonitor.userside.view import GUI

def main():
    app = GUI()  # No need to pass root
    app.mainloop()

if __name__ == "__main__":
    main()
def exception_handler(exception_type, exception, traceback):
    print(f"{exception_type.__name__}: {exception}")

sys.excepthook = exception_handler


