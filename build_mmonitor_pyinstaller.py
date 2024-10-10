import os
import sys
import shutil
import PyInstaller.__main__

# Define the root directory of your project
ROOT = os.path.abspath(os.path.dirname(__file__))
IMAGES_PATH = os.path.join(ROOT, "src", "resources")

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(ROOT, "src"))

# Define the entry point of your application
ENTRY_POINT = os.path.join(ROOT, "src", "mmonitor", "__main__.py")

# Define additional data files to include
additional_data = [
    (os.path.join(ROOT, "src", "resources", "*"), "resources"),
    (os.path.join(ROOT, "src", "mmonitor", "dashapp", "assets"), "mmonitor/dashapp/assets"),
    (os.path.join(ROOT, "src", "mmonitor", "dashapp", "components"), "mmonitor/dashapp/components"),
    (os.path.join(ROOT, "src", "mmonitor", "dashapp", "pages"), "mmonitor/dashapp/pages"),
]

# Define PyInstaller options
pyinstaller_options = [
    ENTRY_POINT,
    "--name=MMonitor",
    "--windowed",
    "--onefile",
    f"--add-data={os.path.join(ROOT, 'src', 'resources')}:resources",
    "--icon=" + os.path.join(ROOT, "src", "resources", "logo.tiff"),
    "--hidden-import=tkinter",
    "--hidden-import=customtkinter",
    "--hidden-import=PIL",
    "--hidden-import=pandas",
    "--hidden-import=dash",
    "--hidden-import=plotly",
]

# Add additional data files to PyInstaller options
for src, dst in additional_data:
    pyinstaller_options.append(f"--add-data={src}:{dst}")

# Run PyInstaller
PyInstaller.__main__.run(pyinstaller_options)

# Copy additional files to the dist directory
dist_dir = os.path.join(ROOT, "dist")
for src, dst in additional_data:
    dst_path = os.path.join(dist_dir, "MMonitor", dst)
    if not os.path.exists(dst_path):
        os.makedirs(dst_path)
    for file in os.listdir(src):
        shutil.copy2(os.path.join(src, file), dst_path)

print("Build completed successfully!")