import os
import json

# Setup paths
project_root = os.path.dirname(os.path.abspath(__file__))
main_path = os.path.join(project_root, 'desktop', 'src', 'mmonitor', '__main__.py')
hooks_path = os.path.join(project_root, 'hooks')
runtime_hooks_path = os.path.join(project_root, 'runtime_hooks.py')
src_dir = os.path.join(project_root, 'desktop', 'src')
resources_dir = os.path.join(project_root, 'desktop', 'src', 'resources')
emu_db_path = os.path.join(resources_dir, 'emu_db')

# Create themes directory if it doesn't exist
themes_dir = os.path.join(resources_dir, 'themes')
os.makedirs(themes_dir, exist_ok=True)

# Theme path - create if it doesn't exist
theme_path = os.path.join(themes_dir, 'grey_theme.json')
if not os.path.exists(theme_path):
    print(f"Creating theme file at {theme_path}")
    with open(theme_path, 'w') as f:
        json.dump({
            "name": "grey",
            "type": "dark",
            "colors": {
                "primary": "#2b2b2b",
                "secondary": "#3c3f41",
                "accent": "#4b6eaf"
            }
        }, f)
    print(f"Created default theme file at {theme_path}")

# Define required files with correct paths
required_files = [
    (emu_db_path, 'resources/emu_db'),
    (theme_path, 'resources/themes'),
    (os.path.join(resources_dir, 'db_config.json'), 'resources'),
]

# Add data files with absolute paths and correct destinations
for src, dst in required_files:
    if os.path.exists(src):
        cmd.extend([f'--add-data={src}:{dst}'])
        print(f"Adding data file: {src} -> {dst}")
    else:
        print(f"Warning: Required file not found: {src}") 