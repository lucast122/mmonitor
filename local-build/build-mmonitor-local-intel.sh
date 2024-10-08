#!/bin/bash
set -e

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check for required tools
for cmd in python3 git brew conda; do
    if ! command_exists $cmd; then
        echo "$cmd is required but not installed. Please install it and try again."
        exit 1
    fi
done

# Install Python dependencies
python3 -m venv pyinstaller-env
source pyinstaller-env/bin/activate
python3 -m pip install --upgrade pip
sed '/^tkinter/d' desktop/requirements.txt > desktop/requirements_modified.txt
pip install -r desktop/requirements_modified.txt
pip install pyinstaller pysam numpy

# Prepare MMonitor build environment
mkdir -p desktop/src/mmonitor/bin
mkdir -p desktop/src/mmonitor/emu
mkdir -p desktop/src/mmonitor/KEGGCharter
mkdir -p desktop/src/mmonitor/lib/python3.10/site-packages

# Clone and build Flye 2.9.5-b1801
git clone https://github.com/fenderglass/Flye.git
cd Flye
git checkout 2.9.5-b1801
make
cd ..

# Copy the built Flye to the appropriate directory
cp Flye/bin/flye desktop/src/mmonitor/bin/

# Find and copy libhts.so (for Linux) or libhts.dylib (for macOS)
LIBHTS=$(find /usr/local/lib -name "libhts.*" | head -n 1)
if [ -n "$LIBHTS" ]; then
  cp "$LIBHTS" desktop/src/mmonitor/lib/
else
  echo "libhts not found. Samtools might not work correctly."
fi

# Adjust shebang lines
sed -i '' '1s|^#!.*|#!/usr/bin/env python3|' desktop/src/mmonitor/bin/bakta
sed -i '' '1s|^#!.*|#!/usr/bin/env python3|' desktop/src/mmonitor/bin/medaka
sed -i '' '1s|^#!.*|#!/usr/bin/env python3|' desktop/src/mmonitor/bin/SemiBin2

# Set PYTHONPATH and DYLD_LIBRARY_PATH (use DYLD_LIBRARY_PATH for macOS)
export PYTHONPATH="$PWD/desktop:$PWD/desktop/src:$PWD/desktop/src/mmonitor:$PWD/desktop/src/mmonitor/lib/python3.10/site-packages"
export DYLD_LIBRARY_PATH="$PWD/desktop/src/mmonitor/lib:$DYLD_LIBRARY_PATH"

# Install Tcl/Tk
brew install tcl-tk

# Build MMonitor
cd desktop/src/mmonitor
TCL_LIBRARY=$(brew --prefix tcl-tk)/lib/tcl8.6
TK_LIBRARY=$(brew --prefix tcl-tk)/lib/tk8.6
pyinstaller --onefile --add-data "bin:bin" --add-data "KEGGCharter:KEGGCharter" --add-data "emu:emu" \
  --add-data "lib:lib" \
  --add-binary "/usr/local/opt/tcl-tk/lib/libtcl8.6.dylib:." \
  --add-binary "/usr/local/opt/tcl-tk/lib/libtk8.6.dylib:." \
  --add-data "$TCL_LIBRARY:tcl" \
  --add-data "$TK_LIBRARY:tk" \
  --hidden-import=pysam \
  --hidden-import=build_mmonitor_pyinstaller \
  __main__.py --distpath ../../output

echo "Build complete. The MMonitor binary is in desktop/output/__main__"

# Test individual binaries
echo "Testing Flye:"
FLYE_VERSION=$($PWD/desktop/src/mmonitor/bin/flye --version)
if [[ "$FLYE_VERSION" == *"2.9.5-b1801"* ]]; then
    echo "Flye version is correct: $FLYE_VERSION"
else
    echo "Error: Flye version is incorrect. Expected 2.9.5-b1801, got $FLYE_VERSION"
    exit 1
fi
echo "Testing Samtools:"
$PWD/desktop/src/mmonitor/bin/samtools --version || true
echo "Testing Centrifuger:"
$PWD/desktop/src/mmonitor/bin/centrifuger -v || true
echo "Testing KEGGCharter:"
python $PWD/desktop/src/mmonitor/KEGGCharter/KEGGCharter.py --help || true
echo "Testing Bakta:"
$PWD/desktop/src/mmonitor/bin/bakta --version || true
echo "Testing EMU:"
$PWD/desktop/src/mmonitor/emu/emu --version || true
echo "Testing Medaka:"
$PWD/desktop/src/mmonitor/bin/medaka --version || true
echo "Testing SemiBin2:"
$PWD/desktop/src/mmonitor/bin/SemiBin2 --version || true

# Test MMonitor binary
cd ../../output
chmod +x __main__
echo "Testing MMonitor version:"
./__main__ --version || true
echo "Testing MMonitor help:"
./__main__ --help || true

# Deactivate the virtual environment
deactivate