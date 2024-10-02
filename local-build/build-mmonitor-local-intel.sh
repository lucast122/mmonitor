# Clean up previous build artifacts
echo "Cleaning up previous build artifacts..."
rm -rf desktop/output
rm -rf desktop/src/mmonitor/bin
rm -rf desktop/src/mmonitor/emu
rm -rf desktop/src/mmonitor/KEGGCharter
rm -rf desktop/src/mmonitor/lib

# Remove previously cloned repositories
echo "Removing previously cloned repositories..."
rm -rf desktop/lib/Flye
rm -rf htslib
rm -rf samtools
rm -rf KEGGCharter
rm -rf emu
rm -rf centrifuger

# Remove Python virtual environments
echo "Removing Python virtual environments..."
rm -rf bakta-env
rm -rf medaka-env
rm -rf pyinstaller-env

# Clean conda environment (if using conda)
if command -v conda >/dev/null 2>&1; then
    echo "Removing conda environment..."
    conda remove --name mmonitor_build --all -y
fi

# Clean pip cache
echo "Cleaning pip cache..."
pip cache purge

echo "Cleanup completed. Ready for a fresh build."

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

# Create and activate a conda environment
# conda create -n mmonitor_build python=3.10 -y

# Install Python dependencies
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
pip install pyinstaller

# Build Flye
if [ ! -d "desktop/lib/Flye" ]; then
    git clone https://github.com/fenderglass/Flye.git desktop/lib/Flye
fi
cd desktop/lib/Flye
pip install -e .
cd ../../..

# Check if htslib directory exists and is not empty
if [ -d "htslib" ] && [ "$(ls -A htslib)" ]; then
    echo "htslib directory already exists and is not empty."
    echo "Removing existing htslib directory..."
    rm -rf htslib
fi

# Now clone htslib
git clone --recursive https://github.com/samtools/htslib.git

# Build Samtools
brew install autoconf automake libtool
cd htslib
autoheader && autoconf && autoreconf -i && ./configure
make
sudo make install
sudo chmod 644 /usr/local/lib/pkgconfig/htslib.pc
sudo install -p annot-tsv bgzip htsfile tabix /usr/local/bin
sudo install -p -m 644 htslib/*.h /usr/local/include/htslib
cd ..

# Check if samtools directory exists and is not empty
if [ -d "samtools" ] && [ "$(ls -A samtools)" ]; then
    echo "samtools directory already exists and is not empty."
    echo "Removing existing samtools directory..."
    rm -rf samtools
fi

git clone https://github.com/samtools/samtools.git
cd samtools
autoreconf -i && ./configure
make
sudo make install
cd ..

# Build Centrifuger
brew install gcc
git clone https://github.com/mourisl/centrifuger.git
cd centrifuger
make
cd ..

# Install KEGGCharter
conda install -y -c conda-forge -c bioconda keggcharter

# Install Bakta
pip install bakta

# Clone EMU
mkdir emu && cd emu
git init
git remote add origin https://github.com/treangenlab/emu.git
git config core.sparseCheckout true
echo '/*' > .git/info/sparse-checkout
echo '!emu_database' >> .git/info/sparse-checkout
echo '!example' >> .git/info/sparse-checkout
echo '!example_customdb' >> .git/info/sparse-checkout
git fetch --depth 1 origin --tags
latest_tag=$(git describe --tags `git rev-list --tags --max-count=1`)
git checkout $latest_tag
cd ..

# Build Medaka
brew install openssl
LDFLAGS="-L/usr/local/opt/openssl@3/lib" CPPFLAGS="-I/usr/local/opt/openssl@3/include" pip install medaka

# Install SemiBin2
conda install -c bioconda bedtools hmmer samtools -y
conda install -c conda-forge -c bioconda semibin -y
conda install -c pytorch pytorch -y

# Prepare MMonitor build environment
mkdir -p desktop/src/mmonitor/bin
mkdir -p desktop/src/mmonitor/lib
mkdir -p desktop/src/mmonitor/emu
mkdir -p desktop/src/mmonitor/KEGGCharter

# Copy binaries and libraries
cp $(which flye) desktop/src/mmonitor/bin/
cp $(which samtools) desktop/src/mmonitor/bin/
cp centrifuger/centrifuger desktop/src/mmonitor/bin/
cp $(which KEGGCharter) desktop/src/mmonitor/bin/
cp $(which bakta) desktop/src/mmonitor/bin/
cp -r emu desktop/src/mmonitor/
cp $(which medaka) desktop/src/mmonitor/bin/
cp $(which SemiBin2) desktop/src/mmonitor/bin/

# Install Tcl/Tk
brew install tcl-tk

# Build MMonitor
cd desktop/src/mmonitor
TCL_LIBRARY=$(brew --prefix tcl-tk)/lib/tcl8.6
TK_LIBRARY=$(brew --prefix tcl-tk)/lib/tk8.6
pyinstaller --onefile --add-data "bin:bin" --add-data "KEGGCharter:KEGGCharter" --add-data "emu:emu" \
  --add-binary "/usr/local/opt/tcl-tk/lib/libtcl8.6.dylib:." \
  --add-binary "/usr/local/opt/tcl-tk/lib/libtk8.6.dylib:." \
  --add-data "$TCL_LIBRARY:tcl" \
  --add-data "$TK_LIBRARY:tk" \
  __main__.py --distpath ../../output

echo "Build complete. The MMonitor binary is in desktop/output/__main__"

# Deactivate and remove the conda environment
conda deactivate
conda env remove -n mmonitor_build -y