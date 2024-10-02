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

# Create and activate a virtual environment
python3 -m venv venv
source venv/bin/activate

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

# Build Samtools
brew install autoconf automake libtool
git clone --recursive https://github.com/samtools/htslib.git
cd htslib
autoheader && autoconf && autoreconf -i && ./configure && make && sudo make install
cd ..
git clone https://github.com/samtools/samtools.git
cd samtools
autoreconf -i && ./configure && make && sudo make install
cd ..

# Build Centrifuger
brew install gcc
git clone https://github.com/mourisl/centrifuger.git
cd centrifuger
make
cd ..

# Install KEGGCharter
conda create -n keggcharter python=3.10 -y
conda activate keggcharter
conda install -y -c conda-forge -c bioconda keggcharter
conda deactivate

# Install Bakta
python3 -m venv bakta-env
source bakta-env/bin/activate
pip install bakta
deactivate

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
python3 -m venv medaka-env
source medaka-env/bin/activate
LDFLAGS="-L/usr/local/opt/openssl@3/lib" CPPFLAGS="-I/usr/local/opt/openssl@3/include" pip install medaka
deactivate

# Install SemiBin2
conda create -n SemiBin python=3.10 -y
conda activate SemiBin
conda install -c bioconda bedtools hmmer samtools -y
conda install -c conda-forge -c bioconda semibin -y
conda install -c pytorch pytorch -y
conda deactivate

# Prepare MMonitor build environment
mkdir -p desktop/src/mmonitor/bin
mkdir -p desktop/src/mmonitor/lib
mkdir -p desktop/src/mmonitor/emu
mkdir -p desktop/src/mmonitor/KEGGCharter

# Copy binaries and libraries
cp $(which flye) desktop/src/mmonitor/bin/
cp $(which samtools) desktop/src/mmonitor/bin/
cp centrifuger/centrifuger desktop/src/mmonitor/bin/
cp $(conda run -n keggcharter which KEGGCharter) desktop/src/mmonitor/bin/
cp bakta-env/bin/bakta desktop/src/mmonitor/bin/
cp -r emu desktop/src/mmonitor/
cp medaka-env/bin/medaka* desktop/src/mmonitor/bin/
cp $(conda run -n SemiBin which SemiBin2) desktop/src/mmonitor/bin/

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