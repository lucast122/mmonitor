#!/bin/bash

# Exit on error
set -e

# Create conda environment
echo "Creating CheckM2 environment..."
conda create -n checkm2 python=3.9 -y

# Activate environment and install dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate checkm2

# Install dependencies
echo "Installing dependencies..."
conda install -c conda-forge -c bioconda \
    numpy \
    scipy \
    scikit-learn \
    pytorch \
    diamond \
    hmmer \
    prodigal \
    mmseqs2 \
    -y

# Install CheckM2
echo "Installing CheckM2..."
pip install checkm2

# Verify installation
echo "Verifying installation..."
which checkm2

# Print path for Python to capture
echo "CHECKM2_PATH=$(which checkm2)"
