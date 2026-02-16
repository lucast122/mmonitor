#!/bin/bash
# Test script to run single genome annotation pipeline
# Downloads a small reference genome and runs Bakta annotation

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Configuration
GENOME_ID="Staph_aureus_NCTC8325"
SAMPLE_NAME="test_genome_run"
RESULTS_DIR="results/${SAMPLE_NAME}"
NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz"

echo "==================================================="
echo "MMonitor Single Genome Annotation Test"
echo "==================================================="
echo "Genome: ${GENOME_ID}"
echo "Sample: ${SAMPLE_NAME}"
echo ""

# Create directories
mkdir -p "${RESULTS_DIR}/genome"
mkdir -p "${RESULTS_DIR}/logs"

# Step 1: Download genome from NCBI
echo "[1/4] Downloading genome from NCBI..."
if [ ! -f "${RESULTS_DIR}/genome/${GENOME_ID}.fasta" ]; then
    echo "  Downloading from: ${NCBI_URL}"
    wget -q -O "${RESULTS_DIR}/genome/${GENOME_ID}.fasta.gz" "${NCBI_URL}" || {
        echo "ERROR: Failed to download genome!"
        exit 1
    }
    gunzip -f "${RESULTS_DIR}/genome/${GENOME_ID}.fasta.gz"
    echo "  Downloaded: ${RESULTS_DIR}/genome/${GENOME_ID}.fasta"
else
    echo "  Already downloaded"
fi

# Verify download
if [ ! -f "${RESULTS_DIR}/genome/${GENOME_ID}.fasta" ]; then
    echo "ERROR: Genome file not found!"
    exit 1
fi

echo "  File size: $(ls -lh ${RESULTS_DIR}/genome/${GENOME_ID}.fasta | awk '{print $5}')"
echo "  Sequences: $(grep -c '^>' ${RESULTS_DIR}/genome/${GENOME_ID}.fasta)"

# Activate conda for snakemake
echo ""
echo "[2/4] Setting up environment..."
source ~/miniforge3/etc/profile.d/conda.sh
conda activate base

# Step 3: Run Snakemake pipeline
echo ""
echo "[3/4] Running annotation pipeline via Snakemake..."
echo "  This will create conda environments on first run (may take a while)..."
echo ""

# Run the complete pipeline
snakemake \
    --snakefile Snakefile \
    --configfile config/test_genome.yaml \
    --config input_fasta="${SCRIPT_DIR}/${RESULTS_DIR}/genome/${GENOME_ID}.fasta" \
    --cores 8 \
    --use-conda \
    --conda-frontend mamba \
    --rerun-incomplete \
    annotate_genome_complete \
    2>&1 | tee "${RESULTS_DIR}/logs/snakemake.log"

echo ""
echo "[4/4] Checking results..."

echo ""
echo "==================================================="
echo "Pipeline complete!"
echo "==================================================="
echo ""
echo "Results are in: ${RESULTS_DIR}"
echo ""
echo "Files created:"
echo "  Annotation:"
ls -la "${RESULTS_DIR}/annotation/${GENOME_ID}/" 2>/dev/null || echo "    (not found)"
echo ""
echo "  Upload:"
ls -la "${RESULTS_DIR}/upload/" 2>/dev/null || echo "    (not found)"
echo ""

# Check upload status
if [ -f "${RESULTS_DIR}/upload/${GENOME_ID}_complete.json" ]; then
    echo "Upload status:"
    cat "${RESULTS_DIR}/upload/${GENOME_ID}_complete.json"
    echo ""
fi

echo ""
echo "To view in the web interface:"
echo "  1. Start the Django server: cd ../../../server && python manage.py runserver"
echo "  2. Open http://localhost:8000 and log in as 'ulli'"
echo "  3. Go to the MAGs tab to see the genome browser"
