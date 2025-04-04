#!/bin/bash

# Set default paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../../" &> /dev/null && pwd)"

# Add the script directory to PYTHONPATH
export PYTHONPATH=$SCRIPT_DIR:$PYTHONPATH

# Default parameters
INDEX_PREFIX="$ROOT_DIR/src/resources/custom_centrifuger_db/ncbi_build_20250401_110630/cfr_ncbi"
WATCH_DIR="$ROOT_DIR/desktop/src/resources/pipeline_out/test_reads/watch_dir"
OUTPUT_FILE="$ROOT_DIR/desktop/src/resources/pipeline_out/test_reads/test_reads_centrifuger_classifications.txt"
THREADS=12
MIN_FILES=5
MODE="auto"  # Options: auto, realtime, batch

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --index-prefix)
      INDEX_PREFIX="$2"
      shift 2
      ;;
    --watch-dir)
      WATCH_DIR="$2"
      shift 2
      ;;
    --output)
      OUTPUT_FILE="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --min-files)
      MIN_FILES="$2"
      shift 2
      ;;
    --mode)
      MODE="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Verify that index prefix directory exists
if [ ! -d "$(dirname "$INDEX_PREFIX")" ]; then
  echo "Error: Index prefix directory doesn't exist: $(dirname "$INDEX_PREFIX")"
  echo "Please provide a valid index prefix with --index-prefix"
  exit 1
fi

# Print parameters
echo "Running with parameters:"
echo "  Index prefix: $INDEX_PREFIX"
echo "  Watch directory: $WATCH_DIR"
echo "  Output file: $OUTPUT_FILE"
echo "  Threads: $THREADS"
echo "  Min files: $MIN_FILES"
echo "  Mode: $MODE"

# Create directories if they don't exist
mkdir -p "$WATCH_DIR"
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Confirm Python version
echo "Using Python: $(which python3)"
python3 --version

# Check for centrifuger executables
echo "Checking for centrifuger executables..."
if which centrifuger > /dev/null; then
  echo "Found centrifuger: $(which centrifuger)"
else
  echo "Warning: centrifuger not found in PATH"
fi

if which centrifuger-realtime > /dev/null; then
  echo "Found centrifuger-realtime: $(which centrifuger-realtime)"
else
  echo "Warning: centrifuger-realtime not found in PATH"
fi

# Run the test script
cd "$SCRIPT_DIR"
echo "Running test script from: $SCRIPT_DIR"
python3 test_centrifuger.py \
  --index-prefix "$INDEX_PREFIX" \
  --watch-dir "$WATCH_DIR" \
  --output "$OUTPUT_FILE" \
  --threads "$THREADS" \
  --min-files "$MIN_FILES" \
  --mode "$MODE"

# Check exit status
STATUS=$?
if [ $STATUS -eq 0 ]; then
  echo "Test completed successfully"
elif [ $STATUS -eq 130 ]; then
  echo "Test terminated by user"
else
  echo "Test failed with exit code $STATUS"
  echo "Check test_centrifuger.log for details"
fi

exit $STATUS 