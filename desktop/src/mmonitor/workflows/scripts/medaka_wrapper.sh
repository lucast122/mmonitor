#!/usr/bin/env bash
# ==============================================================================
# Medaka wrapper: auto-detect model, fall back to user config, warn on mismatch
# ==============================================================================
#
# Usage (called by Snakemake):
#   medaka_wrapper.sh <reads> <assembly> <output_dir> <threads> <user_model> <logfile>
#
set -uo pipefail

READS="$1"
ASSEMBLY="$2"
OUTPUT_DIR="$3"
THREADS="$4"
USER_MODEL="$5"
LOG="$6"

echo "============================================" | tee "${LOG}"
echo "  Medaka Polishing (auto-detect + fallback)" | tee -a "${LOG}"
echo "============================================" | tee -a "${LOG}"
echo "User-specified model: ${USER_MODEL}" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"

# ---- Step 1: Try auto-detect (no -m flag) ------------------------------------
echo "[medaka] Attempting automatic model selection..." | tee -a "${LOG}"

AUTO_LOG="${OUTPUT_DIR}/medaka_auto.log"
mkdir -p "${OUTPUT_DIR}"

medaka_consensus \
    -i "${READS}" \
    -d "${ASSEMBLY}" \
    -o "${OUTPUT_DIR}" \
    -t "${THREADS}" \
    -f \
    2>&1 | tee "${AUTO_LOG}"

AUTO_EXIT=$?

# ---- Step 2: Extract auto-detected model from log ----------------------------
DETECTED_MODEL=""
if [ -f "${AUTO_LOG}" ]; then
    # medaka prints "Selecting model: <model>" or "Using model: <model>" in logs
    DETECTED_MODEL=$(grep -oP '(?:Selecting|Using|Autoselected) model[:\s]+\K\S+' "${AUTO_LOG}" | head -1 || true)
    if [ -z "${DETECTED_MODEL}" ]; then
        # Alternative pattern: "Model <model>" or just the model name after "model"
        DETECTED_MODEL=$(grep -oP 'model\s*[:=]\s*\K\S+' "${AUTO_LOG}" | head -1 || true)
    fi
fi

# ---- Step 3: Handle auto-detect result ---------------------------------------
if [ "${AUTO_EXIT}" -eq 0 ] && [ -f "${OUTPUT_DIR}/consensus.fasta" ]; then
    echo "" | tee -a "${LOG}"
    echo "[medaka] Auto-detection succeeded." | tee -a "${LOG}"

    if [ -n "${DETECTED_MODEL}" ]; then
        echo "[medaka] Auto-detected model: ${DETECTED_MODEL}" | tee -a "${LOG}"

        # Compare with user-specified model
        if [ -n "${USER_MODEL}" ] && [ "${DETECTED_MODEL}" != "${USER_MODEL}" ]; then
            echo "" | tee -a "${LOG}"
            echo "WARNING: Auto-detected model '${DETECTED_MODEL}' differs from" | tee -a "${LOG}"
            echo "         user-specified model '${USER_MODEL}'." | tee -a "${LOG}"
            echo "         Using auto-detected model. If results look wrong," | tee -a "${LOG}"
            echo "         re-run with the correct model in your config." | tee -a "${LOG}"
            echo "" | tee -a "${LOG}"
        fi
    else
        echo "[medaka] Could not determine which model was auto-selected." | tee -a "${LOG}"
    fi

    # Append auto-detect log into main log
    cat "${AUTO_LOG}" >> "${LOG}" 2>/dev/null
    exit 0
fi

# ---- Step 4: Auto-detect failed — fall back to user model --------------------
echo "" | tee -a "${LOG}"
echo "[medaka] Auto-detection failed (exit code: ${AUTO_EXIT})." | tee -a "${LOG}"

if [ -z "${USER_MODEL}" ]; then
    echo "[medaka] ERROR: No fallback model specified. Set medaka.model in config." | tee -a "${LOG}"
    cat "${AUTO_LOG}" >> "${LOG}" 2>/dev/null
    exit 1
fi

echo "[medaka] Falling back to user-specified model: ${USER_MODEL}" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"

# Clean partial output before retry
rm -f "${OUTPUT_DIR}/consensus.fasta" "${OUTPUT_DIR}/calls_to_draft.bam"

FALLBACK_LOG="${OUTPUT_DIR}/medaka_fallback.log"

medaka_consensus \
    -i "${READS}" \
    -d "${ASSEMBLY}" \
    -o "${OUTPUT_DIR}" \
    -t "${THREADS}" \
    -m "${USER_MODEL}" \
    -f \
    2>&1 | tee "${FALLBACK_LOG}"

FALLBACK_EXIT=$?

cat "${FALLBACK_LOG}" >> "${LOG}" 2>/dev/null

if [ "${FALLBACK_EXIT}" -eq 0 ] && [ -f "${OUTPUT_DIR}/consensus.fasta" ]; then
    echo "" | tee -a "${LOG}"
    echo "[medaka] Fallback with model '${USER_MODEL}' succeeded." | tee -a "${LOG}"
    exit 0
else
    echo "" | tee -a "${LOG}"
    echo "[medaka] ERROR: Both auto-detect and fallback model '${USER_MODEL}' failed." | tee -a "${LOG}"
    exit 1
fi
