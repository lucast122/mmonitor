#!/usr/bin/env bash
# ==============================================================================
# Medaka wrapper: auto-detect model, fall back to user config, GPU -> CPU retry
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

mkdir -p "${OUTPUT_DIR}"

echo "============================================" | tee "${LOG}"
echo "  Medaka Polishing (auto-detect + fallback)" | tee -a "${LOG}"
echo "============================================" | tee -a "${LOG}"
echo "User-specified model: ${USER_MODEL}" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"

# Helper: run medaka_consensus with optional -m flag and GPU/CPU mode
# Usage: run_medaka <gpu|cpu> [model_flag...]
run_medaka() {
    local device="$1"
    shift
    local extra_args=("$@")
    local device_log="${OUTPUT_DIR}/medaka_${device}.log"

    if [ "${device}" = "cpu" ]; then
        echo "[medaka] Running on CPU (GPU disabled)..." | tee -a "${LOG}"
        CUDA_VISIBLE_DEVICES="" medaka_consensus \
            -i "${READS}" \
            -d "${ASSEMBLY}" \
            -o "${OUTPUT_DIR}" \
            -t "${THREADS}" \
            -f \
            "${extra_args[@]}" \
            2>&1 | tee "${device_log}"
    else
        echo "[medaka] Running on GPU..." | tee -a "${LOG}"
        medaka_consensus \
            -i "${READS}" \
            -d "${ASSEMBLY}" \
            -o "${OUTPUT_DIR}" \
            -t "${THREADS}" \
            -f \
            "${extra_args[@]}" \
            2>&1 | tee "${device_log}"
    fi

    local exit_code=$?
    cat "${device_log}" >> "${LOG}" 2>/dev/null
    return ${exit_code}
}

# Check if result is valid
check_result() {
    [ -f "${OUTPUT_DIR}/consensus.fasta" ] && [ -s "${OUTPUT_DIR}/consensus.fasta" ]
}

# Clean partial output before retry
clean_partial() {
    rm -f "${OUTPUT_DIR}/consensus.fasta" "${OUTPUT_DIR}/calls_to_draft.bam"
}

# ---- Step 1: Try auto-detect model on GPU ------------------------------------
echo "[medaka] Step 1: Auto-detect model, GPU mode..." | tee -a "${LOG}"

run_medaka gpu
if [ $? -eq 0 ] && check_result; then
    echo "[medaka] Success: auto-detect model on GPU." | tee -a "${LOG}"
    exit 0
fi

# ---- Step 2: Auto-detect model on CPU (GPU may have failed) ------------------
echo "" | tee -a "${LOG}"
echo "[medaka] Step 2: Auto-detect model, CPU fallback..." | tee -a "${LOG}"
clean_partial

run_medaka cpu
if [ $? -eq 0 ] && check_result; then
    echo "[medaka] Success: auto-detect model on CPU." | tee -a "${LOG}"
    exit 0
fi

# ---- Step 3: User model on GPU -----------------------------------------------
if [ -z "${USER_MODEL}" ]; then
    echo "[medaka] ERROR: Auto-detect failed and no fallback model specified." | tee -a "${LOG}"
    echo "         Set medaka.model in your config." | tee -a "${LOG}"
    exit 1
fi

echo "" | tee -a "${LOG}"
echo "[medaka] Step 3: User model '${USER_MODEL}', GPU mode..." | tee -a "${LOG}"
clean_partial

run_medaka gpu -m "${USER_MODEL}"
if [ $? -eq 0 ] && check_result; then
    echo "[medaka] Success: user model '${USER_MODEL}' on GPU." | tee -a "${LOG}"
    exit 0
fi

# ---- Step 4: User model on CPU (last resort) ---------------------------------
echo "" | tee -a "${LOG}"
echo "[medaka] Step 4: User model '${USER_MODEL}', CPU fallback..." | tee -a "${LOG}"
clean_partial

run_medaka cpu -m "${USER_MODEL}"
if [ $? -eq 0 ] && check_result; then
    echo "[medaka] Success: user model '${USER_MODEL}' on CPU." | tee -a "${LOG}"
    exit 0
fi

# ---- All attempts failed ------------------------------------------------------
echo "" | tee -a "${LOG}"
echo "[medaka] ERROR: All attempts failed (auto+GPU, auto+CPU, user+GPU, user+CPU)." | tee -a "${LOG}"
exit 1
