#!/usr/bin/env bash
set -euo pipefail

# ==============================
# Basic settings
# ==============================
BATCH_SIZE="${BATCH_SIZE:-25}"
THREADS="${THREADS:-1}"
PVAL_THRESHOLD="${PVAL_THRESHOLD:-0.05}"
BIOLEARN_IMPUTATION_METHOD="${BIOLEARN_IMPUTATION_METHOD:-default}"

RUN_EPICV2="${RUN_EPICV2:-1}"
RUN_MSA="${RUN_MSA:-1}"
RUN_FINAL_R_ANALYSIS="${RUN_FINAL_R_ANALYSIS:-1}"

WORK_ROOT="${WORK_ROOT:-/work}"
CONDA_ENV="${CONDA_ENV:-/opt/conda/envs/methylation-beta-pipeline}"
PYTHON_BIN="${CONDA_ENV}/bin/python"
RSCRIPT_BIN="${CONDA_ENV}/bin/Rscript"

SESAME_CACHE="${SESAME_CACHE:-${WORK_ROOT}/vendor/sesame_cache}"
MANIFEST_DIR="${MANIFEST_DIR:-${WORK_ROOT}/manifest}"

# ==============================
# Force serial execution
# ==============================
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1
export THREADS=1
export N_THREADS=1
export NUM_THREADS=1

export EXPERIMENT_HUB_CACHE="${EXPERIMENT_HUB_CACHE:-${SESAME_CACHE}}"
export BFC_CACHE="${BFC_CACHE:-${SESAME_CACHE}}"
export SESAME_CACHE="${SESAME_CACHE}"
export MANIFEST_DIR="${MANIFEST_DIR}"

echo "========================================"
echo "Methylation IDAT pipeline"
echo "BATCH_SIZE=${BATCH_SIZE}"
echo "THREADS=${THREADS}"
echo "PVAL_THRESHOLD=${PVAL_THRESHOLD}"
echo "BIOLEARN_IMPUTATION_METHOD=${BIOLEARN_IMPUTATION_METHOD}"
echo "RUN_EPICV2=${RUN_EPICV2}"
echo "RUN_MSA=${RUN_MSA}"
echo "RUN_FINAL_R_ANALYSIS=${RUN_FINAL_R_ANALYSIS}"
echo "WORK_ROOT=${WORK_ROOT}"
echo "SESAME_CACHE=${SESAME_CACHE}"
echo "MANIFEST_DIR=${MANIFEST_DIR}"
echo "========================================"

# ==============================
# Helper functions
# ==============================
require_file() {
  local f="$1"
  if [ ! -s "${f}" ]; then
    echo "ERROR: required file not found or empty: ${f}"
    exit 1
  fi
}

require_dir() {
  local d="$1"
  if [ ! -d "${d}" ]; then
    echo "ERROR: required directory not found: ${d}"
    exit 1
  fi
}

# ==============================
# Common checks
# ==============================
require_file "${WORK_ROOT}/meta.csv"
require_dir "${SESAME_CACHE}"
require_dir "${MANIFEST_DIR}"

if [ "${RUN_EPICV2}" = "1" ]; then
  require_dir "${WORK_ROOT}/idat_EPICv2"
  require_file "${WORK_ROOT}/name_table/EPICv2_blood.txt"
  require_file "${MANIFEST_DIR}/EPICv2_manifest.csv"
fi

if [ "${RUN_MSA}" = "1" ]; then
  require_dir "${WORK_ROOT}/idat_MSA"
  require_file "${WORK_ROOT}/name_table/MSA_blood.txt"
  require_file "${MANIFEST_DIR}/MSA_manifest.csv"
fi

# These are required for the final EPICv2-vs-MSA comparison.
if [ "${RUN_EPICV2}" = "1" ] && [ "${RUN_MSA}" = "1" ]; then
  require_file "${WORK_ROOT}/templates/training.csv"
  require_file "${WORK_ROOT}/templates/testing.csv"
fi

# ==============================
# EPICv2
# ==============================
if [ "${RUN_EPICV2}" = "1" ]; then
  echo "========================================"
  echo "EPICv2: make IDAT batches"
  echo "========================================"

  "${PYTHON_BIN}" "${WORK_ROOT}/bin/make_idat_batches.py" \
    --src "${WORK_ROOT}/idat_EPICv2" \
    --out "${WORK_ROOT}/idat_batches_EPICv2" \
    --batch-size "${BATCH_SIZE}"

  echo "========================================"
  echo "EPICv2: SeSAMe IDAT to beta"
  echo "========================================"

  THREADS="${THREADS}" \
  PVAL_THRESHOLD="${PVAL_THRESHOLD}" \
  SESAME_CACHE="${SESAME_CACHE}" \
  MANIFEST_DIR="${MANIFEST_DIR}" \
    "${WORK_ROOT}/bin/run_sesame_batches_in_container.sh" EPICv2

  echo "========================================"
  echo "EPICv2: combine collapsed beta matrices"
  echo "========================================"

  "${PYTHON_BIN}" "${WORK_ROOT}/scripts_postprocess/prepare_collapsed_beta.py" \
    --batch-root "${WORK_ROOT}/results/EPICv2_batches" \
    --sample-sheet "${WORK_ROOT}/name_table/EPICv2_blood.txt" \
    --out-dir "${WORK_ROOT}/postprocess/EPICv2" \
    --batch-file-name collapsed_beta_matrix.txt \
    --output-name collapsed_beta_matrix.txt \
    --orientation cpg_rows \
    --drop-unmapped

  echo "========================================"
  echo "EPICv2: Biolearn clocks with default imputation"
  echo "========================================"

  "${PYTHON_BIN}" -u "${WORK_ROOT}/scripts_postprocess/run_clocks_biolearn_default_meta_grimagev1_lowmem.py" \
    --input "${WORK_ROOT}/postprocess/EPICv2/collapsed_beta_matrix.txt" \
    --meta "${WORK_ROOT}/meta.csv" \
    --out-dir "${WORK_ROOT}/postprocess/EPICv2_biolearn_default_grimagev1" \
    --biolearn-imputation-method "${BIOLEARN_IMPUTATION_METHOD}" \
    --allow-missing-age \
    --allow-missing-sex
fi

# ==============================
# MSA
# ==============================
if [ "${RUN_MSA}" = "1" ]; then
  echo "========================================"
  echo "MSA: make IDAT batches"
  echo "========================================"

  "${PYTHON_BIN}" "${WORK_ROOT}/bin/make_idat_batches.py" \
    --src "${WORK_ROOT}/idat_MSA" \
    --out "${WORK_ROOT}/idat_batches_MSA" \
    --batch-size "${BATCH_SIZE}"

  echo "========================================"
  echo "MSA: SeSAMe IDAT to beta"
  echo "========================================"

  THREADS="${THREADS}" \
  PVAL_THRESHOLD="${PVAL_THRESHOLD}" \
  SESAME_CACHE="${SESAME_CACHE}" \
  MANIFEST_DIR="${MANIFEST_DIR}" \
    "${WORK_ROOT}/bin/run_sesame_batches_in_container.sh" MSA

  echo "========================================"
  echo "MSA: combine collapsed beta matrices"
  echo "========================================"

  "${PYTHON_BIN}" "${WORK_ROOT}/scripts_postprocess/prepare_collapsed_beta.py" \
    --batch-root "${WORK_ROOT}/results/MSA_batches" \
    --sample-sheet "${WORK_ROOT}/name_table/MSA_blood.txt" \
    --out-dir "${WORK_ROOT}/postprocess/MSA" \
    --batch-file-name collapsed_beta_matrix.txt \
    --output-name collapsed_beta_matrix.txt \
    --orientation cpg_rows \
    --drop-unmapped

  echo "========================================"
  echo "MSA: Biolearn clocks with default imputation"
  echo "========================================"

  "${PYTHON_BIN}" -u "${WORK_ROOT}/scripts_postprocess/run_clocks_biolearn_default_meta_grimagev1_lowmem.py" \
    --input "${WORK_ROOT}/postprocess/MSA/collapsed_beta_matrix.txt" \
    --meta "${WORK_ROOT}/meta.csv" \
    --out-dir "${WORK_ROOT}/postprocess/MSA_biolearn_default_grimagev1" \
    --biolearn-imputation-method "${BIOLEARN_IMPUTATION_METHOD}" \
    --allow-missing-age \
    --allow-missing-sex
fi

# ==============================
# EPICv2 vs MSA comparison
# ==============================
if [ "${RUN_EPICV2}" = "1" ] && [ "${RUN_MSA}" = "1" ]; then
  echo "========================================"
  echo "Create 12 R-input files from EPICv2/MSA clock.csv"
  echo "========================================"

  "${PYTHON_BIN}" "${WORK_ROOT}/scripts_postprocess/create_clock_training_testing_files.py" \
    --epicv2-clock "${WORK_ROOT}/postprocess/EPICv2_biolearn_default_grimagev1/clock.csv" \
    --msa-clock "${WORK_ROOT}/postprocess/MSA_biolearn_default_grimagev1/clock.csv" \
    --training-template "${WORK_ROOT}/templates/training.csv" \
    --testing-template "${WORK_ROOT}/templates/testing.csv" \
    --out-dir "${WORK_ROOT}/clock_compare_outputs_grimagev1" \
    --clocks Horvath Hannum PhenoAge GrimAgeV2 DunedinPACE GrimAgeV1

  if [ "${RUN_FINAL_R_ANALYSIS}" = "1" ]; then
    echo "========================================"
    echo "Run final R analysis script"
    echo "========================================"

    "${WORK_ROOT}/bin/run_final_r_analysis.sh"
  fi
else
  echo "========================================"
  echo "Skip EPICv2-vs-MSA comparison"
  echo "Reason: RUN_EPICV2=${RUN_EPICV2}, RUN_MSA=${RUN_MSA}"
  echo "========================================"
fi

echo "========================================"
echo "Pipeline finished."
echo "Outputs:"
echo "  ${WORK_ROOT}/results"
echo "  ${WORK_ROOT}/postprocess"
echo "  ${WORK_ROOT}/clock_compare_outputs_grimagev1"
echo "  ${WORK_ROOT}/results_grimagev1"
echo "========================================"
