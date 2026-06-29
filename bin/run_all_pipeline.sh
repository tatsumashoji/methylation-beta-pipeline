#!/usr/bin/env bash
set -euo pipefail

# ==============================
# User-configurable settings
# ==============================
BATCH_SIZE="${BATCH_SIZE:-25}"
THREADS=1
PVAL_THRESHOLD="${PVAL_THRESHOLD:-0.05}"
BIOLEARN_IMPUTATION_METHOD="${BIOLEARN_IMPUTATION_METHOD:-default}"

RUN_EPICV2="${RUN_EPICV2:-1}"
RUN_MSA="${RUN_MSA:-1}"
RUN_QC_REPORT="${RUN_QC_REPORT:-1}"
RUN_BIOLEARN_CLOCKS="${RUN_BIOLEARN_CLOCKS:-1}"
RUN_FINAL_R_ANALYSIS="${RUN_FINAL_R_ANALYSIS:-1}"

# Large optional outputs are disabled by default to reduce disk usage.
# Set any of these to 1 only when the corresponding intermediate matrix is needed.
SAVE_RAW_BETA="${SAVE_RAW_BETA:-0}"
SAVE_NONCOLLAPSED_BETA="${SAVE_NONCOLLAPSED_BETA:-0}"
SAVE_DETECTION_PVALS="${SAVE_DETECTION_PVALS:-0}"
SAVE_INTENSITY="${SAVE_INTENSITY:-0}"

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
export SESAME_CACHE
export MANIFEST_DIR
export RSCRIPT_BIN
export SAVE_RAW_BETA
export SAVE_NONCOLLAPSED_BETA
export SAVE_DETECTION_PVALS
export SAVE_INTENSITY

# ==============================
# Helpers
# ==============================
require_file() {
  local path="$1"
  if [ ! -s "${path}" ]; then
    echo "ERROR: required file not found or empty: ${path}" >&2
    exit 1
  fi
}

require_dir() {
  local path="$1"
  if [ ! -d "${path}" ]; then
    echo "ERROR: required directory not found: ${path}" >&2
    exit 1
  fi
}

require_executable() {
  local path="$1"
  if [ ! -x "${path}" ]; then
    echo "ERROR: required executable not found: ${path}" >&2
    exit 1
  fi
}

is_enabled() {
  [ "$1" = "1" ]
}

combine_collapsed_beta() {
  local dataset="$1"
  local batch_root="$2"
  local sample_sheet="$3"
  local out_dir="$4"

  echo "========================================"
  echo "${dataset}: combine collapsed beta matrices"
  echo "========================================"

  "${PYTHON_BIN}" "${WORK_ROOT}/scripts_postprocess/prepare_collapsed_beta.py" \
    --batch-root "${batch_root}" \
    --sample-sheet "${sample_sheet}" \
    --out-dir "${out_dir}" \
    --batch-file-name collapsed_beta_matrix.txt \
    --output-name collapsed_beta_matrix.txt \
    --orientation cpg_rows \
    --drop-unmapped

  require_file "${out_dir}/collapsed_beta_matrix.txt"
}

run_qc_report() {
  local dataset="$1"
  local batch_root="$2"
  local final_beta="$3"
  local out_dir="$4"

  if ! is_enabled "${RUN_QC_REPORT}"; then
    echo "Skipping ${dataset} QC probe-loss report: RUN_QC_REPORT=0"
    return
  fi

  echo "========================================"
  echo "${dataset}: QC probe-loss and clock-overlap report"
  echo "========================================"

  "${PYTHON_BIN}" \
    "${WORK_ROOT}/scripts_postprocess/report_qc_probe_loss.py" \
    --dataset "${dataset}" \
    --batch-root "${batch_root}" \
    --final-beta "${final_beta}" \
    --out-dir "${out_dir}"

  require_file "${out_dir}/qc_probe_loss_summary.csv"
  require_file "${out_dir}/qc_affected_cpgs.tsv.gz"
  require_file "${out_dir}/qc_completely_lost_cpgs.tsv.gz"
  require_file "${out_dir}/clock_qc_overlap_summary.csv"
  require_file "${out_dir}/clock_qc_overlap_cpgs.tsv.gz"
}

run_biolearn_clocks() {
  local dataset="$1"
  local input_beta="$2"
  local out_dir="$3"

  if ! is_enabled "${RUN_BIOLEARN_CLOCKS}"; then
    echo "Skipping ${dataset} Biolearn clocks: RUN_BIOLEARN_CLOCKS=0"
    return
  fi

  echo "========================================"
  echo "${dataset}: Biolearn clocks with ${BIOLEARN_IMPUTATION_METHOD} imputation"
  echo "========================================"

  "${PYTHON_BIN}" -u \
    "${WORK_ROOT}/scripts_postprocess/run_clocks_biolearn_default_meta_grimagev1_lowmem.py" \
    --input "${input_beta}" \
    --meta "${WORK_ROOT}/meta.csv" \
    --out-dir "${out_dir}" \
    --biolearn-imputation-method "${BIOLEARN_IMPUTATION_METHOD}" \
    --allow-missing-age \
    --allow-missing-sex

  require_file "${out_dir}/clock.csv"
}

# ==============================
# Startup checks and settings log
# ==============================
echo "========================================"
echo "Methylation IDAT pipeline"
echo "BATCH_SIZE=${BATCH_SIZE}"
echo "THREADS=${THREADS}"
echo "PVAL_THRESHOLD=${PVAL_THRESHOLD}"
echo "BIOLEARN_IMPUTATION_METHOD=${BIOLEARN_IMPUTATION_METHOD}"
echo "RUN_EPICV2=${RUN_EPICV2}"
echo "RUN_MSA=${RUN_MSA}"
echo "RUN_QC_REPORT=${RUN_QC_REPORT}"
echo "RUN_BIOLEARN_CLOCKS=${RUN_BIOLEARN_CLOCKS}"
echo "RUN_FINAL_R_ANALYSIS=${RUN_FINAL_R_ANALYSIS}"
echo "SAVE_RAW_BETA=${SAVE_RAW_BETA}"
echo "SAVE_NONCOLLAPSED_BETA=${SAVE_NONCOLLAPSED_BETA}"
echo "SAVE_DETECTION_PVALS=${SAVE_DETECTION_PVALS}"
echo "SAVE_INTENSITY=${SAVE_INTENSITY}"
echo "WORK_ROOT=${WORK_ROOT}"
echo "SESAME_CACHE=${SESAME_CACHE}"
echo "MANIFEST_DIR=${MANIFEST_DIR}"
echo "========================================"

if ! is_enabled "${RUN_EPICV2}" && ! is_enabled "${RUN_MSA}"; then
  echo "ERROR: at least one of RUN_EPICV2 or RUN_MSA must be 1." >&2
  exit 1
fi

require_executable "${PYTHON_BIN}"
require_executable "${RSCRIPT_BIN}"
require_file "${WORK_ROOT}/bin/make_idat_batches.py"
require_file "${WORK_ROOT}/bin/run_sesame_batches_in_container.sh"
require_file "${WORK_ROOT}/scripts/sesame.R"
require_file "${WORK_ROOT}/scripts_postprocess/prepare_collapsed_beta.py"
require_file "${WORK_ROOT}/scripts_postprocess/report_qc_probe_loss.py"
require_dir "${SESAME_CACHE}"
require_dir "${MANIFEST_DIR}"

mkdir -p \
  "${WORK_ROOT}/results" \
  "${WORK_ROOT}/postprocess" \
  "${WORK_ROOT}/clock_compare_outputs_grimagev1" \
  "${WORK_ROOT}/results_grimagev1"

if is_enabled "${RUN_BIOLEARN_CLOCKS}"; then
  require_file "${WORK_ROOT}/meta.csv"
  require_file \
    "${WORK_ROOT}/scripts_postprocess/run_clocks_biolearn_default_meta_grimagev1_lowmem.py"
fi

# ==============================
# EPICv2
# ==============================
if is_enabled "${RUN_EPICV2}"; then
  require_dir "${WORK_ROOT}/idat_EPICv2"
  require_file "${WORK_ROOT}/name_table/EPICv2_blood.txt"
  require_file "${MANIFEST_DIR}/EPICv2_manifest.csv"

  echo "========================================"
  echo "EPICv2: make IDAT batches"
  echo "========================================"

  "${PYTHON_BIN}" "${WORK_ROOT}/bin/make_idat_batches.py" \
    --src "${WORK_ROOT}/idat_EPICv2" \
    --out "${WORK_ROOT}/idat_batches_EPICv2" \
    --batch-size "${BATCH_SIZE}"

  echo "========================================"
  echo "EPICv2: SeSAMe IDAT to per-batch beta/QC outputs"
  echo "========================================"

  PVAL_THRESHOLD="${PVAL_THRESHOLD}" \
  SESAME_CACHE="${SESAME_CACHE}" \
  MANIFEST_DIR="${MANIFEST_DIR}" \
  RSCRIPT_BIN="${RSCRIPT_BIN}" \
  SAVE_RAW_BETA="${SAVE_RAW_BETA}" \
  SAVE_NONCOLLAPSED_BETA="${SAVE_NONCOLLAPSED_BETA}" \
  SAVE_DETECTION_PVALS="${SAVE_DETECTION_PVALS}" \
  SAVE_INTENSITY="${SAVE_INTENSITY}" \
    "${WORK_ROOT}/bin/run_sesame_batches_in_container.sh" EPICv2

  combine_collapsed_beta \
    EPICv2 \
    "${WORK_ROOT}/results/EPICv2_batches" \
    "${WORK_ROOT}/name_table/EPICv2_blood.txt" \
    "${WORK_ROOT}/postprocess/EPICv2"

  run_qc_report \
    EPICv2 \
    "${WORK_ROOT}/results/EPICv2_batches" \
    "${WORK_ROOT}/postprocess/EPICv2/collapsed_beta_matrix.txt" \
    "${WORK_ROOT}/postprocess/EPICv2/qc_probe_report"

  run_biolearn_clocks \
    EPICv2 \
    "${WORK_ROOT}/postprocess/EPICv2/collapsed_beta_matrix.txt" \
    "${WORK_ROOT}/postprocess/EPICv2_biolearn_default_grimagev1"
fi

# ==============================
# MSA
# ==============================
if is_enabled "${RUN_MSA}"; then
  require_dir "${WORK_ROOT}/idat_MSA"
  require_file "${WORK_ROOT}/name_table/MSA_blood.txt"
  require_file "${MANIFEST_DIR}/MSA_manifest.csv"

  echo "========================================"
  echo "MSA: make IDAT batches"
  echo "========================================"

  "${PYTHON_BIN}" "${WORK_ROOT}/bin/make_idat_batches.py" \
    --src "${WORK_ROOT}/idat_MSA" \
    --out "${WORK_ROOT}/idat_batches_MSA" \
    --batch-size "${BATCH_SIZE}"

  echo "========================================"
  echo "MSA: SeSAMe IDAT to per-batch beta/QC outputs"
  echo "========================================"

  PVAL_THRESHOLD="${PVAL_THRESHOLD}" \
  SESAME_CACHE="${SESAME_CACHE}" \
  MANIFEST_DIR="${MANIFEST_DIR}" \
  RSCRIPT_BIN="${RSCRIPT_BIN}" \
  SAVE_RAW_BETA="${SAVE_RAW_BETA}" \
  SAVE_NONCOLLAPSED_BETA="${SAVE_NONCOLLAPSED_BETA}" \
  SAVE_DETECTION_PVALS="${SAVE_DETECTION_PVALS}" \
  SAVE_INTENSITY="${SAVE_INTENSITY}" \
    "${WORK_ROOT}/bin/run_sesame_batches_in_container.sh" MSA

  combine_collapsed_beta \
    MSA \
    "${WORK_ROOT}/results/MSA_batches" \
    "${WORK_ROOT}/name_table/MSA_blood.txt" \
    "${WORK_ROOT}/postprocess/MSA"

  run_qc_report \
    MSA \
    "${WORK_ROOT}/results/MSA_batches" \
    "${WORK_ROOT}/postprocess/MSA/collapsed_beta_matrix.txt" \
    "${WORK_ROOT}/postprocess/MSA/qc_probe_report"

  run_biolearn_clocks \
    MSA \
    "${WORK_ROOT}/postprocess/MSA/collapsed_beta_matrix.txt" \
    "${WORK_ROOT}/postprocess/MSA_biolearn_default_grimagev1"
fi

# ==============================
# EPICv2-vs-MSA downstream comparison
# ==============================
if is_enabled "${RUN_EPICV2}" \
  && is_enabled "${RUN_MSA}" \
  && is_enabled "${RUN_BIOLEARN_CLOCKS}"; then

  require_file "${WORK_ROOT}/templates/training.csv"
  require_file "${WORK_ROOT}/templates/testing.csv"
  require_file "${WORK_ROOT}/scripts_postprocess/create_clock_training_testing_files.py"

  echo "========================================"
  echo "Create 12 R-input files from EPICv2/MSA clock.csv"
  echo "========================================"

  "${PYTHON_BIN}" \
    "${WORK_ROOT}/scripts_postprocess/create_clock_training_testing_files.py" \
    --epicv2-clock \
      "${WORK_ROOT}/postprocess/EPICv2_biolearn_default_grimagev1/clock.csv" \
    --msa-clock \
      "${WORK_ROOT}/postprocess/MSA_biolearn_default_grimagev1/clock.csv" \
    --training-template "${WORK_ROOT}/templates/training.csv" \
    --testing-template "${WORK_ROOT}/templates/testing.csv" \
    --out-dir "${WORK_ROOT}/clock_compare_outputs_grimagev1" \
    --clocks Horvath Hannum PhenoAge GrimAgeV2 DunedinPACE GrimAgeV1

  if is_enabled "${RUN_FINAL_R_ANALYSIS}"; then
    require_file "${WORK_ROOT}/bin/run_final_r_analysis.sh"
    echo "========================================"
    echo "Run final R analysis script"
    echo "========================================"
    "${WORK_ROOT}/bin/run_final_r_analysis.sh"
  fi
else
  echo "========================================"
  echo "Skip EPICv2-vs-MSA comparison/final R analysis"
  echo "RUN_EPICV2=${RUN_EPICV2}, RUN_MSA=${RUN_MSA}, RUN_BIOLEARN_CLOCKS=${RUN_BIOLEARN_CLOCKS}"
  echo "========================================"
fi

echo "========================================"
echo "Pipeline finished."
echo "Core outputs:"
echo "  ${WORK_ROOT}/results"
echo "  ${WORK_ROOT}/postprocess"
echo "  ${WORK_ROOT}/clock_compare_outputs_grimagev1"
echo "  ${WORK_ROOT}/results_grimagev1"
echo "QC reports:"
if is_enabled "${RUN_EPICV2}"; then
  echo "  ${WORK_ROOT}/postprocess/EPICv2/qc_probe_report"
fi
if is_enabled "${RUN_MSA}"; then
  echo "  ${WORK_ROOT}/postprocess/MSA/qc_probe_report"
fi
echo "========================================"
