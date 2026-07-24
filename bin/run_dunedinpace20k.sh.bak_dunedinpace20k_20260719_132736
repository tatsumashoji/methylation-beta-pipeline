#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/work}"
CONDA_ENV="${CONDA_ENV:-/opt/conda/envs/methylation-beta-pipeline}"
PYTHON_BIN="${CONDA_ENV}/bin/python"
BIOLEARN_IMPUTATION_METHOD="${BIOLEARN_IMPUTATION_METHOD:-default}"
RUN_EPICV2="${RUN_EPICV2:-1}"
RUN_MSA="${RUN_MSA:-1}"
REPLACE_EXISTING_DUNEDINPACE="${REPLACE_EXISTING_DUNEDINPACE:-1}"

run_one() {
  local dataset="$1"
  local input="${WORK_ROOT}/postprocess/${dataset}/collapsed_beta_matrix.txt"
  local out_dir="${WORK_ROOT}/postprocess/${dataset}_biolearn_default_grimagev1"
  local clock_csv="${out_dir}/clock.csv"

  echo "========================================"
  echo "${dataset}: DunedinPACE with 20k gold/background probes"
  echo "========================================"
  echo "Input:     ${input}"
  echo "Out dir:   ${out_dir}"
  echo "Clock CSV: ${clock_csv}"

  if [ ! -s "${input}" ]; then
    echo "ERROR: input collapsed beta matrix not found: ${input}" >&2
    exit 1
  fi

  mkdir -p "${out_dir}"

  local replace_flag=""
  if [ "${REPLACE_EXISTING_DUNEDINPACE}" = "1" ] && [ -s "${clock_csv}" ]; then
    replace_flag="--replace-existing"
  fi

  local clock_args=()
  if [ -s "${clock_csv}" ]; then
    clock_args=(--clock-csv "${clock_csv}")
  fi

  "${PYTHON_BIN}" -u "${WORK_ROOT}/scripts_postprocess/run_dunedinpace_20k_biolearn.py" \
    --input "${input}" \
    --meta "${WORK_ROOT}/meta.csv" \
    --out-dir "${out_dir}" \
    --imputation-method "${BIOLEARN_IMPUTATION_METHOD}" \
    "${clock_args[@]}" \
    ${replace_flag}
}

if [ "${RUN_EPICV2}" = "1" ]; then
  run_one EPICv2
fi

if [ "${RUN_MSA}" = "1" ]; then
  run_one MSA
fi

echo "========================================"
echo "DunedinPACE 20k recalculation finished."
echo "========================================"
