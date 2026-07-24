#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/work}"
CONDA_ENV="${CONDA_ENV:-/opt/conda/envs/methylation-beta-pipeline}"
RSCRIPT_BIN="${CONDA_ENV}/bin/Rscript"

INPUT_DIR="${INPUT_DIR:-${WORK_ROOT}/clock_compare_outputs_grimagev1}"
OUT_DIR="${OUT_DIR:-${WORK_ROOT}/results_grimagev1}"
SCRIPT="${SCRIPT:-${WORK_ROOT}/scripts_postprocess/epigenetic_clock_msa_vs_epicv2_pipeline.R}"

mkdir -p "${OUT_DIR}/figures" "${OUT_DIR}/tables" "${OUT_DIR}/logs"

echo "========================================"
echo "Final R analysis: MSA vs EPICv2 clocks"
echo "Script:    ${SCRIPT}"
echo "Input dir: ${INPUT_DIR}"
echo "Out dir:   ${OUT_DIR}"
echo "========================================"

if [ ! -s "${SCRIPT}" ]; then
  echo "ERROR: R script not found: ${SCRIPT}"
  exit 1
fi

required_files=(
  "${INPUT_DIR}/Horvath.training.csv"
  "${INPUT_DIR}/Horvath.testing.csv"
  "${INPUT_DIR}/Hannum.training.csv"
  "${INPUT_DIR}/Hannum.testing.csv"
  "${INPUT_DIR}/PhenoAge.training.csv"
  "${INPUT_DIR}/PhenoAge.testing.csv"
  "${INPUT_DIR}/GrimAgeV1.training.csv"
  "${INPUT_DIR}/GrimAgeV1.testing.csv"
  "${INPUT_DIR}/GrimAgeV2.training.csv"
  "${INPUT_DIR}/GrimAgeV2.testing.csv"
  "${INPUT_DIR}/DunedinPACE.training.csv"
  "${INPUT_DIR}/DunedinPACE.testing.csv"
)

for f in "${required_files[@]}"; do
  if [ ! -s "${f}" ]; then
    echo "ERROR: required input file not found or empty: ${f}"
    exit 1
  fi
done

cd "${WORK_ROOT}"

"${RSCRIPT_BIN}" "${SCRIPT}" \
  2>&1 | tee "${OUT_DIR}/logs/epigenetic_clock_msa_vs_epicv2.log"

echo "========================================"
echo "Final R analysis finished."
echo "Outputs:"
echo "  ${OUT_DIR}/figures"
echo "  ${OUT_DIR}/tables"
echo "  ${OUT_DIR}/logs"
echo "========================================"
