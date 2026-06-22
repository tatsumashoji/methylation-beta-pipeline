#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="${INPUT_DIR:-/work/clock_compare_outputs_grimagev1}"
OUT_DIR="${OUT_DIR:-/work/results_grimagev1}"
R_SCRIPT="${R_SCRIPT:-/opt/pipeline/scripts_postprocess/epigenetic_clock_msa_vs_epicv2_pipeline.R}"

CLOCKS=(
  Horvath
  Hannum
  PhenoAge
  GrimAgeV1
  GrimAgeV2
  DunedinPACE
)

echo "========================================"
echo "Run final R analysis"
echo "R script : ${R_SCRIPT}"
echo "Input dir: ${INPUT_DIR}"
echo "Output   : ${OUT_DIR}"
echo "Clocks   : ${CLOCKS[*]}"
echo "========================================"

if [ ! -s "${R_SCRIPT}" ]; then
  echo "ERROR: R script not found: ${R_SCRIPT}"
  exit 1
fi

for CLOCK in "${CLOCKS[@]}"
do
  TRAIN="${INPUT_DIR}/${CLOCK}.training.csv"
  TEST="${INPUT_DIR}/${CLOCK}.testing.csv"

  if [ ! -s "${TRAIN}" ]; then
    echo "ERROR: training file not found or empty: ${TRAIN}"
    exit 1
  fi

  if [ ! -s "${TEST}" ]; then
    echo "ERROR: testing file not found or empty: ${TEST}"
    exit 1
  fi
done

mkdir -p "${OUT_DIR}/logs"

/opt/conda/envs/methylation-beta-pipeline/bin/Rscript -e '
pkgs <- c("tidyverse", "ggplot2", "knitr", "gridExtra", "cowplot")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing R packages: ", paste(missing, collapse = ", "))
}
cat("R package check OK\n")
'

LOG="${OUT_DIR}/logs/epigenetic_clock_msa_vs_epicv2.log"

/opt/conda/envs/methylation-beta-pipeline/bin/Rscript "${R_SCRIPT}" "${INPUT_DIR}" "${OUT_DIR}" 2>&1 | tee "${LOG}"

echo "========================================"
echo "Final R analysis finished."
echo "Log: ${LOG}"
echo "Figures:"
ls -lh "${OUT_DIR}/figures" || true
echo "Tables:"
ls -lh "${OUT_DIR}/tables" || true
echo "========================================"
