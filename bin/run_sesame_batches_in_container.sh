#!/usr/bin/env bash
set -euo pipefail

# Force serial execution to avoid Docker/BiocParallel pthread issues.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1
export THREADS=1
export N_THREADS=1
export NUM_THREADS=1

DATASET="${1:?Usage: run_sesame_batches_in_container.sh EPICv2|MSA}"
THREADS=1
PVAL_THRESHOLD="${PVAL_THRESHOLD:-0.05}"
SESAME_CACHE="${SESAME_CACHE:-/work/vendor/sesame_cache}"
MANIFEST_DIR="${MANIFEST_DIR:-/work/manifest}"
CONDA_ENV="${CONDA_ENV:-/opt/conda/envs/methylation-beta-pipeline}"
RSCRIPT_BIN="${CONDA_ENV}/bin/Rscript"

SAVE_RAW_BETA="${SAVE_RAW_BETA:-0}"
SAVE_NONCOLLAPSED_BETA="${SAVE_NONCOLLAPSED_BETA:-0}"
SAVE_DETECTION_PVALS="${SAVE_DETECTION_PVALS:-0}"
SAVE_INTENSITY="${SAVE_INTENSITY:-0}"

case "${DATASET}" in
  EPICv2)
    ARRAY_TYPE="EPICv2"
    BATCH_ROOT="/work/idat_batches_EPICv2"
    RESULT_ROOT="/work/results/EPICv2_batches"
    MANIFEST_FILE="${MANIFEST_DIR}/EPICv2_manifest.csv"
    ;;
  MSA)
    ARRAY_TYPE="MSA"
    BATCH_ROOT="/work/idat_batches_MSA"
    RESULT_ROOT="/work/results/MSA_batches"
    MANIFEST_FILE="${MANIFEST_DIR}/MSA_manifest.csv"
    ;;
  *)
    echo "ERROR: DATASET must be EPICv2 or MSA"
    exit 1
    ;;
esac

if [ ! -d "${BATCH_ROOT}" ]; then
  echo "ERROR: batch root not found: ${BATCH_ROOT}"
  exit 1
fi

mkdir -p "${RESULT_ROOT}"

for B in $(find "${BATCH_ROOT}" -maxdepth 1 -type d -name 'batch_*' | sort)
do
  BN=$(basename "${B}")
  OUT="${RESULT_ROOT}/${BN}"
  LOG="${OUT}/run.log"

  echo "=============================="
  echo "Running ${DATASET} ${BN}"
  echo "Input: ${B}"
  echo "Output: ${OUT}"
  echo "SAVE_RAW_BETA=${SAVE_RAW_BETA}"
  echo "SAVE_NONCOLLAPSED_BETA=${SAVE_NONCOLLAPSED_BETA}"
  echo "SAVE_DETECTION_PVALS=${SAVE_DETECTION_PVALS}"
  echo "SAVE_INTENSITY=${SAVE_INTENSITY}"
  echo "=============================="

  mkdir -p "${OUT}"

  if [ -s "${OUT}/collapsed_beta_matrix.txt" ] && \
     [ -s "${OUT}/qc_probe_failures_by_sample.tsv.gz" ] && \
     [ -s "${OUT}/qc_probe_failures.tsv.gz" ] && \
     [ -s "${OUT}/pre_qc_collapsed_probe_ids.txt.gz" ] && \
     [ -s "${OUT}/post_qc_collapsed_probe_ids.txt.gz" ] && \
     [ -s "${OUT}/qc_batch_summary.csv" ]; then
    echo "Skipping ${BN}: collapsed beta matrix and QC reports already exist."
    continue
  fi

  RED_COUNT=$(find "${B}" -maxdepth 1 -type f ! -name '._*' -name '*_Red.idat' | wc -l | tr -d '[:space:]')
  GRN_COUNT=$(find "${B}" -maxdepth 1 -type f ! -name '._*' -name '*_Grn.idat' | wc -l | tr -d '[:space:]')

  if [ "${RED_COUNT}" -eq 0 ] || [ "${RED_COUNT}" -ne "${GRN_COUNT}" ]; then
    echo "ERROR: invalid IDAT count in ${B}: Red=${RED_COUNT}, Grn=${GRN_COUNT}"
    exit 1
  fi

  set +e
  "${RSCRIPT_BIN}" /work/scripts/sesame.R \
    --idat_dir "${B}" \
    --out_path "${OUT}" \
    --array_type "${ARRAY_TYPE}" \
    --manifest_file "${MANIFEST_FILE}" \
    --pval_threshold "${PVAL_THRESHOLD}" \
    --threads 1 \
    --sesame_data "${SESAME_CACHE}" \
    --skip_annotation true \
    --save_raw_beta "${SAVE_RAW_BETA}" \
    --save_noncollapsed_beta "${SAVE_NONCOLLAPSED_BETA}" \
    --save_detection_pvals "${SAVE_DETECTION_PVALS}" \
    --save_intensity "${SAVE_INTENSITY}" \
    2>&1 | tee "${LOG}"

  STATUS=${PIPESTATUS[0]}
  set -e

  if [ "${STATUS}" -ne 0 ]; then
    echo "ERROR: ${DATASET} ${BN} failed. See ${LOG}"
    exit "${STATUS}"
  fi

  for REQUIRED in \
    collapsed_beta_matrix.txt \
    qc_probe_failures_by_sample.tsv.gz \
    qc_probe_failures.tsv.gz \
    pre_qc_collapsed_probe_ids.txt.gz \
    post_qc_collapsed_probe_ids.txt.gz \
    qc_batch_summary.csv
  do
    if [ ! -s "${OUT}/${REQUIRED}" ]; then
      echo "ERROR: ${OUT}/${REQUIRED} was not created."
      exit 1
    fi
  done

  # Remove only temporary files. Do not delete QC reports or run.log.
  find "${OUT}" -maxdepth 1 -type f \( -name '*.tmp' -o -name '*.tmp.*' \) -delete

  echo "OK: ${DATASET} ${BN}"
done
