#!/usr/bin/env bash

# --- Force serial execution in Docker to avoid pthread_create / BiocParallel errors ---
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1
export THREADS=1
export N_THREADS=1
export NUM_THREADS=1
# --- End serial execution patch ---

set -euo pipefail

DATASET="${1:?Usage: run_sesame_batches_in_container.sh EPICv2|MSA}"

THREADS=1
PVAL_THRESHOLD="${PVAL_THRESHOLD:-0.05}"
SESAME_CACHE="${SESAME_CACHE:-/work/vendor/sesame_cache}"

MANIFEST_DIR="${MANIFEST_DIR:-/work/manifest}"

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
  echo "=============================="

  mkdir -p "${OUT}"

  if [ -s "${OUT}/collapsed_beta_matrix.txt" ] && [ -s "${OUT}/raw_beta_matrix.txt" ]; then
    echo "Skipping ${BN}: required outputs already exist."
    continue
  fi

  RED_COUNT=$(find "${B}" -maxdepth 1 -type f ! -name '._*' -name '*_Red.idat' | wc -l | tr -d '[:space:]')
  GRN_COUNT=$(find "${B}" -maxdepth 1 -type f ! -name '._*' -name '*_Grn.idat' | wc -l | tr -d '[:space:]')

  if [ "${RED_COUNT}" -eq 0 ] || [ "${RED_COUNT}" -ne "${GRN_COUNT}" ]; then
    echo "ERROR: invalid IDAT count in ${B}: Red=${RED_COUNT}, Grn=${GRN_COUNT}"
    exit 1
  fi

  set +e
  /opt/conda/envs/methylation-beta-pipeline/bin/Rscript /opt/pipeline/scripts/sesame.R \
    --idat_dir "${B}" \
    --out_path "${OUT}" \
    --array_type "${ARRAY_TYPE}" \
    --manifest_file "${MANIFEST_FILE}" \
    --pval_threshold "${PVAL_THRESHOLD}" \
    --threads 1 \
    --sesame_data "${SESAME_CACHE}" \
    --skip_annotation true \
    2>&1 | tee "${LOG}"

  STATUS=${PIPESTATUS[0]}
  set -e

  if [ "${STATUS}" -ne 0 ]; then
    echo "ERROR: ${DATASET} ${BN} failed. See ${LOG}"
    exit "${STATUS}"
  fi

  if [ ! -s "${OUT}/collapsed_beta_matrix.txt" ]; then
    echo "ERROR: ${OUT}/collapsed_beta_matrix.txt was not created."
    exit 1
  fi

  if [ ! -s "${OUT}/raw_beta_matrix.txt" ]; then
    echo "ERROR: ${OUT}/raw_beta_matrix.txt was not created."
    exit 1
  fi

  find "${OUT}" \
    -type f \
    ! \( -name "collapsed_beta_matrix.txt" -o -name "raw_beta_matrix.txt" \) \
    -delete

  echo "OK: ${DATASET} ${BN}"
done
