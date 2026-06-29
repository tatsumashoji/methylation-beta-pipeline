#!/usr/bin/env bash
set -euo pipefail

# Force serial execution to avoid nested-thread creation failures in Docker.
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

PVAL_THRESHOLD="${PVAL_THRESHOLD:-0.05}"
SESAME_CACHE="${SESAME_CACHE:-/work/vendor/sesame_cache}"
MANIFEST_DIR="${MANIFEST_DIR:-/work/manifest}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/opt/conda/envs/methylation-beta-pipeline/bin/Rscript}"
SESAME_SCRIPT="${SESAME_SCRIPT:-/work/scripts/sesame.R}"

# Large optional outputs are disabled by default to avoid filling the host disk.
# Set to 1 only when these intermediate matrices are explicitly required.
SAVE_RAW_BETA="${SAVE_RAW_BETA:-0}"
SAVE_NONCOLLAPSED_BETA="${SAVE_NONCOLLAPSED_BETA:-0}"
SAVE_DETECTION_PVALS="${SAVE_DETECTION_PVALS:-0}"
SAVE_INTENSITY="${SAVE_INTENSITY:-0}"

normalize_bool() {
  case "${1}" in
    1|true|TRUE|yes|YES|y|Y) echo "true" ;;
    0|false|FALSE|no|NO|n|N) echo "false" ;;
    *)
      echo "ERROR: invalid boolean value: ${1}" >&2
      exit 2
      ;;
  esac
}

SAVE_RAW_BETA_BOOL="$(normalize_bool "${SAVE_RAW_BETA}")"
SAVE_NONCOLLAPSED_BETA_BOOL="$(normalize_bool "${SAVE_NONCOLLAPSED_BETA}")"
SAVE_DETECTION_PVALS_BOOL="$(normalize_bool "${SAVE_DETECTION_PVALS}")"
SAVE_INTENSITY_BOOL="$(normalize_bool "${SAVE_INTENSITY}")"

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
    echo "ERROR: DATASET must be EPICv2 or MSA" >&2
    exit 2
    ;;
esac

required_paths=(
  "${BATCH_ROOT}"
  "${SESAME_CACHE}"
  "${MANIFEST_FILE}"
  "${RSCRIPT_BIN}"
  "${SESAME_SCRIPT}"
)

for path in "${required_paths[@]}"; do
  if [ ! -e "${path}" ]; then
    echo "ERROR: required path not found: ${path}" >&2
    exit 1
  fi
done

mkdir -p "${RESULT_ROOT}"

mapfile -t batches < <(
  find "${BATCH_ROOT}" -maxdepth 1 -type d -name 'batch_*' -print | sort
)

if [ "${#batches[@]}" -eq 0 ]; then
  echo "ERROR: no batch_* directories found under ${BATCH_ROOT}" >&2
  exit 1
fi

required_outputs=(
  collapsed_beta_matrix.txt
  qc_summary.txt
  qc_probe_failures.tsv.gz
  pre_qc_collapsed_probe_ids.txt.gz
  post_qc_collapsed_probe_ids.txt.gz
  qc_batch_summary.csv
)

if [ "${SAVE_RAW_BETA_BOOL}" = "true" ]; then
  required_outputs+=(raw_beta_matrix.txt)
fi
if [ "${SAVE_NONCOLLAPSED_BETA_BOOL}" = "true" ]; then
  required_outputs+=(beta_matrix.txt)
fi
if [ "${SAVE_DETECTION_PVALS_BOOL}" = "true" ]; then
  required_outputs+=(detection_pvals.txt)
fi
if [ "${SAVE_INTENSITY_BOOL}" = "true" ]; then
  required_outputs+=(methylated_intensity.txt unmethylated_intensity.txt)
fi

for batch_dir in "${batches[@]}"; do
  batch_name="$(basename "${batch_dir}")"
  out_dir="${RESULT_ROOT}/${batch_name}"
  log_file="${out_dir}/run.log"

  echo "=============================="
  echo "Running ${DATASET} ${batch_name}"
  echo "Input: ${batch_dir}"
  echo "Output: ${out_dir}"
  echo "SAVE_RAW_BETA=${SAVE_RAW_BETA}"
  echo "SAVE_NONCOLLAPSED_BETA=${SAVE_NONCOLLAPSED_BETA}"
  echo "SAVE_DETECTION_PVALS=${SAVE_DETECTION_PVALS}"
  echo "SAVE_INTENSITY=${SAVE_INTENSITY}"
  echo "=============================="

  mkdir -p "${out_dir}"

  all_outputs_exist=1
  for output_name in "${required_outputs[@]}"; do
    if [ ! -s "${out_dir}/${output_name}" ]; then
      all_outputs_exist=0
      break
    fi
  done

  if [ "${all_outputs_exist}" -eq 1 ]; then
    echo "Skipping ${batch_name}: beta matrices and QC reports already exist."
    continue
  fi

  red_count="$(
    find "${batch_dir}" -maxdepth 1 -type f ! -name '._*' -name '*_Red.idat' \
      | wc -l | tr -d '[:space:]'
  )"
  grn_count="$(
    find "${batch_dir}" -maxdepth 1 -type f ! -name '._*' -name '*_Grn.idat' \
      | wc -l | tr -d '[:space:]'
  )"

  if [ "${red_count}" -eq 0 ] || [ "${red_count}" -ne "${grn_count}" ]; then
    echo "ERROR: invalid IDAT count in ${batch_dir}: Red=${red_count}, Grn=${grn_count}" >&2
    exit 1
  fi

  set +e
  "${RSCRIPT_BIN}" "${SESAME_SCRIPT}" \
    --idat_dir "${batch_dir}" \
    --out_path "${out_dir}" \
    --array_type "${ARRAY_TYPE}" \
    --manifest_file "${MANIFEST_FILE}" \
    --pval_threshold "${PVAL_THRESHOLD}" \
    --threads 1 \
    --sesame_data "${SESAME_CACHE}" \
    --skip_annotation true \
    --save_raw_beta "${SAVE_RAW_BETA_BOOL}" \
    --save_noncollapsed_beta "${SAVE_NONCOLLAPSED_BETA_BOOL}" \
    --save_detection_pvals "${SAVE_DETECTION_PVALS_BOOL}" \
    --save_intensity "${SAVE_INTENSITY_BOOL}" \
    2>&1 | tee "${log_file}"
  status=${PIPESTATUS[0]}
  set -e

  if [ "${status}" -ne 0 ]; then
    echo "ERROR: ${DATASET} ${batch_name} failed. See ${log_file}" >&2
    exit "${status}"
  fi

  for output_name in "${required_outputs[@]}"; do
    if [ ! -s "${out_dir}/${output_name}" ]; then
      echo "ERROR: ${out_dir}/${output_name} was not created." >&2
      exit 1
    fi
  done

  # Remove optional intermediates unless explicitly requested. They can be
  # regenerated from IDATs and can be very large.
  if [ "${SAVE_RAW_BETA_BOOL}" != "true" ]; then
    rm -f "${out_dir}/raw_beta_matrix.txt"
  fi
  if [ "${SAVE_NONCOLLAPSED_BETA_BOOL}" != "true" ]; then
    rm -f "${out_dir}/beta_matrix.txt" "${out_dir}/beta_matrix_ann.txt"
  fi
  if [ "${SAVE_DETECTION_PVALS_BOOL}" != "true" ]; then
    rm -f "${out_dir}/detection_pvals.txt"
  fi
  if [ "${SAVE_INTENSITY_BOOL}" != "true" ]; then
    rm -f "${out_dir}/methylated_intensity.txt" "${out_dir}/unmethylated_intensity.txt"
  fi
  rm -f "${out_dir}/collapsed_beta_matrix_ann.txt"

  echo "OK: ${DATASET} ${batch_name}"
done
