#!/usr/bin/env bash
set -euo pipefail

BATCH_SIZE="${BATCH_SIZE:-25}"
THREADS="${THREADS:-1}"
PVAL_THRESHOLD="${PVAL_THRESHOLD:-0.05}"
BIOLEARN_IMPUTATION_METHOD="${BIOLEARN_IMPUTATION_METHOD:-default}"

RUN_EPICV2="${RUN_EPICV2:-1}"
RUN_MSA="${RUN_MSA:-1}"
RUN_FINAL_R_ANALYSIS="${RUN_FINAL_R_ANALYSIS:-1}"

SESAME_CACHE="${SESAME_CACHE:-/opt/sesame_cache}"

echo "========================================"
echo "Methylation IDAT pipeline"
echo "BATCH_SIZE=${BATCH_SIZE}"
echo "THREADS=${THREADS}"
echo "PVAL_THRESHOLD=${PVAL_THRESHOLD}"
echo "BIOLEARN_IMPUTATION_METHOD=${BIOLEARN_IMPUTATION_METHOD}"
echo "RUN_EPICV2=${RUN_EPICV2}"
echo "RUN_MSA=${RUN_MSA}"
echo "RUN_FINAL_R_ANALYSIS=${RUN_FINAL_R_ANALYSIS}"
echo "========================================"

required_files=(
  /work/meta.csv
  /work/name_table/EPICv2_blood.txt
  /work/name_table/MSA_blood.txt
  /work/templates/training.csv
  /work/templates/testing.csv
)

for f in "${required_files[@]}"
do
  if [ ! -s "${f}" ]; then
    echo "ERROR: required file not found or empty: ${f}"
    exit 1
  fi
done

if [ "${RUN_EPICV2}" = "1" ]; then
  if [ ! -d /work/idat_EPICv2 ]; then
    echo "ERROR: /work/idat_EPICv2 not found. Please place EPICv2 IDAT files there."
    exit 1
  fi

  echo "========================================"
  echo "EPICv2: make IDAT batches"
  echo "========================================"

  /opt/pipeline/bin/make_idat_batches.py \
    --src /work/idat_EPICv2 \
    --out /work/idat_batches_EPICv2 \
    --batch-size "${BATCH_SIZE}"

  echo "========================================"
  echo "EPICv2: SeSAMe IDAT to beta"
  echo "========================================"

  THREADS="${THREADS}" \
  PVAL_THRESHOLD="${PVAL_THRESHOLD}" \
  SESAME_CACHE="${SESAME_CACHE}" \
    /opt/pipeline/bin/run_sesame_batches_in_container.sh EPICv2

  echo "========================================"
  echo "EPICv2: combine collapsed beta matrices"
  echo "========================================"

  /opt/conda/envs/methylation-beta-pipeline/bin/python /opt/pipeline/scripts_postprocess/prepare_collapsed_beta.py \
    --batch-root /work/results/EPICv2_batches \
    --sample-sheet /work/name_table/EPICv2_blood.txt \
    --out-dir /work/postprocess/EPICv2 \
    --output-name collapsed_beta_matrix.txt \
    --orientation cpg_rows \
    --drop-unmapped

  echo "========================================"
  echo "EPICv2: Biolearn clocks with default imputation"
  echo "========================================"

  /opt/conda/envs/methylation-beta-pipeline/bin/python -u /opt/pipeline/scripts_postprocess/run_clocks_biolearn_default_meta_grimagev1_lowmem.py \
    --input /work/postprocess/EPICv2/collapsed_beta_matrix.txt \
    --meta /work/meta.csv \
    --out-dir /work/postprocess/EPICv2_biolearn_default_grimagev1 \
    --biolearn-imputation-method "${BIOLEARN_IMPUTATION_METHOD}" \
    --allow-missing-age \
    --allow-missing-sex
fi

if [ "${RUN_MSA}" = "1" ]; then
  if [ ! -d /work/idat_MSA ]; then
    echo "ERROR: /work/idat_MSA not found. Please place MSA IDAT files there."
    exit 1
  fi

  echo "========================================"
  echo "MSA: make IDAT batches"
  echo "========================================"

  /opt/pipeline/bin/make_idat_batches.py \
    --src /work/idat_MSA \
    --out /work/idat_batches_MSA \
    --batch-size "${BATCH_SIZE}"

  echo "========================================"
  echo "MSA: SeSAMe IDAT to beta"
  echo "========================================"

  THREADS="${THREADS}" \
  PVAL_THRESHOLD="${PVAL_THRESHOLD}" \
  SESAME_CACHE="${SESAME_CACHE}" \
    /opt/pipeline/bin/run_sesame_batches_in_container.sh MSA

  echo "========================================"
  echo "MSA: combine collapsed beta matrices"
  echo "========================================"

  /opt/conda/envs/methylation-beta-pipeline/bin/python /opt/pipeline/scripts_postprocess/prepare_collapsed_beta.py \
    --batch-root /work/results/MSA_batches \
    --sample-sheet /work/name_table/MSA_blood.txt \
    --out-dir /work/postprocess/MSA \
    --output-name collapsed_beta_matrix.txt \
    --orientation cpg_rows \
    --drop-unmapped

  echo "========================================"
  echo "MSA: Biolearn clocks with default imputation"
  echo "========================================"

  /opt/conda/envs/methylation-beta-pipeline/bin/python -u /opt/pipeline/scripts_postprocess/run_clocks_biolearn_default_meta_grimagev1_lowmem.py \
    --input /work/postprocess/MSA/collapsed_beta_matrix.txt \
    --meta /work/meta.csv \
    --out-dir /work/postprocess/MSA_biolearn_default_grimagev1 \
    --biolearn-imputation-method "${BIOLEARN_IMPUTATION_METHOD}" \
    --allow-missing-age \
    --allow-missing-sex
fi

echo "========================================"
echo "Create 12 R-input files from EPICv2/MSA clock.csv"
echo "========================================"

/opt/conda/envs/methylation-beta-pipeline/bin/python /opt/pipeline/scripts_postprocess/create_clock_training_testing_files.py \
  --epicv2-clock /work/postprocess/EPICv2_biolearn_default_grimagev1/clock.csv \
  --msa-clock /work/postprocess/MSA_biolearn_default_grimagev1/clock.csv \
  --training-template /work/templates/training.csv \
  --testing-template /work/templates/testing.csv \
  --out-dir /work/clock_compare_outputs_grimagev1 \
  --clocks Horvath Hannum PhenoAge GrimAgeV2 DunedinPACE GrimAgeV1

if [ "${RUN_FINAL_R_ANALYSIS}" = "1" ]; then
  echo "========================================"
  echo "Run final R analysis script"
  echo "========================================"

  /opt/pipeline/bin/run_final_r_analysis.sh
fi

echo "========================================"
echo "Pipeline finished."
echo "Outputs:"
echo "  /work/results"
echo "  /work/postprocess"
echo "  /work/clock_compare_outputs_grimagev1"
echo "  /work/results_grimagev1"
echo "========================================"
