# methylation-beta-pipeline

Dockerized pipeline for processing Illumina DNA methylation IDAT files, generating suffix-collapsed beta matrices, calculating Biolearn epigenetic clocks, comparing MSA and EPICv2 results, and reporting CpGs affected by probe-level QC.

## Main workflow

For each enabled array, the pipeline performs:

1. Split raw IDAT pairs into batches.
2. Read IDATs with SeSAMe.
3. Generate pre-filtering and post-QCDPB beta values.
4. Collapse suffixed probe IDs to CpG prefixes with `collapseToPfx = TRUE`.
5. Record CpGs whose non-missing values were lost between unmasked and masked beta vectors derived from the same QCDPB-processed data.
6. Combine batch-level collapsed beta matrices.
7. Aggregate QC losses across batches and quantify overlap with CpGs required by six Biolearn clocks.
8. Optionally calculate clocks and run the EPICv2-vs-MSA comparison analysis.

The software environment is pinned by `environment.yml` and `conda-lock.yml` and is built for `linux/amd64`.

---

## Repository contents

Version-controlled files include:

```text
Dockerfile
environment.yml
conda-lock.yml
README.md
bin/
scripts/
scripts_postprocess/
templates/
package_versions/
```

Project-specific metadata files expected by the pipeline are:

```text
meta.csv
name_table/
├── EPICv2_blood.txt
└── MSA_blood.txt
```

The following large or generated resources are not stored in GitHub:

```text
manifest/EPICv2_manifest.csv
manifest/MSA_manifest.csv
vendor/sesame_cache/
idat_EPICv2/
idat_MSA/
idat_batches_EPICv2/
idat_batches_MSA/
results/
postprocess/
clock_compare_outputs_grimagev1/
results_grimagev1/
```

---

## Required external resources

Large manifest files and the SeSAMe/ExperimentHub cache are shared through Google Drive:

```text
https://drive.google.com/drive/folders/1Jh86-HM5R9FCyHvzW2t-efXkMYIYU1AO?usp=sharing
```

After downloading, the repository root must contain:

```text
manifest/
├── EPICv2_manifest.csv
└── MSA_manifest.csv

vendor/
└── sesame_cache/
    ├── BiocFileCache.sqlite
    └── ...
```

Raw IDAT files must be placed as follows:

```text
idat_EPICv2/
├── *_Red.idat
└── *_Grn.idat

idat_MSA/
├── *_Red.idat
└── *_Grn.idat
```

---

## Clone and select the QC-report branch

```bash
git clone https://github.com/tatsumashoji/methylation-beta-pipeline.git
cd methylation-beta-pipeline
git checkout paper_msa_epicv2_grimagev1
```

If the branch was created under a different name, replace `feature/qc-probe-loss-report` with the actual branch name.

---

## Build the Docker image

```bash
docker build \
  --platform linux/amd64 \
  -t methylation-idat-pipeline:grimagev1-biolearn-default .
```

Do not run `conda update`, `pip install`, `install.packages()`, or `BiocManager::install()` during pipeline execution. Package versions are fixed by `conda-lock.yml`.

---

## Run the complete pipeline

This runs EPICv2 and MSA from IDATs through beta-table generation, QC reporting, Biolearn clock calculation, creation of comparison inputs, and the final R analysis.

```bash
docker run --rm \
  --platform linux/amd64 \
  -e EXPERIMENT_HUB_CACHE=/work/vendor/sesame_cache \
  -e BFC_CACHE=/work/vendor/sesame_cache \
  -e SESAME_CACHE=/work/vendor/sesame_cache \
  -e MANIFEST_DIR=/work/manifest \
  -e OMP_NUM_THREADS=1 \
  -e OPENBLAS_NUM_THREADS=1 \
  -e MKL_NUM_THREADS=1 \
  -e VECLIB_MAXIMUM_THREADS=1 \
  -e NUMEXPR_NUM_THREADS=1 \
  -e RCPP_PARALLEL_NUM_THREADS=1 \
  -v "$PWD":/work \
  -w /work \
  methylation-idat-pipeline:grimagev1-biolearn-default \
  /work/bin/run_all_pipeline.sh
```

---

## Generate beta tables and QC reports only

These commands stop after generation of the integrated collapsed beta matrix and QC/clock-overlap reports. They do not calculate clock values or run the final comparison analysis.

### EPICv2 and MSA

```bash
docker run --rm \
  --platform linux/amd64 \
  -e RUN_EPICV2=1 \
  -e RUN_MSA=1 \
  -e RUN_QC_REPORT=1 \
  -e RUN_BIOLEARN_CLOCKS=0 \
  -e RUN_FINAL_R_ANALYSIS=0 \
  -e EXPERIMENT_HUB_CACHE=/work/vendor/sesame_cache \
  -e BFC_CACHE=/work/vendor/sesame_cache \
  -e SESAME_CACHE=/work/vendor/sesame_cache \
  -e MANIFEST_DIR=/work/manifest \
  -e OMP_NUM_THREADS=1 \
  -e OPENBLAS_NUM_THREADS=1 \
  -e MKL_NUM_THREADS=1 \
  -e VECLIB_MAXIMUM_THREADS=1 \
  -e NUMEXPR_NUM_THREADS=1 \
  -e RCPP_PARALLEL_NUM_THREADS=1 \
  -v "$PWD":/work \
  -w /work \
  methylation-idat-pipeline:grimagev1-biolearn-default \
  /work/bin/run_all_pipeline.sh
```

### EPICv2 only

```bash
docker run --rm \
  --platform linux/amd64 \
  -e RUN_EPICV2=1 \
  -e RUN_MSA=0 \
  -e RUN_QC_REPORT=1 \
  -e RUN_BIOLEARN_CLOCKS=0 \
  -e RUN_FINAL_R_ANALYSIS=0 \
  -e EXPERIMENT_HUB_CACHE=/work/vendor/sesame_cache \
  -e BFC_CACHE=/work/vendor/sesame_cache \
  -e SESAME_CACHE=/work/vendor/sesame_cache \
  -e MANIFEST_DIR=/work/manifest \
  -e OMP_NUM_THREADS=1 \
  -e OPENBLAS_NUM_THREADS=1 \
  -e MKL_NUM_THREADS=1 \
  -e VECLIB_MAXIMUM_THREADS=1 \
  -e NUMEXPR_NUM_THREADS=1 \
  -e RCPP_PARALLEL_NUM_THREADS=1 \
  -v "$PWD":/work \
  -w /work \
  methylation-idat-pipeline:grimagev1-biolearn-default \
  /work/bin/run_all_pipeline.sh
```

### MSA only

```bash
docker run --rm \
  --platform linux/amd64 \
  -e RUN_EPICV2=0 \
  -e RUN_MSA=1 \
  -e RUN_QC_REPORT=1 \
  -e RUN_BIOLEARN_CLOCKS=0 \
  -e RUN_FINAL_R_ANALYSIS=0 \
  -e EXPERIMENT_HUB_CACHE=/work/vendor/sesame_cache \
  -e BFC_CACHE=/work/vendor/sesame_cache \
  -e SESAME_CACHE=/work/vendor/sesame_cache \
  -e MANIFEST_DIR=/work/manifest \
  -e OMP_NUM_THREADS=1 \
  -e OPENBLAS_NUM_THREADS=1 \
  -e MKL_NUM_THREADS=1 \
  -e VECLIB_MAXIMUM_THREADS=1 \
  -e NUMEXPR_NUM_THREADS=1 \
  -e RCPP_PARALLEL_NUM_THREADS=1 \
  -v "$PWD":/work \
  -w /work \
  methylation-idat-pipeline:grimagev1-biolearn-default \
  /work/bin/run_all_pipeline.sh
```

---

## QC definitions

The QC report uses the following definitions in the suffix-collapsed CpG namespace:

- **Pre-QC value:** `getBetas(..., mask = FALSE, collapseToPfx = TRUE)` after applying the same SeSAMe `QCDPB` sequence and p-value threshold as the final matrix. The accumulated Q/P mask is ignored only during beta extraction.
- **Post-QC value:** `getBetas(..., mask = TRUE, collapseToPfx = TRUE)` from that same QCDPB-processed `SigDF`.
- **QC-affected CpG:** at least one sample had a non-missing pre-QC collapsed beta value but a missing post-QC collapsed beta value.
- **Completely lost CpG:** a QC-affected CpG with no non-missing post-QC value in any batch/sample for that array.
- **Partially affected CpG:** a QC-affected CpG that retains at least one non-missing post-QC value.

The broad `qc_affected` definition is appropriate for auditing sample-level data loss. The stricter `qc_completely_lost` definition identifies CpGs unavailable for downstream calculation without imputation.

---

## Batch-level QC outputs

For each batch, the pipeline keeps:

```text
results/<ARRAY>_batches/batch_XXX/
├── raw_beta_matrix.txt
├── collapsed_beta_matrix.txt
├── qc_summary.txt
├── qc_probe_failures.tsv.gz
├── pre_qc_collapsed_probe_ids.txt.gz
├── post_qc_collapsed_probe_ids.txt.gz
├── qc_batch_summary.csv
└── run.log
```

### `qc_probe_failures.tsv.gz`

One row per collapsed CpG affected in at least one sample in that batch. Important columns are:

```text
ProbeID
pre_qc_nonmissing_samples
post_qc_nonmissing_samples
qc_failed_samples
qc_failed_fraction_of_evaluable
qc_failed_fraction_of_all_samples
qc_failed_in_all_evaluable_samples
qc_completely_lost_in_batch
```

### `qc_batch_summary.csv`

One-row batch summary. Important columns are:

```text
n_qc_affected_cpgs
n_qc_completely_lost_cpgs
n_qc_partially_affected_cpgs
n_qc_failed_sample_cpg_pairs
```

---

## Integrated QC and clock-overlap reports

Reports are written to:

```text
postprocess/EPICv2/qc_probe_report/
postprocess/MSA/qc_probe_report/
```

Each directory contains:

```text
qc_probe_loss_summary.csv
qc_affected_cpgs.tsv.gz
qc_completely_lost_cpgs.tsv.gz
clock_qc_overlap_summary.csv
clock_qc_overlap_cpgs.tsv.gz
qc_report_metadata.json
```

### Where to find the CpG list and count lost by QC

#### Broad list: affected in one or more samples

```text
qc_affected_cpgs.tsv.gz
```

- `ProbeID` is the CpG name.
- The number of data rows is the number of distinct QC-affected CpGs.
- `qc_failed_samples` is the number of sample values lost for that CpG.
- `qc_completely_lost` indicates whether the CpG was lost from all samples.

The corresponding total count is also recorded in:

```text
qc_probe_loss_summary.csv
```

Column:

```text
qc_affected_cpg_n
```

#### Strict list: unavailable after QC in every sample

```text
qc_completely_lost_cpgs.tsv.gz
```

The corresponding total count is recorded in:

```text
qc_probe_loss_summary.csv
```

Column:

```text
qc_completely_lost_cpg_n
```

### Where to find overlap with each clock

The summary file is:

```text
clock_qc_overlap_summary.csv
```

It contains one row per clock. Important columns are:

```text
clock
biolearn_model
required_cpg_n
qc_affected_cpg_n
qc_affected_pct_of_required
qc_completely_lost_cpg_n
qc_completely_lost_pct_of_required
present_after_qc_n
missing_after_qc_n
present_in_integrated_matrix_n
absent_from_integrated_matrix_n
```

Therefore:

- `qc_affected_cpg_n` is how many required CpGs were affected in at least one sample.
- `qc_completely_lost_cpg_n` is how many required CpGs were unavailable in all samples after QC.
- The two percentage columns give these counts as percentages of all CpGs required by that clock.

The CpG-level detail is in:

```text
clock_qc_overlap_cpgs.tsv.gz
```

It contains every required CpG for every clock and the columns:

```text
ProbeID
present_before_qc
qc_affected
qc_completely_lost
qc_partially_affected
qc_failed_samples
present_after_qc
present_in_integrated_matrix
```

The following clocks are evaluated using the model definitions in the pinned Biolearn installation:

```text
Horvath / Horvathv1
Hannum
PhenoAge
GrimAgeV1
GrimAgeV2
DunedinPACE
```

The QC-affected and completely-lost CpG list files also include Boolean columns such as `required_by_Horvath`, plus `required_clock_count` and `required_clocks`.

---

## Main integrated beta matrices

```text
postprocess/EPICv2/collapsed_beta_matrix.txt
postprocess/MSA/collapsed_beta_matrix.txt
```

Rows are probe/CpG IDs and columns are renamed samples after batch integration.

---

## Existing results and rerunning

Old batch directories that contain only `raw_beta_matrix.txt` and `collapsed_beta_matrix.txt` do not contain enough information to reconstruct the QC-loss audit exactly. The batch runner skips a batch only when both beta matrices and all required QC files exist. Therefore, after adding this feature, old batches without QC files are automatically rerun once.

To force regeneration manually, remove the relevant batch outputs, for example:

```bash
rm -rf results/MSA_batches
rm -rf postprocess/MSA
```

Use this only when intentional, because SeSAMe processing can take substantial time.

---

## Expected output roots

```text
idat_batches_EPICv2/
idat_batches_MSA/
results/
postprocess/
clock_compare_outputs_grimagev1/
results_grimagev1/
```

---

## Important notes

- `vendor/sesame_cache/` and both manifest CSV files are required.
- Raw IDATs and generated outputs must not be committed to GitHub.
- QC reporting is limited to collapsed CpG IDs matching `cg` followed by digits; controls, `cgBK`, SNP and CH probes are not counted as clock CpGs.
- Clock-required CpG sets are obtained programmatically from `ModelGallery(...).get(..., imputation_method="none").methylation_sites()` in the pinned Biolearn environment rather than being hard-coded in this repository.
- `PVAL_THRESHOLD` defaults to `0.05` and can be changed with `-e PVAL_THRESHOLD=<value>`.
