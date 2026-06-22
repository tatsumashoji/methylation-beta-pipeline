# methylation-beta-pipeline

This repository contains a Dockerized pipeline for generating DNA methylation beta tables and running downstream MSA/EPICv2 comparison analyses.

The repository includes source code, Docker configuration, scripts, templates, and lightweight processed manifest files.
Large required files, including sesame cache files and raw Illumina manifest files, are distributed separately via Google Drive.

---

## Repository contents

The following files and directories are included in this GitHub repository:

```text
Dockerfile
environment.yml
README.md
bin/
scripts/
scripts_postprocess/
templates/
meta.csv
name_table/
```

The following files and directories are not included in GitHub:

```text
vendor/sesame_cache/
manifest/EPICv2_manifest.csv
manifest/MSA_manifest.csv
idat_EPICv2/
idat_MSA/
idat_batches_EPICv2/
idat_batches_MSA/
results/
postprocess/
results_grimagev1/
clock_compare_outputs_grimagev1/
```

These files are excluded because they are large, generated, cache-derived, or may contain sample-related information.

---

## Required external resources

The following external resource archive is required to run the pipeline:

```text
manifest
vendor
```

This file is shared via Google Drive.

Google Drive URL:

```text
https://drive.google.com/drive/folders/1Jh86-HM5R9FCyHvzW2t-efXkMYIYU1AO?usp=sharing
```

The archive should contain:

```text
manifest/
├── EPICv2_manifest.csv
└── MSA_manifest.csv

vendor/
└── sesame_cache/
    ├── BiocFileCache.sqlite
    ├── BiocFileCache.sqlite.LOCK
    └── ...
```

The raw IDAT files and sample metadata are shared separately and should be placed as follows:

```text
idat_EPICv2/
idat_MSA/
```

---

## Setup

Clone this repository:

```bash
git clone https://github.com/tatsumashoji/methylation-beta-pipeline.git
cd methylation-beta-pipeline
```

Switch to the manuscript branch:

```bash
git checkout paper_msa_epicv2_grimagev1
```

Download the external resource archive from Google Drive and place the downloaded archive in the repository root:

```text
methylation-beta-pipeline/
├── manifest
├── vendor
├── Dockerfile
├── bin/
├── scripts/
└── ...
```

Then place the IDAT files and metadata in the following locations:

```text
idat_EPICv2/
├── *_Red.idat
└── *_Grn.idat

idat_MSA/
├── *_Red.idat
└── *_Grn.idat
```

Create output directories:

```bash
mkdir -p idat_batches_EPICv2 idat_batches_MSA
mkdir -p results postprocess results_grimagev1 clock_compare_outputs_grimagev1
```

---

## Build Docker image

Build the Docker image:

```bash
docker build \
  --platform linux/amd64 \
  -t methylation-idat-pipeline:grimagev1-biolearn-default .
```

---

## Run full pipeline

To run the full pipeline, use:

```bash
docker run --rm \
  --platform linux/amd64 \
  -e EXPERIMENT_HUB_CACHE=/work/vendor/sesame_cache \
  -v "$PWD":/work \
  -w /work \
  methylation-idat-pipeline:grimagev1-biolearn-default \
  /work/bin/run_all_pipeline.sh
```

---

## Expected outputs

The pipeline writes outputs to:

```text
idat_batches_EPICv2/
idat_batches_MSA/
results/
postprocess/
results_grimagev1/
clock_compare_outputs_grimagev1/
```

---

## Important notes

* `vendor/sesame_cache/` is required for reproducibility and stable execution.
* Raw Illumina manifest files are required by the methylation array preprocessing scripts.
* These large files are distributed via Google Drive instead of GitHub.
* Do not commit `vendor/sesame_cache/`, raw IDAT files, raw manifest files, metadata, or output directories to GitHub.
* If `vendor/sesame_cache/` is missing or incomplete, sesame/sesameData may attempt to download required resources again.
* The cache location is explicitly set by:

```bash
-e EXPERIMENT_HUB_CACHE=/work/vendor/sesame_cache
```

