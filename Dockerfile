FROM --platform=linux/amd64 condaforge/miniforge3:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Tokyo
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

WORKDIR /opt/pipeline

# Minimal OS packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    git \
    bash \
    && rm -rf /var/lib/apt/lists/*

# Copy dependency definitions
COPY environment.yml /opt/pipeline/environment.yml
COPY conda-lock.yml /opt/pipeline/conda-lock.yml

# Create the locked conda environment
RUN conda install -y -n base -c conda-forge conda-lock && \
    conda-lock install --name methylation-beta-pipeline /opt/pipeline/conda-lock.yml && \
    conda clean -afy

# Use the locked environment by default
ENV PATH="/opt/conda/bin:/opt/conda/envs/methylation-beta-pipeline/bin:${PATH}"

# Copy only code and lightweight files into the image.
# Large manifest CSV files and sesame_cache are mounted at runtime from /work.
COPY bin/ /opt/pipeline/bin/
COPY scripts/ /opt/pipeline/scripts/
COPY scripts_postprocess/ /opt/pipeline/scripts_postprocess/
COPY templates/ /opt/pipeline/templates/
COPY README.md /opt/pipeline/README.md

RUN chmod +x /opt/pipeline/bin/* || true && \
    chmod +x /opt/pipeline/scripts_postprocess/*.py || true

# Runtime resource paths.
# The repository root is mounted to /work at docker run time.
ENV EXPERIMENT_HUB_CACHE=/work/vendor/sesame_cache

# Sanity checks for package availability
RUN /opt/conda/envs/methylation-beta-pipeline/bin/Rscript -e 'suppressPackageStartupMessages({library(sesame); library(sesameData); library(ExperimentHub); library(data.table); library(dplyr); library(RcppTOML); library(this.path)}); cat("R package check OK\n")'

RUN /opt/conda/envs/methylation-beta-pipeline/bin/python - <<'PY'
import pandas
import numpy
import scipy
import sklearn
import seaborn
from biolearn.data_library import GeoData
from biolearn.model_gallery import ModelGallery
print("Python/Biolearn check OK")
PY

ENTRYPOINT ["/bin/bash"]
