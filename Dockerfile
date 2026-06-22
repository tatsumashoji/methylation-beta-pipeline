FROM bioconductor/bioconductor_docker:RELEASE_3_23

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Tokyo
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

WORKDIR /opt/pipeline

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    git \
    python3 \
    python3-pip \
    python3-venv \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# R packages: SeSAMe + final R analysis script dependencies.
RUN R -q -e '\
options(timeout = 1200); \
options(download.file.method = "libcurl"); \
options(repos = c(CRAN = "https://cloud.r-project.org")); \
BiocManager::install( \
  c("sesame", "sesameData", "BiocParallel", "ExperimentHub"), \
  ask = FALSE, update = FALSE \
); \
install.packages( \
  c( \
    "logger", "data.table", "dplyr", "magrittr", "this.path", "tidyr", "RcppTOML", \
    "tidyverse", "ggplot2", "knitr", "gridExtra", "cowplot", "broom", "scales" \
  ), \
  repos = "https://cloud.r-project.org" \
); \
'

# Python environment for post-processing and Biolearn.
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

RUN pip install --upgrade pip setuptools wheel && \
    pip install pandas numpy scipy scikit-learn threadpoolctl matplotlib seaborn biolearn && \
    pip install --index-url https://download.pytorch.org/whl/cpu torch torchvision torchaudio

COPY scripts/ /opt/pipeline/scripts/
COPY scripts_postprocess/ /opt/pipeline/scripts_postprocess/
COPY bin/ /opt/pipeline/bin/
COPY manifest/ /opt/pipeline/manifest/
COPY vendor/sesame_cache/ /opt/sesame_cache/

ENV EXPERIMENT_HUB_CACHE=/opt/sesame_cache

RUN chmod +x /opt/pipeline/bin/* || true && \
    chmod +x /opt/pipeline/scripts_postprocess/*.py || true

# Check SeSAMe cache.
RUN R -q -e '\
suppressPackageStartupMessages({library(ExperimentHub); library(sesameData)}); \
ExperimentHub::setExperimentHubOption("CACHE", "/opt/sesame_cache"); \
for (title in c("idatSignature", "MSA.address", "KYCG.MSA.Mask.20260122")) { \
  cat("Checking sesameData cache:", title, "\n"); \
  sesameData::sesameDataGet(title); \
}; \
cat("SeSAMe cache check OK\n"); \
'

# Check Python/Biolearn.
RUN python - <<'PY'
import pandas, numpy, sklearn, seaborn, torch
from biolearn.data_library import GeoData
from biolearn.model_gallery import ModelGallery
print("Python/Biolearn check OK")
PY

ENTRYPOINT ["/bin/bash"]
