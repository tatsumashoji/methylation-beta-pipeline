FROM condaforge/miniforge3:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Tokyo
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

WORKDIR /opt/pipeline

# Minimal OS packages.
# Conda/R source builds mainly use conda-provided compilers,
# but these system packages are still useful for downloads and basic builds.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    git \
    bash \
    build-essential \
    gfortran \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy dependency definitions
COPY environment.yml /opt/pipeline/environment.yml
COPY conda-lock.yml /opt/pipeline/conda-lock.yml

# Create locked conda environment
RUN conda install -y -n base -c conda-forge conda-lock && \
    conda-lock install --name methylation-beta-pipeline /opt/pipeline/conda-lock.yml && \
    conda clean -afy

# Use the locked environment by default
ENV CONDA_PREFIX=/opt/conda/envs/methylation-beta-pipeline
ENV PATH="/opt/conda/bin:/opt/conda/envs/methylation-beta-pipeline/bin:${PATH}"

# Check that conda compiler toolchain exists.
# These tools must be included in environment.yml / conda-lock.yml.
RUN set -eux; \
    test -x "${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-cc"; \
    test -x "${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-c++"; \
    test -x "${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gfortran"; \
    "${CONDA_PREFIX}/bin/R" CMD config CC; \
    "${CONDA_PREFIX}/bin/R" CMD config CXX; \
    "${CONDA_PREFIX}/bin/R" CMD config FC

# Reinstall preprocessCore from source with threading disabled.
# This avoids:
# ERROR; return code from pthread_create() is 22
# in preprocessCore::normalize.quantiles.use.target().
RUN set -eux; \
    export PATH="${CONDA_PREFIX}/bin:/opt/conda/bin:${PATH}"; \
    export CC="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-cc"; \
    export CXX="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-c++"; \
    export FC="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gfortran"; \
    export F77="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-gfortran"; \
    export CPPFLAGS="-I${CONDA_PREFIX}/include"; \
    export CFLAGS="-O2 -pipe"; \
    export CXXFLAGS="-O2 -pipe"; \
    export FFLAGS="-O2 -pipe"; \
    export LDFLAGS="-L${CONDA_PREFIX}/lib"; \
    R_BIN="${CONDA_PREFIX}/bin/R"; \
    RSCRIPT_BIN="${CONDA_PREFIX}/bin/Rscript"; \
    R_LIB="$("${RSCRIPT_BIN}" -e 'cat(.libPaths()[1])')"; \
    echo "R library path: ${R_LIB}"; \
    curl -L -o /tmp/preprocessCore_1.72.0.tar.gz \
      "https://bioconductor.org/packages/3.22/bioc/src/contrib/preprocessCore_1.72.0.tar.gz"; \
    rm -rf "${R_LIB}/preprocessCore"; \
    "${R_BIN}" CMD INSTALL \
      --configure-args=--disable-threading \
      --library="${R_LIB}" \
      /tmp/preprocessCore_1.72.0.tar.gz; \
    rm -f /tmp/preprocessCore_1.72.0.tar.gz

# Sanity check: preprocessCore must work without pthread_create error.
RUN /opt/conda/envs/methylation-beta-pipeline/bin/Rscript -e 'suppressPackageStartupMessages(library(preprocessCore)); cat("preprocessCore version:", as.character(packageVersion("preprocessCore")), "\n"); x <- matrix(runif(1000), nrow = 100); target <- sort(rowMeans(x)); y <- preprocessCore::normalize.quantiles.use.target(x, target); cat("preprocessCore no-thread check OK\n")'

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
ENV BFC_CACHE=/work/vendor/sesame_cache
ENV SESAME_CACHE=/work/vendor/sesame_cache
ENV MANIFEST_DIR=/work/manifest

# Runtime thread limits
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
ENV MKL_NUM_THREADS=1
ENV VECLIB_MAXIMUM_THREADS=1
ENV NUMEXPR_NUM_THREADS=1
ENV RCPP_PARALLEL_NUM_THREADS=1

# R package check
RUN /opt/conda/envs/methylation-beta-pipeline/bin/Rscript -e 'suppressPackageStartupMessages({library(sesame); library(sesameData); library(ExperimentHub); library(data.table); library(dplyr); library(RcppTOML); library(this.path); library(preprocessCore)}); cat("R package check OK\n")'

# Python package check
RUN /opt/conda/envs/methylation-beta-pipeline/bin/python - <<'PY'
import pandas
import numpy
import scipy
import sklearn
import seaborn
import torch
from biolearn.data_library import GeoData
from biolearn.model_gallery import ModelGallery
print("Python/Biolearn check OK")
PY

ENTRYPOINT ["/bin/bash"]
