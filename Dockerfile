# syntax=docker/dockerfile:1.6
#
# Bulk RNA-seq downstream pipeline — container image.
#
# Multi-arch: linux/amd64 + linux/arm64 (Apple Silicon).
# Per-rule R/Python envs under workflow/envs/*.yaml are created by Snakemake
# on first run (--use-conda); they live in .snakemake/conda/ on the mounted
# project directory and are reused across runs.
#
# Build (single arch, local):
#   docker build -t bulk-rnaseq:latest .
#
# Build (multi-arch, requires buildx):
#   docker buildx create --use
#   docker buildx build --platform linux/amd64,linux/arm64 \
#       -t bulk-rnaseq:latest --load .
#
# Run (mount project at /project, pass snakemake args as CMD):
#   docker run --rm -v "$PWD":/project bulk-rnaseq -c1
#   docker run --rm -v "$PWD":/project bulk-rnaseq --configfile config/config.yaml -n

FROM mambaorg/micromamba:latest

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# System packages that Snakemake and conda env resolution expect.
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        git \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Base env: Snakemake + mamba frontend + Snakefile helpers (pandas, pyyaml).
# Quarto / R / Bioconductor packages are NOT installed here — each per-rule
# conda env under workflow/envs/*.yaml is resolved lazily at runtime.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        snakemake \
        mamba \
        pandas \
        pyyaml \
 && micromamba clean --all --yes

# Bake the pipeline source into the image. Config and data stay on the host
# and are provided via bind-mount at /project.
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow ./workflow
COPY --chown=$MAMBA_USER:$MAMBA_USER report   ./report

# Expected runtime layout (bind-mount your project root here):
#   /project/config/{config.yaml, samples.tsv, contrasts.tsv}
#   /project/<counts>.tsv           (path comes from config.yaml)
#   /project/multiqc_data/
#   /project/.snakemake/conda/      (created automatically; reused on rerun)
WORKDIR /project

ENTRYPOINT ["micromamba", "run", "-n", "base", \
    "snakemake", "--snakefile", "/app/workflow/Snakefile", "--use-conda"]
CMD ["-c1"]
