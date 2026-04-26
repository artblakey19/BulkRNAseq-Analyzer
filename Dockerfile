# Bulk RNA-seq container image.
#
# Single image for three use cases, selected by the first CMD arg:
#   init     → interactive config generator (workflow/scripts/init_project.py)
#   jupyter  → JupyterLab server on port 8888 for notebooks/explore.ipynb
#   default  → Snakemake pipeline (any extra args forwarded to snakemake)
#
# The complete pipeline runtime is baked into the base env at image build
# time. Snakemake still carries per-rule conda directives for native/advanced
# runs, but this Docker image does not pass --use-conda by default, so runtime
# execution does not resolve or create .snakemake/conda environments.
#
# Build (single arch, local):
#   docker build -t bulk-rnaseq:latest .
#
# Build (multi-arch, requires buildx):
#   docker buildx create --use
#   docker buildx build --platform linux/amd64,linux/arm64 \
#       -t bulk-rnaseq:latest --load .
#
# Run examples (bind-mount project at /project):
#   docker run --rm -it -v "$PWD":/project bulk-rnaseq init
#   docker run --rm    -v "$PWD":/project bulk-rnaseq --configfile config/config.yaml --cores all
#   docker run --rm    -v "$PWD":/project -p 8888:8888 bulk-rnaseq jupyter
#
# The entrypoint (docker/entrypoint.sh) picks up the host UID/GID from
# /project's ownership and switches to it via gosu, and sets HOME=/tmp so
# conda/mamba/jupyter have a writable cache. No `-u` or `-e HOME=/tmp` needed.

# Pin the builder/runtime base instead of using :latest. Bump intentionally
# after rebuilding and smoke-testing the image.
FROM mambaorg/micromamba:1.5.10

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# System packages that Snakemake, Quarto, and the UID shim expect.
# tzdata is included so glibc honours `docker run -e TZ=<zone>` at runtime;
# no default zone is baked in.
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        git \
        ca-certificates \
        gosu \
        tzdata \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Base env contains the full Docker runtime: pipeline driver, Python rule
# dependencies, Quarto report stack, JupyterLab, and R/Bioconductor packages.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        snakemake=9 \
        pandas=2 \
        numpy \
        anndata \
        decoupler-py \
        requests \
        pyyaml=6 \
        jupyterlab \
        quarto \
        texlive-core \
        r-base=4.3 \
        r-irkernel \
        r-rmarkdown \
        r-knitr \
        r-quarto \
        r-dplyr \
        r-tidyr \
        r-tibble \
        r-ggplot2 \
        r-scales \
        r-dt \
        r-plotly \
        r-htmltools \
        r-matrixstats \
        r-yaml \
        r-readr \
        r-msigdbr \
        r-data.table \
        bioconductor-deseq2 \
        bioconductor-apeglm \
        bioconductor-summarizedexperiment \
        bioconductor-fgsea \
        bioconductor-clusterprofiler \
        bioconductor-org.hs.eg.db \
 && micromamba clean --all --yes

# Bake the pipeline source into the image. Config and data stay on the host
# and are provided via bind-mount at /project. Notebooks live on the host
# side too (mounted as /project/notebooks/) so users can edit without a
# rebuild.
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER VERSION ./VERSION
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow ./workflow
COPY --chown=$MAMBA_USER:$MAMBA_USER report   ./report
COPY --chown=$MAMBA_USER:$MAMBA_USER config/config.template.yaml ./config/config.template.yaml

# Install entrypoint + sub-command dispatch. entrypoint.sh runs as root (the
# final USER below) so gosu can drop privileges to the host user derived
# from /project ownership. bulk-rnaseq-entry.sh then branches on the first
# arg: init / jupyter / default → snakemake.
USER root
COPY docker/entrypoint.sh          /usr/local/bin/entrypoint.sh
COPY docker/bulk-rnaseq-entry.sh   /usr/local/bin/bulk-rnaseq-entry.sh
RUN chmod +x /usr/local/bin/entrypoint.sh /usr/local/bin/bulk-rnaseq-entry.sh

# Expected runtime layout (bind-mount your project root here):
#   /project/config/{config.yaml, samples.tsv, contrasts.tsv}
#   /project/<counts>.tsv           (path comes from config.yaml)
#   /project/multiqc_data/
#   /project/notebooks/             (for the `jupyter` sub-command)
WORKDIR /project

EXPOSE 8888

ENTRYPOINT ["/usr/local/bin/entrypoint.sh", "/usr/local/bin/bulk-rnaseq-entry.sh"]
CMD ["--cores", "all"]
