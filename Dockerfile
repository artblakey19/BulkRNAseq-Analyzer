# Bulk RNA-seq container image.
#
# Single image for three use cases, selected by the first CMD arg:
#   init     → interactive config generator (workflow/scripts/init_project.py)
#   jupyter  → JupyterLab server on port 8888 for notebooks/explore.ipynb
#   default  → Snakemake pipeline (any extra args forwarded to snakemake)
#
# The R/Bioconductor stack that notebooks read is baked into the base env
# so the Jupyter use case is zero-setup. Pipeline rules still resolve their
# own conda envs under workflow/envs/*.yaml via --use-conda on first run
# (cached in .snakemake/conda/ on the bind-mounted project directory).
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

FROM mambaorg/micromamba:latest

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# System packages that Snakemake, conda env resolution, and the UID shim
# expect.
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        git \
        ca-certificates \
        gosu \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Base env combines:
#   - pipeline drivers (snakemake, mamba, pandas, pyyaml)
#   - JupyterLab + R kernel + the R/Bioconductor stack that explore.ipynb
#     loads (versions aligned with workflow/envs/r-*.yaml so .rds files
#     written by pipeline rules are readable by the notebook)
# Python per-rule deps (decoupler-py, requests) stay out of base — those
# rules still resolve their own envs under workflow/envs/py-*.yaml.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        snakemake \
        conda \
        mamba \
        pandas \
        pyyaml \
        jupyterlab \
        r-base=4.3 \
        r-irkernel \
        r-tidyverse \
        r-matrixstats \
        r-pheatmap \
        r-ggrepel \
        r-yaml \
        r-readr \
        r-msigdbr \
        bioconductor-deseq2 \
        bioconductor-apeglm \
        bioconductor-summarizedexperiment \
        bioconductor-fgsea \
        bioconductor-clusterprofiler \
        bioconductor-enrichplot \
        bioconductor-org.hs.eg.db \
        bioconductor-complexheatmap \
 && micromamba clean --all --yes

# Bake the pipeline source into the image. Config and data stay on the host
# and are provided via bind-mount at /project. Notebooks live on the host
# side too (mounted as /project/notebooks/) so users can edit without a
# rebuild.
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow ./workflow
COPY --chown=$MAMBA_USER:$MAMBA_USER report   ./report

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
#   /project/.snakemake/conda/      (created automatically; reused on rerun)
WORKDIR /project

EXPOSE 8888

ENTRYPOINT ["/usr/local/bin/entrypoint.sh", "/usr/local/bin/bulk-rnaseq-entry.sh"]
CMD ["--cores", "all"]
