#!/bin/bash
# Sub-command dispatch for the unified bulk-rnaseq image.
#
#   init      → runs workflow/scripts/init_project.py interactively to write
#               config/config.yaml, config/samples.tsv, config/contrasts.tsv
#               into /project based on host-mounted counts + sample info.
#   jupyter   → launches JupyterLab on :8888 rooted at /project, so
#               /project/notebooks/explore.ipynb can load pipeline outputs.
#   *         → default: runs the Snakemake pipeline. Positional args are
#               forwarded to snakemake (e.g. --cores all, --configfile ...,
#               --dry-run).
set -euo pipefail

case "${1:-}" in
    init)
        shift
        cd /project
        exec micromamba run -n base \
            python /app/workflow/scripts/init_project.py "$@"
        ;;
    jupyter)
        shift
        cd /project
        exec micromamba run -n base jupyter lab \
            --ip=0.0.0.0 --port=8888 --no-browser \
            --ServerApp.root_dir=/project \
            --ServerApp.allow_origin='*' \
            "$@"
        ;;
    *)
        exec micromamba run -n base \
            snakemake --snakefile /app/workflow/Snakefile --use-conda "$@"
        ;;
esac
