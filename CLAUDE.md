# Project file map — what to touch together

Internal guide for making coherent changes. The pipeline spans Snakefile
rules, conda env yamls, R/Python scripts, Quarto report partials, two
explore notebooks, and two Colab notebooks. A "simple" change to one layer
almost always implies mirrored edits in the others; this file catalogs
those pairings.

---

## Layer map

| Concern | Where it lives |
| --- | --- |
| DAG / rule definitions | `workflow/Snakefile` + `workflow/rules/*.smk` |
| Per-rule conda envs | `workflow/envs/*.yaml` |
| Rule implementations | `workflow/scripts/*.{R,py}` |
| Config schema & defaults | `config/config.template.yaml` (tracked), `config/config.yaml` (gitignored, user-generated) |
| CI test fixture config | `tests/test_data/config_test.yaml` |
| HTML report assembly | `report/template.qmd` (entry) + `report/sections/_*.qmd` (per-stage partials) |
| Local interactive notebook | `notebooks/explore.ipynb` + `notebooks/explore.ko.ipynb` |
| Colab end-to-end notebook | `notebooks/colab_pipeline.ipynb` + `notebooks/colab_pipeline.ko.ipynb` |
| Docker image definition | `Dockerfile`, `docker/entrypoint.sh`, `docker/bulk-rnaseq-entry.sh` |
| CI | `.github/workflows/docker.yml` (image build), `.github/workflows/test.yml` (snakemake dry-run) |
| User-facing docs | `README.md` + `README.ko.md` |
| Workflow diagram | `docs/rulegraph.svg` (generated from `snakemake --rulegraph`) |

---

## Interdependency recipes

These describe the minimal set of files to touch for common change types.
When editing, grep before committing — `grep -rn "<old-name>" --include='*.smk' --include='*.py' --include='*.R' --include='*.qmd' --include='*.yaml' --include='*.ipynb' --include='*.md'` catches strays.

### Adding a new analysis step (e.g. a new per-contrast tool)

1. **Snakefile**
   - Add rule declaration: `workflow/rules/<area>.smk` (enrichment / qc / de / report)
   - Add output to `rule all` in `workflow/Snakefile`
2. **Conda env** — new `workflow/envs/<x>.yaml` or reuse an existing one (r-enrichment, py-decoupler, etc.)
3. **Script** — `workflow/scripts/<x>.{R,py}`; bind via the rule's `script:` directive
4. **Config**
   - `config/config.template.yaml` — add section with default values
   - `tests/test_data/config_test.yaml` — mirror with values valid for the fixture
   - `workflow/scripts/init_project.py` — only if the step needs user-provided params (otherwise defaults from template suffice)
5. **Report**
   - `report/sections/_<x>.qmd` — new partial
   - `report/template.qmd` — add `{{< include sections/_<x>.qmd >}}`
6. **Notebooks** (four!)
   - `notebooks/explore.ipynb` — add a cell that reads the new output and plots
   - `notebooks/explore.ko.ipynb` — mirror of above, Korean markdown
   - `notebooks/colab_pipeline.ipynb` — add a cell in section 11 (Explore) with `%%R` magic
   - `notebooks/colab_pipeline.ko.ipynb` — Korean mirror
7. **README**
   - `README.md` §"Report sections" table — add row
   - `README.ko.md` §"리포트 섹션" table — add row
8. **Docker image** — Docker's default execution path uses the baked base env, not `--use-conda`. Any package needed by a pipeline rule or notebook must be mirrored in the base-env install in `Dockerfile`.
9. **Rulegraph** — regenerate `docs/rulegraph.svg` with `snakemake --rulegraph | dot -Tsvg > docs/rulegraph.svg` so README's Workflow diagram reflects the new node.
10. **Smoke test** — the dry-run workflow (`.github/workflows/test.yml`) picks up the new rule automatically, but verify `snakemake --dry-run --configfile tests/test_data/config_test.yaml` still resolves locally before pushing.

### Renaming a rule

Rule names appear as directory names, config keys, and labels. Touch:

- `workflow/rules/*.smk` — rule declaration + `output:` directory
- `workflow/Snakefile` — `rule all` output paths
- Any rule whose `input:` references the old path
- `workflow/scripts/*` — rename file if it mirrored the rule name; update internal logger/label strings
- `config/config.template.yaml`, `config/config.yaml`, `tests/test_data/config_test.yaml` — config key
- `report/template.qmd` + `report/sections/_*.qmd` — path strings, config key reads (`cfg$<key>`)
- `notebooks/*.ipynb` (all four) — path strings in `rp("results", "<rule>", ...)` or equivalents
- `docs/rulegraph.svg` — regenerate

### Changing DE / enrichment cutoffs or adding an MSigDB collection / ORA database

- `config/config.template.yaml` is the **source of truth** for the defaults shipped to new users. Update here.
- `tests/test_data/config_test.yaml` — mirror if the schema changed; otherwise only touch if fixture needs new keys.
- The report and notebooks read values from the user's live `config.yaml`, so they don't need edits — they just surface whatever is set.
- If adding a new **MSigDB collection id** (e.g. C8), also check `workflow/scripts/gsea.R` dispatches it correctly.
- ORA databases use MSigDB collection strings (e.g. `C5:GO:BP`, `C2:CP:REACTOME`, `H`) parsed by `workflow/scripts/ora.R` the same way `gsea.R` parses its collection ids. The literal id `KEGG` is the one special case, dispatched through `enrichKEGG`. Adding a new MSigDB-backed db only needs a config entry; adding non-MSigDB sources requires a new branch in `ora.R`.

### Changing an environment's package set

- Edit `workflow/envs/<env>.yaml`. Snakemake hashes this file; on next `--use-conda` run it rebuilds the env at a new hash (old `.snakemake/conda/<old-hash>` lingers until user cleans up — safe to ignore).
- `Dockerfile` base env — mirror package changes in `Dockerfile`'s `micromamba install` line. Docker runs use the baked base env by default, so rule dependencies and notebook dependencies both belong there.
- `notebooks/colab_pipeline*.ipynb` — the Colab Explore section installs `fgsea/msigdbr/httr/jsonlite` into the r-deseq2 env on demand. If a notebook cell starts requiring a new package, extend the `need = [...]` list in that install cell.

### Renaming display terminology (TFEA, cMap, PROGENy, etc.)

Separate from rule identifiers. Display strings live in:

- `report/sections/_*.qmd` — section headings
- `README.md` / `README.ko.md` — report-sections table
- `notebooks/*.ipynb` — markdown cells
- `docs/rulegraph.svg` — regenerate (labels derive from rule names; if display != rule name, update `docs/_rulegraph_tidy.py` LABEL_RENAMES **or** apply manually in the SVG source)

### Adding a new Colab notebook cell (Explore section)

Because there are **four notebooks** (explore EN/KO + colab EN/KO) with overlapping plot logic, adding an analysis means four cells. Current recommendation: keep one as canonical (`notebooks/explore.ipynb`), copy to the others, translate markdown for `.ko` variants, adjust for Colab's `%%R` magic syntax in `colab_pipeline.*`. There is **no automated sync** — manual discipline.

---

## Things that don't need mirrored updates (intentional)

- **`.snakemake/conda/<hash>/`** — purely cache. Never commit; regenerated on demand.
- **`results/`** — gitignored pipeline outputs.
- **`config/config.yaml`, `config/samples.tsv`, `config/contrasts.tsv`** — user-generated via `init_project.py`, gitignored. Only the `.template` versions are tracked.
- **Old Docker image tags** (`bulk-rnaseq-jupyter`) — the image was merged into the main `bulk-rnaseq` tag. Any leftover references in external docs or local copies can be ignored; this repo has none.

---

## Red flags when reviewing a PR

- Change touches `workflow/rules/` but not `config/config.template.yaml` → likely missing a new config key
- Change renames a rule but `docs/rulegraph.svg` is untouched → stale diagram
- Change adds an output path but no `notebooks/` cell reads it → the local Jupyter + Colab Explore flows will silently miss this plot
- Change edits only `README.md` (no `README.ko.md`) → the two drift; they should move together
- Change edits only one of `notebooks/explore.ipynb` / `notebooks/colab_pipeline.ipynb` (or their `.ko` siblings) → the two Jupyter surfaces drift
