"""B6: TF activity inference via decoupler-py (ULM on CollecTRI network).

Ported from the prior R implementation (decoupleR) to Python so that all
decoupler-family analyses (this + PROGENy pathway activity) share a single
tool, language, and conda env.

Logic structure preserved from the original R implementation:
  1. log sink setup
  2. param defaults for missing values
  3. method must be "ulm" (only mode supported)
  4. read DE CSV, keep finite stat + non-empty gene_name
  5. collapse duplicate gene_name by mean
  6. empty-output helper: on any missing network / empty ULM result we write
     header-only TSVs and exit 0 (downstream report degrades gracefully)
  7. fetch CollecTRI with error handling
  8. run ULM with error handling
  9. BH-adjusted p (padj_internal) used for ordering and top_n filter
 10. top_n selection by abs(score) desc, ranked, with `rank` column
"""
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import decoupler as dc

snake = snakemake  # noqa: F821  (injected by snakemake script: directive)

# --- Log sink -------------------------------------------------------------
log_path = Path(snake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=log_path,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger("tfea")
sys.stderr = open(log_path, "a")

np.random.seed(42)

# --- Inputs / outputs / params --------------------------------------------
de_path = Path(snake.input["de"])
scores_out = Path(snake.output["scores"])
top_out = Path(snake.output["top"])


def _param(key, default):
    val = snake.params.get(key)
    return default if val is None else val


method = _param("method", "ulm")
min_size = int(_param("min_size", 5))
padj_cutoff = float(_param("padj_cutoff", 0.05))
top_n = int(_param("top_n", 50))
split_complexes = bool(_param("split_complexes", False))

if method != "ulm":
    raise ValueError(
        f"Unsupported TF activity method: {method!r}. Only 'ulm' is supported."
    )

assert de_path.exists(), f"DE CSV not found: {de_path}"
scores_out.parent.mkdir(parents=True, exist_ok=True)
top_out.parent.mkdir(parents=True, exist_ok=True)

SCORE_COLS = ["source", "condition", "score", "p_value", "statistic"]
TOP_COLS = SCORE_COLS + ["rank"]


def write_empty_outputs(reason: str) -> None:
    logger.info(reason)
    pd.DataFrame(columns=SCORE_COLS).to_csv(scores_out, sep="\t", index=False)
    pd.DataFrame(columns=TOP_COLS).to_csv(top_out, sep="\t", index=False)


def bh_adjust(pvals) -> np.ndarray:
    """Benjamini-Hochberg FDR. Preserves NaN positions; monotonic from largest."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    valid = ~np.isnan(p)
    n = int(valid.sum())
    if n == 0:
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    adj = np.minimum(1.0, pv * n / ranks)
    sorted_adj = adj[order]
    sorted_adj = np.minimum.accumulate(sorted_adj[::-1])[::-1]
    adj_out = np.empty(n)
    adj_out[order] = sorted_adj
    out[valid] = adj_out
    return out


# --- Load DE results and build input matrix -------------------------------
de_res = pd.read_csv(de_path)
required_cols = ["gene_name", "stat"]
missing = [c for c in required_cols if c not in de_res.columns]
if missing:
    raise ValueError(f"DE results missing required columns: {missing}")

mask = (
    de_res["gene_name"].notna()
    & (de_res["gene_name"].astype(str).str.strip() != "")
    & de_res["stat"].notna()
    & np.isfinite(de_res["stat"].astype(float))
)
stat_tbl = (
    de_res.loc[mask, ["gene_name", "stat"]]
    .assign(gene_name=lambda d: d["gene_name"].astype(str))
    .groupby("gene_name", as_index=False)["stat"]
    .mean()
)

if stat_tbl.empty:
    write_empty_outputs(
        "No finite gene_name/stat pairs available for TF activity inference."
    )
    sys.exit(0)

# decoupler expects samples x features. One "condition" pseudo-sample.
mat = (
    stat_tbl.set_index("gene_name")["stat"]
    .astype(float)
    .to_frame()
    .T
)
mat.index = ["condition"]
mat.index.name = "sample"
logger.info("Built decoupler input: %d genes x 1 condition", mat.shape[1])

# --- Fetch CollecTRI ------------------------------------------------------
# Newer decoupler-py (>=2.x) dropped the `split_complexes` keyword from
# dc.op.collectri(); older versions still accept it. Try the legacy signature
# first, fall back to the minimal one on TypeError.
try:
    try:
        network = dc.op.collectri(organism="human",
                                  split_complexes=split_complexes)
    except TypeError:
        if split_complexes:
            logger.info(
                "dc.op.collectri() does not accept split_complexes in this "
                "decoupler version; falling back to default behavior."
            )
        network = dc.op.collectri(organism="human")
except Exception as exc:  # noqa: BLE001
    write_empty_outputs(f"CollecTRI network fetch failed: {exc}")
    sys.exit(0)

if network is None or len(network) == 0:
    write_empty_outputs("CollecTRI network unavailable; wrote empty TF activity outputs.")
    sys.exit(0)

# decoupler-py <2 exposed the edge sign as `mor`; >=2 renamed it to `weight`.
# Accept either so the script works across both generations.
weight_col = next((c for c in ("weight", "mor") if c in network.columns), None)
for col in ("source", "target"):
    if col not in network.columns:
        write_empty_outputs(
            f"CollecTRI network missing required column '{col}'; wrote empty outputs."
        )
        sys.exit(0)
if weight_col is None:
    write_empty_outputs(
        "CollecTRI network missing edge-weight column ('weight' or 'mor'); wrote empty outputs."
    )
    sys.exit(0)
if weight_col != "weight":
    network = network.rename(columns={weight_col: "weight"})

logger.info(
    "CollecTRI network: %d TFs, %d interactions",
    network["source"].nunique(),
    len(network),
)

# --- Run ULM ---------------------------------------------------------------
try:
    result = dc.mt.ulm(data=mat, net=network, tmin=min_size, verbose=False)
except Exception as exc:  # noqa: BLE001
    write_empty_outputs(f"decoupler.mt.ulm failed: {exc}")
    sys.exit(0)

if isinstance(result, tuple):
    est_df = result[0]
    pval_df = result[1] if len(result) > 1 else None
else:
    est_df = result
    pval_df = None

if est_df is None or est_df.empty:
    write_empty_outputs("No TF activities returned by decoupler.mt.ulm.")
    sys.exit(0)

# Reshape wide (1 condition x n_TFs) → long (n_TFs rows)
scores = (
    est_df.iloc[0]
    .rename("score")
    .rename_axis("source")
    .reset_index()
)
scores["condition"] = "condition"
if pval_df is not None and not pval_df.empty:
    scores["p_value"] = pval_df.iloc[0].reindex(scores["source"]).to_numpy()
else:
    scores["p_value"] = np.nan
scores["statistic"] = method

scores["padj_internal"] = bh_adjust(scores["p_value"])

# Ordering matches the R implementation: padj asc, abs(score) desc, source asc.
scores = (
    scores.assign(_abs=scores["score"].abs())
    .sort_values(
        by=["padj_internal", "_abs", "source"],
        ascending=[True, False, True],
        na_position="last",
    )
    .drop(columns="_abs")
    .reset_index(drop=True)
)

# Top filter + rank (abs(score) desc, then padj asc, then source).
top = scores.loc[
    scores["padj_internal"].notna() & (scores["padj_internal"] < padj_cutoff)
].copy()
top = (
    top.assign(_abs=top["score"].abs())
    .sort_values(
        by=["_abs", "padj_internal", "source"],
        ascending=[False, True, True],
    )
    .drop(columns="_abs")
    .head(top_n)
    .reset_index(drop=True)
)
top["rank"] = np.arange(1, len(top) + 1, dtype=int)

scores[SCORE_COLS].to_csv(scores_out, sep="\t", index=False, na_rep="NA")
top[TOP_COLS].to_csv(top_out, sep="\t", index=False, na_rep="NA")

logger.info(
    "TF activity analysis complete. Scores: %d; top: %d",
    len(scores),
    len(top),
)
