"""PROGENy pathway activity via decoupler-py (MLM on DESeq2 Wald stat).

Per-contrast: builds a 1-condition x genes matrix from the DE `stat` column
and runs decoupler.mt.mlm against the PROGENy network. Output is a long
table of pathway activity scores + p-values, one row per pathway.

Mirrors the saezlab vignette (pw_bk.html) and the TFEA implementation —
single decoupler call on per-contrast statistics rather than per-sample VST.
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
logger = logging.getLogger("progeny")
sys.stderr = open(log_path, "a")

np.random.seed(42)

# --- Inputs / outputs / params --------------------------------------------
de_path = Path(snake.input["de"])
out_path = Path(snake.output["scores"])
top_targets = int(snake.params["top_targets"])
out_path.parent.mkdir(parents=True, exist_ok=True)

SCORE_COLS = ["pathway", "score", "p_value", "padj"]


def write_empty_outputs(reason: str) -> None:
    logger.info(reason)
    pd.DataFrame(columns=SCORE_COLS).to_csv(out_path, sep="\t", index=False)


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
    write_empty_outputs("No finite gene_name/stat pairs available for PROGENy.")
    sys.exit(0)

mat = (
    stat_tbl.set_index("gene_name")["stat"]
    .astype(float)
    .to_frame()
    .T
)
mat.index = ["condition"]
mat.index.name = "sample"
logger.info("Built decoupler input: %d genes x 1 condition", mat.shape[1])

# --- Fetch PROGENy network ------------------------------------------------
try:
    net = dc.op.progeny(organism="human", top=top_targets)
except Exception as exc:  # noqa: BLE001
    write_empty_outputs(f"PROGENy network fetch failed: {exc}")
    sys.exit(0)

if net is None or len(net) == 0:
    write_empty_outputs("PROGENy network unavailable.")
    sys.exit(0)

logger.info("progeny network: %d pathways, %d interactions",
            net["source"].nunique(), len(net))

# --- Run MLM --------------------------------------------------------------
try:
    result = dc.mt.mlm(data=mat, net=net, tmin=5, verbose=False)
except Exception as exc:  # noqa: BLE001
    write_empty_outputs(f"decoupler.mt.mlm failed: {exc}")
    sys.exit(0)

# decoupler-py <2 returns (estimate, pvalue) tuple; >=2 returns DataFrame
if isinstance(result, tuple):
    est_df = result[0]
    pval_df = result[1] if len(result) > 1 else None
else:
    est_df = result
    pval_df = None

if est_df is None or est_df.empty:
    write_empty_outputs("No pathway activities returned by decoupler.mt.mlm.")
    sys.exit(0)

# Reshape wide (1 condition x n_pathways) → long (n_pathways rows)
scores = (
    est_df.iloc[0]
    .rename("score")
    .rename_axis("pathway")
    .reset_index()
)
if pval_df is not None and not pval_df.empty:
    scores["p_value"] = pval_df.iloc[0].reindex(scores["pathway"]).to_numpy()
else:
    scores["p_value"] = np.nan

scores["padj"] = bh_adjust(scores["p_value"])

scores = (
    scores.assign(_abs=scores["score"].abs())
    .sort_values(by=["_abs", "pathway"], ascending=[False, True])
    .drop(columns="_abs")
    .reset_index(drop=True)
)

scores[SCORE_COLS].to_csv(out_path, sep="\t", index=False, na_rep="NA",
                          float_format="%.6g")
logger.info("wrote %s: %d pathways", out_path, len(scores))
