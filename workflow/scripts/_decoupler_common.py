"""Shared helpers for decoupler-py per-contrast scripts (tfea.py, progeny.py).

Both scripts:
- write logs through the same per-contrast log sink
- collapse a DESeq2 results CSV into a 1xN "condition" matrix on `stat`
- BH-adjust p-values
- unpack decoupler.mt.* outputs that may be either a (estimate, pvalue) tuple
  (decoupler-py <2 and current 2.x) or a single DataFrame (some 2.x builds)

Centralizing here so the two scripts cannot drift in subtle ways
(filtering rules, NaN handling, log formatting, decoupler version handling).
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def setup_logger(log_path: Path, name: str) -> logging.Logger:
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=str(log_path),
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    sys.stderr = open(log_path, "a")
    return logging.getLogger(name)


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


def load_de_stat_matrix(
    de_csv_path: Path, logger: logging.Logger
) -> pd.DataFrame:
    """Read DESeq2 results CSV into a 1xN samples-by-genes matrix on `stat`.

    Drops rows with empty/NaN gene_name or non-finite stat, collapses duplicate
    gene_name by mean, returns a one-row DataFrame indexed `["condition"]`
    suitable as `data` for `decoupler.mt.*`.
    """
    de_res = pd.read_csv(de_csv_path)
    required = ["gene_name", "stat"]
    missing = [c for c in required if c not in de_res.columns]
    if missing:
        raise ValueError(f"DE results missing required columns: {missing}")

    raw_gene = de_res["gene_name"]
    gene_str = raw_gene.astype(str).str.strip()
    stat = pd.to_numeric(de_res["stat"], errors="coerce")
    mask = raw_gene.notna() & gene_str.ne("") & np.isfinite(stat)

    stat_tbl = (
        pd.DataFrame({"gene_name": gene_str[mask], "stat": stat[mask]})
        .groupby("gene_name", as_index=False)["stat"]
        .mean()
    )

    mat = (
        stat_tbl.set_index("gene_name")["stat"]
        .astype(float)
        .to_frame()
        .T
    )
    mat.index = ["condition"]
    mat.index.name = "sample"
    logger.info("Built decoupler input: %d genes x 1 condition", mat.shape[1])
    return mat


def unpack_decoupler_result(result):
    """decoupler.mt.* may return (estimate, pvalue) tuple or a single DataFrame."""
    if isinstance(result, tuple):
        est_df = result[0]
        pval_df = result[1] if len(result) > 1 else None
    else:
        est_df = result
        pval_df = None
    return est_df, pval_df
