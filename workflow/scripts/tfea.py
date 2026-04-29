"""TF activity inference via decoupler-py (ULM on CollecTRI network).

Per-contrast: builds a 1-condition x genes matrix from the DE `stat` column
and runs decoupler.mt.ulm against the CollecTRI network. Output is a long
table of TF activity scores + p-values, one row per TF.

Mirrors the saezlab vignette (tf_bk.html) and the PROGENy implementation —
single decoupler call on per-contrast statistics rather than per-sample VST.
"""
from pathlib import Path

import numpy as np

import decoupler as dc

from _decoupler_common import (
    bh_adjust,
    load_de_stat_matrix,
    setup_logger,
    unpack_decoupler_result,
)

snake = snakemake  # noqa: F821  (injected by snakemake script: directive)

logger = setup_logger(Path(snake.log[0]), "tfea")

de_path = Path(snake.input["de"])
out_path = Path(snake.output["scores"])
out_path.parent.mkdir(parents=True, exist_ok=True)

# decoupler-py official tutorial recommends DESeq2 Wald stat (t-value) input
# and ULM for CollecTRI.
# https://decoupler.readthedocs.io/en/latest/notebooks/bulk/rna.html#enrichment-analysis
mat = load_de_stat_matrix(de_path, logger)
net = dc.op.collectri(organism="human")
logger.info("collectri network: %d TFs, %d interactions",
            net["source"].nunique(), len(net))

result = dc.mt.ulm(data=mat, net=net)
est_df, pval_df = unpack_decoupler_result(result)

scores = (
    est_df.iloc[0]
    .rename("score")
    .rename_axis("source")
    .reset_index()
)
if pval_df is not None and not pval_df.empty:
    scores["p_value"] = pval_df.iloc[0].reindex(scores["source"]).to_numpy()
else:
    scores["p_value"] = np.nan
scores["padj"] = bh_adjust(scores["p_value"])

# Ordering follows the decoupleR R vignette: padj asc, abs(score) desc, source asc.
scores = (
    scores.assign(_abs=scores["score"].abs())
    .sort_values(
        by=["padj", "_abs", "source"],
        ascending=[True, False, True],
        na_position="last",
    )
    .drop(columns="_abs")
    .reset_index(drop=True)
)

scores[["source", "score", "p_value", "padj"]].to_csv(
    out_path, sep="\t", index=False, na_rep="NA", float_format="%.6g"
)
logger.info("wrote %s: %d TFs", out_path, len(scores))
