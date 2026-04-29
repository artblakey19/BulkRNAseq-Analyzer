"""PROGENy pathway activity via decoupler-py (MLM on DESeq2 Wald stat).

Per-contrast: builds a 1-condition x genes matrix from the DE `stat` column
and runs decoupler.mt.mlm against the PROGENy network. Output is a long
table of pathway activity scores + p-values, one row per pathway.

Mirrors the saezlab vignette (pw_bk.html) and the TFEA implementation —
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

logger = setup_logger(Path(snake.log[0]), "progeny")

de_path = Path(snake.input["de"])
out_path = Path(snake.output["scores"])
out_path.parent.mkdir(parents=True, exist_ok=True)

# PROGENy vignette recommands DESeq2(stat) as input.
# https://saezlab.github.io/decoupleR/articles/pw_bk.html
mat = load_de_stat_matrix(de_path, logger)
net = dc.op.progeny(organism="human")
logger.info("progeny network: %d pathways, %d interactions",
            net["source"].nunique(), len(net))

result = dc.mt.mlm(data=mat, net=net)
est_df, pval_df = unpack_decoupler_result(result)

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

scores[["pathway", "score", "p_value", "padj"]].to_csv(
    out_path, sep="\t", index=False, na_rep="NA", float_format="%.6g"
)
logger.info("wrote %s: %d pathways", out_path, len(scores))
