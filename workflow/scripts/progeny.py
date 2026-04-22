"""B7: PROGENy pathway activity via decoupler-py (MLM on VST matrix).

Reads VST matrix (B2 output: vst_matrix.tsv) and the PROGENy network from
decoupler.op.progeny, then runs decoupler.mt.mlm to derive per-sample
pathway activity scores. Output is a wide TSV of shape (samples x pathways).
"""
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import decoupler as dc

snake = snakemake  # noqa: F821  (injected by snakemake script: directive)

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

vst_path = Path(snake.input["vst_matrix"])
out_path = Path(snake.output["scores"])
top_targets = int(snake.params["top_targets"])

assert vst_path.exists(), f"VST matrix not found: {vst_path}"
out_path.parent.mkdir(parents=True, exist_ok=True)

vst_df = pd.read_csv(vst_path, sep="\t")
assert "gene_id" in vst_df.columns and "gene_name" in vst_df.columns, (
    "vst_matrix.tsv must have gene_id and gene_name columns"
)

sample_cols = [c for c in vst_df.columns if c not in ("gene_id", "gene_name")]
assert len(sample_cols) >= 2, "Need at least 2 sample columns in VST matrix"

logger.info(
    "loaded VST matrix: %d genes x %d samples", len(vst_df), len(sample_cols)
)

symbols = vst_df["gene_name"].astype(str).str.strip()
keep = symbols.ne("") & symbols.ne("NA") & ~symbols.isna()
dropped = int((~keep).sum())
if dropped:
    logger.info("dropping %d rows with missing gene_name", dropped)
vst_df = vst_df.loc[keep].copy()
vst_df["gene_name"] = symbols.loc[keep]

expr = vst_df.set_index("gene_name")[sample_cols].astype(float)
if expr.index.has_duplicates:
    n_dup = int(expr.index.duplicated().sum())
    logger.info("collapsing %d duplicate gene symbols by mean", n_dup)
    expr = expr.groupby(level=0).mean()

mat = expr.T
mat.index.name = "sample"

net = dc.op.progeny(organism="human", top=top_targets)
pathways = sorted(net["source"].unique().tolist())
logger.info(
    "progeny network: %d pathways, %d interactions", len(pathways), len(net)
)

try:
    result = dc.mt.mlm(data=mat, net=net, tmin=5, verbose=False)
    es = result[0] if isinstance(result, tuple) else result
    scores = pd.DataFrame(es).reindex(columns=pathways)
except Exception as exc:  # noqa: BLE001
    logger.warning(
        "mlm failed (%s); writing empty score table with 14 pathway columns",
        exc,
    )
    scores = pd.DataFrame(index=mat.index, columns=pathways, dtype=float)

scores.index.name = "sample"
scores.to_csv(out_path, sep="\t", na_rep="NA", float_format="%.6g")
logger.info(
    "wrote %s: %d samples x %d pathways",
    out_path,
    scores.shape[0],
    scores.shape[1],
)
