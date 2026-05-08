"""
L2S2 connectivity query

Extracts up/down gene signatures from a DESeq2 results CSV (top-N selected
by |wald stat| within each LFC sign) and queries the L2S2 paired-enrichment
GraphQL API. Per-drug MoA is fetched in a follow-up batched fdaCount lookup
(L2S2 returns drug names only in pairedEnrich; MoA lives on the FdaCount
type and the web UI joins them client-side). On any transient error
(timeout / 5xx) retries with exponential backoff. On permanent failure
writes an empty TSV with the contracted header and records the reason in
the log.

API: https://l2s2.maayanlab.cloud (GraphQL, public, no auth).
Docs: https://github.com/MaayanLab/L2S2
"""
import csv
import json
import logging
import random
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests


OUTPUT_COLUMNS = ["rank", "perturbagen", "moa", "score", "pvalue", "fdr"]

ENDPOINT = "https://l2s2.maayanlab.cloud/graphql"

# L2S2 GraphQL gateway rejects single-request bodies above ~tens of KB
# (HTTP 413). Aliased fdaCount queries grow linearly with N; chunk to stay
# safely under the limit and isolate per-chunk failures.
MOA_CHUNK_SIZE = 200

PAIRED_ENRICH_QUERY = """
query PairEnrichmentQuery($genesUp: [String]!, $genesDown: [String]!) {
  currentBackground {
    pairedEnrich(genesUp: $genesUp, genesDown: $genesDown) {
      consensus {
        drug
        oddsRatio
        pvalue
        adjPvalue
      }
    }
  }
}
""".strip()


def _setup_logger(log_path: str) -> logging.Logger:
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("l2s2_query")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.addHandler(sh)
    return logger


def _write_empty(out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(OUTPUT_COLUMNS)


def build_signature(
    de_csv: str, top_up: int, top_down: int, logger: logging.Logger
) -> tuple[list[str], list[str]]:
    df = pd.read_csv(de_csv)
    required = {"gene_name", "log2FoldChange", "stat"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"DE CSV missing required columns: {sorted(missing)}")

    df = df.dropna(subset=["gene_name", "log2FoldChange", "stat"])
    df = df[df["gene_name"].astype(str).str.len() > 0]

    up_df = df.loc[df["log2FoldChange"] > 0].nlargest(top_up, columns="stat")
    down_df = df.loc[df["log2FoldChange"] < 0].nsmallest(top_down, columns="stat")
    up = up_df["gene_name"].astype(str).drop_duplicates().tolist()
    down = down_df["gene_name"].astype(str).drop_duplicates().tolist()
    logger.info("signature built: up=%d down=%d (from %s)", len(up), len(down), de_csv)
    return up, down


def post_with_retry(
    endpoint: str,
    payload: dict[str, Any],
    logger: logging.Logger,
    max_attempts: int = 5,
    base_delay: float = 2.0,
    timeout: float = 120.0,
) -> dict[str, Any]:
    last_error: str = ""
    for attempt in range(1, max_attempts + 1):
        try:
            resp = requests.post(
                endpoint,
                json=payload,
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=timeout,
            )
        except (requests.Timeout, requests.ConnectionError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            logger.warning("attempt %d/%d network error: %s", attempt, max_attempts, last_error)
        else:
            if resp.status_code == 200:
                try:
                    return resp.json()
                except ValueError as exc:
                    raise RuntimeError(f"non-JSON response: {exc}") from exc
            if 500 <= resp.status_code < 600 or resp.status_code == 429:
                last_error = f"HTTP {resp.status_code}: {resp.text[:200]}"
                logger.warning(
                    "attempt %d/%d transient HTTP error: %s", attempt, max_attempts, last_error
                )
            else:
                raise RuntimeError(f"HTTP {resp.status_code}: {resp.text[:500]}")

        if attempt < max_attempts:
            delay = base_delay * (2 ** (attempt - 1)) + random.uniform(0, 1)
            logger.info("sleeping %.1fs before retry", delay)
            time.sleep(delay)

    raise RuntimeError(f"exhausted {max_attempts} retries; last error: {last_error}")


def parse_results(response: dict[str, Any], logger: logging.Logger) -> list[dict[str, Any]]:
    if response.get("errors"):
        raise RuntimeError(f"GraphQL errors: {response['errors']}")
    data = response.get("data") or {}
    bg = data.get("currentBackground") or {}
    paired = bg.get("pairedEnrich") or {}
    consensus = paired.get("consensus") or []

    consensus_sorted = sorted(
        consensus,
        key=lambda r: (
            r.get("adjPvalue") if r.get("adjPvalue") is not None else float("inf"),
            -(r.get("oddsRatio") or 0.0),
        ),
    )

    rows: list[dict[str, Any]] = []
    for i, entry in enumerate(consensus_sorted, start=1):
        rows.append(
            {
                "rank": i,
                "perturbagen": entry.get("drug") or "",
                "moa": "",
                "score": entry.get("oddsRatio"),
                "pvalue": entry.get("pvalue"),
                "fdr": entry.get("adjPvalue"),
            }
        )
    logger.info("parsed %d consensus perturbagens", len(rows))
    return rows


def fetch_moas(
    drugs: list[str], endpoint: str, logger: logging.Logger
) -> dict[str, str]:
    """Batch-lookup MoA for each drug via aliased fdaCount queries, in chunks.

    L2S2's pairedEnrich returns drug names only; MoA lives on FdaCount and
    the web UI joins them client-side. Per-chunk failures are isolated:
    surviving chunks still populate the returned mapping.
    """
    if not drugs:
        return {}
    out: dict[str, str] = {}
    chunks = [drugs[i:i + MOA_CHUNK_SIZE] for i in range(0, len(drugs), MOA_CHUNK_SIZE)]
    for ci, chunk in enumerate(chunks, start=1):
        parts = [
            f"a{i}: fdaCount(perturbation: {json.dumps(d)}) {{ moa }}"
            for i, d in enumerate(chunk)
        ]
        query = "query MoaLookup { " + " ".join(parts) + " }"
        try:
            resp = post_with_retry(endpoint, {"query": query}, logger)
        except Exception as exc:
            logger.warning(
                "MoA chunk %d/%d failed: %s; continuing with remaining chunks",
                ci, len(chunks), exc,
            )
            continue
        if resp.get("errors"):
            logger.warning("MoA chunk %d/%d GraphQL errors: %s", ci, len(chunks), resp["errors"])
        data = resp.get("data") or {}
        for i, drug in enumerate(chunk):
            node = data.get(f"a{i}") or {}
            moa = node.get("moa")
            if moa:
                out[drug] = moa
    logger.info("MoA lookup: %d/%d drugs annotated (%d chunks)", len(out), len(drugs), len(chunks))
    return out


def write_tsv(out_path: Path, rows: list[dict[str, Any]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


snake = snakemake  # noqa: F821  (injected by snakemake script: directive)

log_path = snake.log[0] if snake.log else "logs/l2s2_query.log"
logger = _setup_logger(log_path)

de_csv = snake.input["de"]
out_path = Path(snake.output["hits"])
top_up = int(snake.params["top_up"])
top_down = int(snake.params["top_down"])

logger.info("L2S2 endpoint: %s", ENDPOINT)

try:
    genes_up, genes_down = build_signature(de_csv, top_up, top_down, logger)
except Exception as exc:
    logger.error("failed to build signature: %s", exc)
    _write_empty(out_path)
    sys.exit(0)

if not genes_up or not genes_down:
    logger.error(
        "signature incomplete (up=%d down=%d); paired enrichment needs both — writing empty TSV",
        len(genes_up), len(genes_down),
    )
    _write_empty(out_path)
    sys.exit(0)

payload = {
    "query": PAIRED_ENRICH_QUERY,
    "variables": {"genesUp": genes_up, "genesDown": genes_down},
}

try:
    response = post_with_retry(ENDPOINT, payload, logger)
    rows = parse_results(response, logger)
except Exception as exc:
    logger.error("L2S2 query failed, writing empty TSV: %s", exc)
    _write_empty(out_path)
    sys.exit(0)

drugs = sorted({r["perturbagen"] for r in rows if r["perturbagen"]})
moa_map = fetch_moas(drugs, ENDPOINT, logger)
for r in rows:
    r["moa"] = moa_map.get(r["perturbagen"], "")

write_tsv(out_path, rows)
logger.info("wrote %d rows to %s", len(rows), out_path)
