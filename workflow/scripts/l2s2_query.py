"""L2S2 connectivity query for drug repositioning (B8).

Extracts up/down gene signatures from a DESeq2 results CSV (ranked by padj,
filtered by sign of log2FoldChange) and queries the L2S2 paired-enrichment
GraphQL API. On any transient error (timeout / 5xx) retries with exponential
backoff. On permanent failure writes an empty TSV with the contracted header
and records the reason in the log — never raises.

API: https://l2s2.maayanlab.cloud (GraphQL, public, no auth).
Docs: https://github.com/MaayanLab/L2S2
"""
import csv
import logging
import os
import random
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests


OUTPUT_COLUMNS = ["rank", "perturbagen", "moa", "target", "score", "pvalue", "fdr"]

DEFAULT_ENDPOINT = "https://l2s2.maayanlab.cloud/graphql"

PAIRED_ENRICH_QUERY = """
query PairEnrichmentQuery(
  $genesUp: [String]!
  $genesDown: [String]!
  $first: Int = 250
  $topN: Int = 10000
  $pvalueLe: Float = 0.05
) {
  currentBackground {
    pairedEnrich(
      genesUp: $genesUp
      genesDown: $genesDown
      first: $first
      topN: $topN
      pvalueLe: $pvalueLe
    ) {
      consensusCount
      consensus {
        drug
        oddsRatio
        pvalue
        adjPvalue
        approved
        countSignificant
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
    required = {"gene_name", "log2FoldChange", "padj"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"DE CSV missing required columns: {sorted(missing)}")

    df = df.dropna(subset=["gene_name", "log2FoldChange", "padj"])
    df = df[df["gene_name"].astype(str).str.len() > 0]
    df = df.sort_values("padj", ascending=True, kind="mergesort")

    up = (
        df[df["log2FoldChange"] > 0]["gene_name"]
        .astype(str)
        .drop_duplicates()
        .head(top_up)
        .tolist()
    )
    down = (
        df[df["log2FoldChange"] < 0]["gene_name"]
        .astype(str)
        .drop_duplicates()
        .head(top_down)
        .tolist()
    )
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
        drug = entry.get("drug") or ""
        rows.append(
            {
                "rank": i,
                "perturbagen": drug,
                "moa": "",
                "target": "",
                "score": entry.get("oddsRatio"),
                "pvalue": entry.get("pvalue"),
                "fdr": entry.get("adjPvalue"),
            }
        )
    logger.info("parsed %d consensus perturbagens", len(rows))
    return rows


def write_tsv(out_path: Path, rows: list[dict[str, Any]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    sm = snakemake  # type: ignore[name-defined]
    log_path = sm.log[0] if sm.log else "logs/l2s2_query.log"
    logger = _setup_logger(log_path)

    de_csv = sm.input["de"]
    out_path = Path(sm.output["hits"])
    params = sm.params
    top_up = int(params["top_up"])
    top_down = int(params["top_down"])

    endpoint = os.environ.get("L2S2_API_URL", DEFAULT_ENDPOINT)
    logger.info("L2S2 endpoint: %s", endpoint)

    try:
        genes_up, genes_down = build_signature(de_csv, top_up, top_down, logger)
    except Exception as exc:
        logger.error("failed to build signature: %s", exc)
        _write_empty(out_path)
        return

    if not genes_up and not genes_down:
        logger.error("signature is empty (no genes with padj/log2FC); writing empty TSV")
        _write_empty(out_path)
        return

    payload = {
        "query": PAIRED_ENRICH_QUERY,
        "variables": {
            "genesUp": genes_up,
            "genesDown": genes_down,
            "first": max(top_up, top_down),
            "topN": 10000,
            "pvalueLe": 0.05,
        },
    }

    try:
        response = post_with_retry(endpoint, payload, logger)
        rows = parse_results(response, logger)
    except Exception as exc:
        logger.error("L2S2 query failed, writing empty TSV: %s", exc)
        _write_empty(out_path)
        return

    write_tsv(out_path, rows)
    logger.info("wrote %d rows to %s", len(rows), out_path)


if __name__ == "__main__" or "snakemake" in globals():
    main()
