#!/usr/bin/env python3
"""Fail-fast validation for the Bulk RNA-seq workflow config and metadata."""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")

SUPPORTED_CONFIG_SCHEMA: dict[str, Any] = {
    "project": {
        "name": None,
        "analyst": None,
        "experiment": {
            "cell_line": None,
            "compound": None,
            "dose": None,
            "duration": None,
            "notes": None,
        },
    },
    "input": {
        "counts_tsv": None,
        "multiqc_data_dir": None,
        "samples_tsv": None,
        "contrasts_tsv": None,
    },
    "de": {
        "prefilter_min_count": None,
        "primary": {"padj": None, "abs_lfc": None},
        "secondary": {"padj": None, "abs_lfc": None},
    },
    "enrichment": {
        "enabled": None,
        "gsea": {
            "ranking": None,
            "min_size": None,
            "max_size": None,
            "seed": None,
            "collections": [{"id": None, "label": None}],
        },
        "ora": {
            "min_input_genes": None,
            "max_input_genes": None,
            "databases": [{"id": None, "label": None}],
        },
    },
    "tfea": {
        "enabled": None,
        "split_complexes": None,
        "min_size": None,
        "padj_cutoff": None,
        "top_n": None,
    },
    "progeny": {
        "enabled": None,
        "top_targets": None,
    },
    "cmap": {
        "enabled": None,
        "top_up": None,
        "top_down": None,
    },
    "report": {
        "formats": None,
    },
}

REQUIRED_SAMPLE_COLUMNS = ["sample", "condition", "replicate", "batch"]
REQUIRED_CONTRAST_COLUMNS = [
    "contrast_id",
    "factor",
    "numerator",
    "denominator",
    "description",
]


def _format_path(path: tuple[str, ...]) -> str:
    return ".".join(path) if path else "<root>"


def _check_unknown_keys(obj: Any, schema: Any, path: tuple[str, ...], errors: list[str]) -> None:
    if isinstance(schema, dict):
        if not isinstance(obj, dict):
            errors.append(f"{_format_path(path)} must be a mapping.")
            return
        for key in obj:
            if key not in schema:
                errors.append(f"Unsupported config key: {_format_path(path + (str(key),))}")
        for key, child_schema in schema.items():
            if key in obj and child_schema is not None:
                _check_unknown_keys(obj[key], child_schema, path + (key,), errors)
        return

    if isinstance(schema, list):
        if not isinstance(obj, list):
            errors.append(f"{_format_path(path)} must be a list.")
            return
        child_schema = schema[0]
        for idx, item in enumerate(obj):
            _check_unknown_keys(item, child_schema, path + (f"[{idx}]",), errors)


def _require_sections(config: dict[str, Any], sections: list[str], errors: list[str]) -> None:
    for section in sections:
        if section not in config:
            errors.append(f"Missing required config section: {section}")


def _require_keys(config: dict[str, Any], section: str, keys: list[str], errors: list[str]) -> None:
    value = config.get(section)
    if not isinstance(value, dict):
        return
    for key in keys:
        if key not in value:
            errors.append(f"Missing required config key: {section}.{key}")


def _require_nested_keys(
    config: dict[str, Any],
    path: tuple[str, ...],
    keys: list[str],
    errors: list[str],
) -> None:
    value: Any = config
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return
        value = value[key]
    if not isinstance(value, dict):
        return
    for key in keys:
        if key not in value:
            errors.append(f"Missing required config key: {_format_path(path + (key,))}")


def _require_bool(config: dict[str, Any], section: str, key: str, errors: list[str]) -> None:
    value = config.get(section, {}).get(key)
    if not isinstance(value, bool):
        errors.append(f"{section}.{key} must be true or false.")


def _nonempty_str(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _check_id_series(values: list[Any], label: str, errors: list[str]) -> None:
    bad = [
        str(v)
        for v in values
        if pd.isna(v) or not _nonempty_str(str(v)) or not ID_RE.fullmatch(str(v))
    ]
    if bad:
        errors.append(
            f"{label} values must match {ID_RE.pattern}; invalid values: {', '.join(bad)}"
        )


def _check_duplicates(values: list[Any], label: str, errors: list[str]) -> None:
    seen: set[str] = set()
    dupes: list[str] = []
    for value in map(str, values):
        if value in seen and value not in dupes:
            dupes.append(value)
        seen.add(value)
    if dupes:
        errors.append(f"Duplicate {label}: {', '.join(dupes)}")


def _read_counts_header(path: str | Path) -> list[str]:
    counts_path = Path(path)
    if not counts_path.exists():
        raise ValueError(f"input.counts_tsv does not exist: {counts_path}")
    with counts_path.open() as handle:
        return handle.readline().rstrip("\n").split("\t")


def load_config_tables(config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Load metadata files referenced by config for validation and Snakefile use."""
    input_cfg = config.get("input", {})
    samples_path = Path(input_cfg.get("samples_tsv", ""))
    contrasts_path = Path(input_cfg.get("contrasts_tsv", ""))
    if not samples_path.exists():
        raise ValueError(f"input.samples_tsv does not exist: {samples_path}")
    if not contrasts_path.exists():
        raise ValueError(f"input.contrasts_tsv does not exist: {contrasts_path}")

    samples = pd.read_csv(samples_path, sep="\t")
    contrasts = pd.read_csv(contrasts_path, sep="\t")
    counts_header = _read_counts_header(input_cfg.get("counts_tsv", ""))
    return samples, contrasts, counts_header


def validate_config(
    config: dict[str, Any],
    samples: pd.DataFrame,
    contrasts: pd.DataFrame,
    counts_header: list[str],
) -> None:
    """Raise ValueError if config/metadata would lead to ambiguous workflow behavior."""
    errors: list[str] = []

    _require_sections(
        config,
        ["project", "input", "de", "enrichment", "tfea", "progeny", "cmap", "report"],
        errors,
    )
    _check_unknown_keys(config, SUPPORTED_CONFIG_SCHEMA, tuple(), errors)

    _require_keys(config, "input", ["counts_tsv", "samples_tsv", "contrasts_tsv"], errors)
    _require_keys(
        config,
        "de",
        ["prefilter_min_count", "primary", "secondary"],
        errors,
    )
    _require_nested_keys(config, ("de", "primary"), ["padj", "abs_lfc"], errors)
    _require_nested_keys(config, ("de", "secondary"), ["padj", "abs_lfc"], errors)
    _require_nested_keys(config, ("enrichment",), ["enabled", "gsea", "ora"], errors)
    _require_nested_keys(
        config,
        ("enrichment", "gsea"),
        ["ranking", "min_size", "max_size", "seed", "collections"],
        errors,
    )
    _require_nested_keys(
        config,
        ("enrichment", "ora"),
        ["min_input_genes", "max_input_genes", "databases"],
        errors,
    )
    _require_keys(
        config,
        "tfea",
        ["enabled", "split_complexes", "min_size", "padj_cutoff", "top_n"],
        errors,
    )
    _require_keys(config, "progeny", ["enabled", "top_targets"], errors)
    _require_keys(config, "cmap", ["enabled", "top_up", "top_down"], errors)
    _require_keys(config, "report", ["formats"], errors)
    for section in ["enrichment", "tfea", "progeny", "cmap"]:
        _require_bool(config, section, "enabled", errors)

    report_formats = config.get("report", {}).get("formats")
    if not isinstance(report_formats, list) or not report_formats:
        errors.append("report.formats must be a non-empty list.")
    else:
        bad_formats = [str(fmt) for fmt in report_formats if fmt not in {"html", "pdf"}]
        if bad_formats:
            errors.append(f"Unsupported report.formats values: {', '.join(bad_formats)}")

    if config.get("enrichment", {}).get("gsea", {}).get("ranking") not in {"stat", "signed_p"}:
        errors.append("enrichment.gsea.ranking must be 'stat' or 'signed_p'.")
    missing_sample_cols = [c for c in REQUIRED_SAMPLE_COLUMNS if c not in samples.columns]
    extra_sample_cols = [c for c in samples.columns if c not in REQUIRED_SAMPLE_COLUMNS]
    if missing_sample_cols:
        errors.append(f"samples.tsv missing required columns: {', '.join(missing_sample_cols)}")
    if extra_sample_cols:
        errors.append(f"samples.tsv has unsupported columns: {', '.join(extra_sample_cols)}")

    missing_contrast_cols = [c for c in REQUIRED_CONTRAST_COLUMNS if c not in contrasts.columns]
    extra_contrast_cols = [c for c in contrasts.columns if c not in REQUIRED_CONTRAST_COLUMNS]
    if missing_contrast_cols:
        errors.append(f"contrasts.tsv missing required columns: {', '.join(missing_contrast_cols)}")
    if extra_contrast_cols:
        errors.append(f"contrasts.tsv has unsupported columns: {', '.join(extra_contrast_cols)}")

    if not errors and len(counts_header) < 3:
        errors.append("counts TSV must contain gene_id, gene_name, and at least one sample column.")
    if not errors and counts_header[:2] != ["gene_id", "gene_name"]:
        errors.append("counts TSV first two columns must be gene_id and gene_name.")

    if "sample" in samples.columns:
        sample_ids = samples["sample"].astype(str).tolist()
        _check_duplicates(sample_ids, "sample ids", errors)
        _check_id_series(sample_ids, "sample", errors)
        if len(counts_header) >= 3:
            count_samples = counts_header[2:]
            missing = sorted(set(sample_ids) - set(count_samples))
            extra = sorted(set(count_samples) - set(sample_ids))
            if missing:
                errors.append(f"samples.tsv samples missing from counts TSV: {', '.join(missing)}")
            if extra:
                errors.append(f"counts TSV sample columns not declared in samples.tsv: {', '.join(extra)}")

    if "condition" in samples.columns:
        conditions = samples["condition"].tolist()
        empty_conditions = [v for v in conditions if pd.isna(v) or not _nonempty_str(str(v))]
        if empty_conditions:
            errors.append("samples.tsv condition values must be non-empty.")

    if "contrast_id" in contrasts.columns:
        contrast_ids = contrasts["contrast_id"].astype(str).tolist()
        _check_duplicates(contrast_ids, "contrast ids", errors)
        _check_id_series(contrast_ids, "contrast_id", errors)

    if set(REQUIRED_CONTRAST_COLUMNS).issubset(contrasts.columns) and "condition" in samples.columns:
        condition_set = set(samples["condition"].astype(str))
        for _, row in contrasts.iterrows():
            cid = str(row["contrast_id"])
            factor = str(row["factor"])
            numerator = str(row["numerator"])
            denominator = str(row["denominator"])
            if factor != "condition":
                errors.append(
                    f"contrast_id '{cid}' uses factor '{factor}', but only factor 'condition' is supported."
                )
            if numerator not in condition_set:
                errors.append(f"contrast_id '{cid}' numerator '{numerator}' not found in samples.condition.")
            if denominator not in condition_set:
                errors.append(
                    f"contrast_id '{cid}' denominator '{denominator}' not found in samples.condition."
                )
            if numerator == denominator:
                errors.append(f"contrast_id '{cid}' numerator and denominator must differ.")

    if errors:
        joined = "\n  - ".join(errors)
        raise ValueError(f"Configuration validation failed:\n  - {joined}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Bulk RNA-seq workflow config.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config YAML.")
    args = parser.parse_args()

    with Path(args.config).open() as handle:
        cfg = yaml.safe_load(handle)
    samples, contrasts, counts_header = load_config_tables(cfg)
    validate_config(cfg, samples, contrasts, counts_header)
    print("Configuration validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
