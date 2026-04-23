#!/usr/bin/env python3
"""Interactive init CLI that generates config.yaml / samples.tsv / contrasts.tsv.

Reads the count matrix header to auto-extract sample names, then prompts for
project metadata, per-sample condition/replicate/batch, and contrast definitions.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from textwrap import dedent

import yaml

# PROJECT_ROOT env var lets the Docker wrapper point REPO_ROOT at the
# bind-mounted /project instead of the baked-in /app, so relative counts
# paths resolve against the host project dir. The config template is a
# pipeline asset that lives next to the script (/app in Docker, repo root
# native), so it's resolved separately from REPO_ROOT.
REPO_ROOT = Path(os.environ.get("PROJECT_ROOT") or Path(__file__).resolve().parents[2])
SCRIPT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "config" / "config.yaml"
DEFAULT_SAMPLES = REPO_ROOT / "config" / "samples.tsv"
DEFAULT_CONTRASTS = REPO_ROOT / "config" / "contrasts.tsv"
CONFIG_TEMPLATE = SCRIPT_ROOT / "config" / "config.template.yaml"


def prompt(message: str, default: str | None = None) -> str:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        answer = input(f"{message}{suffix}: ").strip()
        if answer:
            return answer
        if default is not None:
            return default
        print("  (required)")


def prompt_choice(message: str, choices: list[str], default: str | None = None) -> str:
    choice_str = "/".join(choices)
    while True:
        answer = prompt(f"{message} ({choice_str})", default=default)
        if answer in choices:
            return answer
        print(f"  pick one of: {choice_str}")


def read_sample_names(counts_path: Path) -> list[str]:
    with counts_path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
    if len(header) < 3:
        raise SystemExit(f"counts TSV header has <3 cols: {header}")
    return header[2:]


def load_existing_config(path: Path) -> dict:
    # Prefer the user's existing config (re-running init preserves prior values);
    # fall back to the committed template for a first-time run.
    source = path if path.exists() else CONFIG_TEMPLATE
    if not source.exists():
        raise SystemExit(f"config template missing: {CONFIG_TEMPLATE}")
    with source.open() as f:
        return yaml.safe_load(f)


def collect_project_meta(existing: dict) -> dict:
    proj = existing.get("project", {})
    print("\n=== Project metadata ===")
    proj["name"] = prompt("Project name", default=proj.get("name"))
    proj["analyst"] = prompt("Analyst name", default=proj.get("analyst"))

    print("\n--- Experiment (reported only, not used in analysis) ---")
    exp = proj.get("experiment", {}) or {}
    exp["cell_line"] = prompt("Cell line", default=exp.get("cell_line") or None)
    exp["compound"] = prompt("Compound / treatment", default=exp.get("compound") or None)
    exp["dose"] = prompt("Dose (e.g. 10 uM)", default=exp.get("dose") or None)
    exp["duration"] = prompt("Duration (e.g. 24 h)", default=exp.get("duration") or None)
    exp["notes"] = prompt("Notes (optional)", default=exp.get("notes") or "")
    proj["experiment"] = exp
    return proj


def collect_input_paths(existing: dict) -> dict:
    inp = existing.get("input", {})
    print("\n=== Input paths ===")
    inp["counts_tsv"] = prompt(
        "Path to salmon.merged.gene_counts_length_scaled.tsv",
        default=inp.get("counts_tsv"),
    )
    inp["multiqc_data_dir"] = prompt(
        "Path to multiqc_data/ directory",
        default=inp.get("multiqc_data_dir"),
    )
    inp["samples_tsv"] = inp.get("samples_tsv", "config/samples.tsv")
    inp["contrasts_tsv"] = inp.get("contrasts_tsv", "config/contrasts.tsv")
    return inp


def collect_sample_table(sample_ids: list[str]) -> list[dict]:
    print(f"\n=== Sample metadata ({len(sample_ids)} samples detected) ===")
    print("Samples from counts header:")
    for s in sample_ids:
        print(f"  - {s}")

    print(dedent("""
        For each sample, enter condition, replicate number, batch.
        Tip: press Enter on condition to reuse the previous sample's condition.
    """).strip())

    rows: list[dict] = []
    prev_condition: str | None = None
    replicate_counter: dict[str, int] = {}
    for sid in sample_ids:
        print(f"\n  [{sid}]")
        cond_default = prev_condition
        condition = prompt("    condition", default=cond_default)
        prev_condition = condition
        replicate_counter[condition] = replicate_counter.get(condition, 0) + 1
        replicate = prompt("    replicate", default=str(replicate_counter[condition]))
        batch = prompt("    batch", default="1")
        rows.append(
            {"sample": sid, "condition": condition, "replicate": replicate, "batch": batch}
        )
    return rows


def collect_contrasts(conditions: list[str]) -> list[dict]:
    unique_conds = sorted(set(conditions))
    print(f"\n=== Contrasts (unique conditions: {unique_conds}) ===")
    rows: list[dict] = []
    idx = 0
    while True:
        more = prompt_choice(
            "Add a contrast?" if idx == 0 else "Add another contrast?",
            ["y", "n"],
            default="y" if idx == 0 else "n",
        )
        if more == "n":
            if idx == 0:
                print("  at least one contrast is required")
                continue
            break
        print(f"\n  Contrast #{idx + 1}")
        numerator = prompt("    numerator (treatment) condition")
        denominator = prompt("    denominator (reference) condition")
        if numerator not in unique_conds or denominator not in unique_conds:
            print(f"  both must be one of {unique_conds}")
            continue
        if numerator == denominator:
            print("  numerator and denominator must differ")
            continue
        default_id = f"{numerator}_vs_{denominator}"
        contrast_id = prompt("    contrast_id", default=default_id)
        description = prompt("    description", default=f"{numerator} vs {denominator}")
        rows.append(
            {
                "contrast_id": contrast_id,
                "factor": "condition",
                "numerator": numerator,
                "denominator": denominator,
                "description": description,
            }
        )
        idx += 1
    return rows


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def confirm_overwrite(paths: list[Path]) -> bool:
    existing = [p for p in paths if p.exists()]
    if not existing:
        return True
    print("\nThe following files will be overwritten:")
    for p in existing:
        print(f"  - {p.relative_to(REPO_ROOT)}")
    return prompt_choice("Proceed?", ["y", "n"], default="n") == "y"


def main() -> int:
    parser = argparse.ArgumentParser(description="Interactive project init")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--samples", type=Path, default=DEFAULT_SAMPLES)
    parser.add_argument("--contrasts", type=Path, default=DEFAULT_CONTRASTS)
    parser.add_argument(
        "--counts",
        type=Path,
        help="Override counts TSV path (otherwise read from config or prompted)",
    )
    args = parser.parse_args()

    existing = load_existing_config(args.config)
    existing["project"] = collect_project_meta(existing)
    existing["input"] = collect_input_paths(existing)

    counts_path = Path(args.counts) if args.counts else Path(existing["input"]["counts_tsv"])
    if not counts_path.is_absolute():
        counts_path = REPO_ROOT / counts_path
    if not counts_path.exists():
        print(f"\nWARNING: counts file not found at {counts_path}")
        print("  You can still proceed but sample names must be entered manually.")
        manual = prompt_choice("Enter sample names manually?", ["y", "n"], default="n")
        if manual == "y":
            raw = prompt("Comma-separated sample names")
            sample_ids = [s.strip() for s in raw.split(",") if s.strip()]
        else:
            return 1
    else:
        sample_ids = read_sample_names(counts_path)

    sample_rows = collect_sample_table(sample_ids)
    contrast_rows = collect_contrasts([r["condition"] for r in sample_rows])

    if not confirm_overwrite([args.config, args.samples, args.contrasts]):
        print("aborted.")
        return 1

    with args.config.open("w") as f:
        yaml.safe_dump(existing, f, sort_keys=False, default_flow_style=False)
    write_tsv(args.samples, sample_rows, ["sample", "condition", "replicate", "batch"])
    write_tsv(
        args.contrasts,
        contrast_rows,
        ["contrast_id", "factor", "numerator", "denominator", "description"],
    )

    print("\nDone.")
    print(f"  {args.config.relative_to(REPO_ROOT)}")
    print(f"  {args.samples.relative_to(REPO_ROOT)}")
    print(f"  {args.contrasts.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
