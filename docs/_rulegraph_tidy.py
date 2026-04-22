"""Tidy Snakemake's --rulegraph dot output for README embedding.

Reads raw dot on stdin, writes cleaned dot to stdout:
  - drops the synthetic `all` target (star-shaped hub)
  - renames rule labels to match display terminology (TFEA, cMap)
  - overrides graph/node/edge styling (Helvetica, polyline splines, spacing)
"""

import re
import sys

LABEL_RENAMES = {
    "tf_activity": "TFEA",
    "drug_repositioning": "cMap",
    "exploratory_analysis": "exploratory",
    "qc_summary": "QC",
    "pathway_activity": "PROGENy",
    "gsea": "GSEA",
    "ora": "ORA",
    "deseq2": "DESeq2",
    "render_report": "Quarto report",
}


def main() -> None:
    src = sys.stdin.read()

    m = re.search(r'(\d+)\[label = "all"', src)
    if m:
        all_id = m.group(1)
        src = re.sub(rf'^\s*{all_id}\[label = "all".*\n', "", src, flags=re.M)
        src = re.sub(rf"^\s*\d+ -> {all_id}\s*\n", "", src, flags=re.M)
        src = re.sub(rf"^\s*{all_id} -> \d+\s*\n", "", src, flags=re.M)

    for old, new in LABEL_RENAMES.items():
        src = re.sub(rf'label = "{old}"', f'label = "{new}"', src)

    # Snakemake assigns each rule a per-node `color="H S V"` border hue plus
    # `style="rounded"`, which overrides the global node defaults and blocks
    # the fillcolor. Drop them so every node inherits the uniform styling.
    src = re.sub(r',\s*color = "[^"]*"', "", src)
    src = re.sub(r',\s*style="rounded"', "", src)

    src = re.sub(
        r"graph\[.*?\];\s*node\[.*?\];\s*edge\[.*?\];",
        (
            "graph[bgcolor=white, margin=0, rankdir=TB, splines=polyline, "
            "ranksep=0.9, nodesep=0.5, fontname=Helvetica];\n"
            '    node[shape=box, style="rounded,filled", fillcolor="#F5F5F5", '
            "fontname=Helvetica, fontsize=11, penwidth=1.2, margin=\"0.2,0.1\"];\n"
            "    edge[penwidth=1.2, color=\"#888888\", arrowsize=0.8, "
            "fontname=Helvetica];"
        ),
        src,
        count=1,
        flags=re.S,
    )

    sys.stdout.write(src)


if __name__ == "__main__":
    main()
