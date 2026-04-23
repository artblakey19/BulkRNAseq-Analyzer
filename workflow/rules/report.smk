rule cmap:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        hits = RESULTS / "cmap" / "{contrast}" / "l2s2_hits.tsv",
    params:
        top_up = config["cmap"]["top_up"],
        top_down = config["cmap"]["top_down"],
        service = config["cmap"]["service"],
    log:
        "logs/cmap/{contrast}.log",
    conda:
        "../envs/py-requests.yaml"
    script:
        "../scripts/l2s2_query.py"


rule render_report:
    input:
        qc_summary = RESULTS / "qc" / "qc_summary.tsv",
        pca = RESULTS / "exploratory" / "pca.rds",
        dds_vst = RESULTS / "exploratory" / "dds_vst.rds",
        de = expand(RESULTS / "de" / "{contrast}" / "deseq2_results.csv", contrast=CONTRAST_IDS),
        gsea = expand(RESULTS / "enrichment" / "{contrast}" / "gsea_combined.csv", contrast=CONTRAST_IDS),
        ora = expand(RESULTS / "enrichment" / "{contrast}" / "ora_combined.csv", contrast=CONTRAST_IDS),
        tf = expand(RESULTS / "tfea" / "{contrast}" / "tf_scores.tsv", contrast=CONTRAST_IDS),
        progeny = expand(RESULTS / "progeny" / "{contrast}" / "progeny_scores.tsv", contrast=CONTRAST_IDS),
        l2s2 = expand(RESULTS / "cmap" / "{contrast}" / "l2s2_hits.tsv", contrast=CONTRAST_IDS),
        # Resolved against the workflow source (repo root native, /app in
        # Docker) rather than the Snakemake workdir (/project bind-mount in
        # Docker), so the baked template is found regardless of CWD.
        template = str(Path(workflow.basedir).parent / "report" / "template.qmd"),
    output:
        html = RESULTS / "report" / "report.html",
    params:
        config_path = "config/config.yaml",
        formats = config["report"]["formats"],
    log:
        "logs/report/render.log",
    conda:
        "../envs/quarto.yaml"
    script:
        "../scripts/render_report.R"
