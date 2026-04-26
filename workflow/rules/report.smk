rule cmap:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        hits = RESULTS / "cmap" / "{contrast}" / "l2s2_hits.tsv",
    params:
        top_up = config["cmap"]["top_up"],
        top_down = config["cmap"]["top_down"],
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
        de = DE_OUTPUTS,
        gsea = GSEA_OUTPUTS,
        ora = ORA_OUTPUTS,
        tf = TFEA_OUTPUTS,
        progeny = PROGENY_OUTPUTS,
        l2s2 = CMAP_OUTPUTS,
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
