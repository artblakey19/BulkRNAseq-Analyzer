rule deseq2:
    input:
        counts = config["input"]["counts_tsv"],
        samples = config["input"]["samples_tsv"],
        contrasts = config["input"]["contrasts_tsv"],
    output:
        results = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
        rds = RESULTS / "de" / "{contrast}" / "deseq2_object.rds",
        summary = RESULTS / "de" / "{contrast}" / "deg_summary.tsv",
    params:
        contrast_id = "{contrast}",
        prefilter_min_count = config["de"]["prefilter_min_count"],
        lfc_shrink = config["de"]["lfc_shrink"],
        primary = config["de"]["primary"],
        secondary = config["de"]["secondary"],
    log:
        "logs/de/{contrast}.log",
    conda:
        "../envs/r-deseq2.yaml"
    script:
        "../scripts/deseq2.R"
