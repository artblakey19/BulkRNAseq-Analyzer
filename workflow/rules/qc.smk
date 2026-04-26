rule qc:
    input:
        counts = config["input"]["counts_tsv"],
        samples = config["input"]["samples_tsv"],
    output:
        summary = RESULTS / "qc" / "qc_summary.tsv",
    params:
        multiqc_dir = config["input"]["multiqc_data_dir"],
    log:
        "logs/qc/qc_summary.log",
    conda:
        "../envs/r-qc.yaml"
    script:
        "../scripts/qc.R"


rule exploratory:
    input:
        counts = config["input"]["counts_tsv"],
        samples = config["input"]["samples_tsv"],
    output:
        dds = RESULTS / "exploratory" / "dds_vst.rds",
        pca = RESULTS / "exploratory" / "pca.rds",
        cor = RESULTS / "exploratory" / "sample_correlation.rds",
        vst_matrix = RESULTS / "exploratory" / "vst_matrix.tsv",
    log:
        "logs/exploratory/exploratory.log",
    conda:
        "../envs/r-deseq2.yaml"
    script:
        "../scripts/exploratory.R"
