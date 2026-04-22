GSEA_COLLECTION_IDS = [c["id"] for c in config["enrichment"]["gsea"]["collections"]]
ORA_DATABASE_IDS = [d["id"] for d in config["enrichment"]["ora"]["databases"]]


rule gsea:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "gsea_combined.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "gsea_combined.rds",
    params:
        collections = GSEA_COLLECTION_IDS,
        ranking = config["enrichment"]["gsea"]["ranking"],
        min_size = config["enrichment"]["gsea"]["min_size"],
        max_size = config["enrichment"]["gsea"]["max_size"],
        seed = config["enrichment"]["gsea"]["seed"],
    log:
        "logs/enrichment/{contrast}_gsea.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/gsea.R"


rule ora:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "ora_combined.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "ora_combined.rds",
    params:
        databases = ORA_DATABASE_IDS,
        primary = config["de"]["primary"],
        secondary = config["de"]["secondary"],
        min_input_genes = config["enrichment"]["ora"]["min_input_genes"],
        max_input_genes = config["enrichment"]["ora"]["max_input_genes"],
    log:
        "logs/enrichment/{contrast}_ora.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/ora.R"


rule tf_activity:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        scores = RESULTS / "tf_activity" / "{contrast}" / "tf_scores.tsv",
        top = RESULTS / "tf_activity" / "{contrast}" / "tf_top.tsv",
    params:
        method = config["tf_activity"]["method"],
        min_size = config["tf_activity"]["min_size"],
        padj_cutoff = config["tf_activity"]["padj_cutoff"],
        top_n = config["tf_activity"]["top_n"],
        split_complexes = config["tf_activity"]["split_complexes"],
    log:
        "logs/tf_activity/{contrast}.log",
    conda:
        "../envs/py-decoupler.yaml"
    script:
        "../scripts/tf_activity.py"


rule pathway_activity:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
        vst_matrix = RESULTS / "exploratory" / "vst_matrix.tsv",
    output:
        scores = RESULTS / "pathway_activity" / "{contrast}" / "progeny_scores.tsv",
    params:
        top_targets = config["pathway_activity"]["top_targets"],
    log:
        "logs/pathway_activity/{contrast}.log",
    conda:
        "../envs/py-decoupler.yaml"
    script:
        "../scripts/pathway_activity.py"
