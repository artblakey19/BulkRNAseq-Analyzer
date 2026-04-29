GSEA_COLLECTION_IDS = [c["id"] for c in config["enrichment"]["gsea"]["collections"]]
ORA_DATABASE_IDS = [d["id"] for d in config["enrichment"]["ora"]["databases"]]


rule gsea:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "gsea_combined.csv",
    params:
        collections = GSEA_COLLECTION_IDS,
        ranking = config["enrichment"]["gsea"]["ranking"],
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
        input_genes = RESULTS / "enrichment" / "{contrast}" / "ora_input_genes.csv",
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


rule tfea:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        scores = RESULTS / "tfea" / "{contrast}" / "tf_scores.tsv",
        top = RESULTS / "tfea" / "{contrast}" / "tf_top.tsv",
    params:
        min_size = config["tfea"]["min_size"],
        padj_cutoff = config["tfea"]["padj_cutoff"],
        top_n = config["tfea"]["top_n"],
        split_complexes = config["tfea"]["split_complexes"],
    log:
        "logs/tfea/{contrast}.log",
    conda:
        "../envs/py-decoupler.yaml"
    script:
        "../scripts/tfea.py"


rule progeny:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
        vst_matrix = RESULTS / "exploratory" / "vst_matrix.tsv",
    output:
        scores = RESULTS / "progeny" / "{contrast}" / "progeny_scores.tsv",
    params:
        top_targets = config["progeny"]["top_targets"],
    log:
        "logs/progeny/{contrast}.log",
    conda:
        "../envs/py-decoupler.yaml"
    script:
        "../scripts/progeny.py"
