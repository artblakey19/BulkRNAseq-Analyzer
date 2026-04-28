log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(msigdbr)
})

set.seed(42)

# --- Params ---------------------------------------------------------------
ranking <- snakemake@params[["ranking"]]
colls <- snakemake@params[["collections"]]

# --- Input / Output -------------------------------------------------------
de_path <- snakemake@input[["de"]]
table_out <- snakemake@output[["table"]]

de_res <- read.csv(de_path, header = TRUE, stringsAsFactors = FALSE)

if (ranking == "stat") {
  ranks_df <- de_res %>%
    filter(!is.na(gene_name), !is.na(stat)) %>%
    group_by(gene_name) %>%
    summarize(metric = mean(stat)) %>%
    arrange(desc(metric))
} else if (ranking == "signed_p") {
  ranks_df <- de_res %>%
    filter(!is.na(gene_name), !is.na(log2FoldChange), !is.na(pvalue)) %>%
    mutate(metric = sign(log2FoldChange) * (-log10(pmax(pvalue, 1e-300)))) %>%
    filter(is.finite(metric)) %>%
    group_by(gene_name) %>%
    summarize(metric = mean(metric)) %>%
    arrange(desc(metric))
} else {
  stop("Unsupported ranking method: ", ranking)
}
ranks <- setNames(ranks_df$metric, ranks_df$gene_name)

# --- MSigDB Data Prep & GSEA ----------------------------------------------
get_pathways <- function(coll_str) {
  parts <- unlist(strsplit(coll_str, ":"))
  coll <- parts[1]
  if (length(parts) > 1) {
    subcoll <- paste(parts[-1], collapse = ":")
    m <- msigdbr(species = "Homo sapiens", collection = coll, subcollection = subcoll)
  } else {
    m <- msigdbr(species = "Homo sapiens", collection = coll)
  }
  split(m$gene_symbol, m$gs_name)
}

all_res <- list()

for (coll in colls) {
  message("Running GSEA for collection: ", coll)
  pathways <- get_pathways(coll)

  fgsea_res <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 15,
    maxSize = 500
  )

  if (nrow(fgsea_res) > 0) {
    res_df <- as.data.frame(fgsea_res)
    res_df$collection <- coll
    res_df$leadingEdge <- sapply(res_df$leadingEdge, paste, collapse = ";")
    all_res[[coll]] <- res_df
  }
}

col_order <- c("collection", "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")

if (length(all_res) > 0) {
  combined_df <- bind_rows(all_res)
  combined_df <- combined_df[, intersect(col_order, colnames(combined_df))]
} else {
  message("No significant enrichment found or no pathways evaluated.")
  combined_df <- data.frame(
    collection = character(), pathway = character(),
    pval = numeric(), padj = numeric(), log2err = numeric(),
    ES = numeric(), NES = numeric(), size = integer(),
    leadingEdge = character(), stringsAsFactors = FALSE
  )
}

dir.create(dirname(table_out), showWarnings = FALSE, recursive = TRUE)
write.csv(combined_df, table_out, row.names = FALSE)
