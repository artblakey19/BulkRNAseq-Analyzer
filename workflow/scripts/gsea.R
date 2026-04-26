log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(data.table)
})

# --- Params ---------------------------------------------------------------
seed     <- snakemake@params[["seed"]]
set.seed(seed)

ranking  <- snakemake@params[["ranking"]]
min_size <- snakemake@params[["min_size"]]
max_size <- snakemake@params[["max_size"]]
colls    <- snakemake@params[["collections"]]

# --- Input / Output -------------------------------------------------------
de_path    <- snakemake@input[["de"]]
table_out  <- snakemake@output[["table"]]
rds_out    <- snakemake@output[["rds"]]

# --- Load DE Results ------------------------------------------------------
de_res <- read.csv(de_path, header = TRUE, stringsAsFactors = FALSE)

# Generate Ranking Vector
if (ranking == "stat") {
  ranks_df <- de_res %>%
    filter(!is.na(gene_name), !is.na(stat)) %>%
    group_by(gene_name) %>%
    summarize(metric = mean(stat)) %>%
    arrange(desc(metric))
  
  ranks <- setNames(ranks_df$metric, ranks_df$gene_name)
} else if (ranking == "signed_p") {
  ranks_df <- de_res %>%
    filter(!is.na(gene_name), !is.na(log2FoldChange), !is.na(pvalue)) %>%
    mutate(
      metric = sign(log2FoldChange) * (-log10(pmax(pvalue, 1e-300)))
    ) %>%
    filter(is.finite(metric)) %>%
    group_by(gene_name) %>%
    summarize(metric = mean(metric)) %>%
    arrange(desc(metric))
    
  ranks <- setNames(ranks_df$metric, ranks_df$gene_name)
} else {
  stop("Unsupported ranking method: ", ranking)
}

# --- MSigDB Data Prep & GSEA ----------------------------------------------
# Function to parse collection string and fetch pathways
get_pathways <- function(coll_str) {
  parts   <- unlist(strsplit(coll_str, ":"))
  coll    <- parts[1]
  if (length(parts) > 1) {
    subcoll <- paste(parts[-1], collapse = ":")
    m <- msigdbr(species = "Homo sapiens", collection = coll, subcollection = subcoll)
  } else {
    m <- msigdbr(species = "Homo sapiens", collection = coll)
  }

  if (nrow(m) == 0) {
    warning("No MSigDB entries found for collection: ", coll_str)
    return(list())
  }
  split(m$gene_symbol, m$gs_name)
}

all_res <- list()
all_res_dt <- list()

for (coll in colls) {
  message("Running GSEA for collection: ", coll)
  pathways <- get_pathways(coll)
  if (length(pathways) == 0) next
  
  # run fgsea
  fgsea_res <- fgsea(pathways = pathways, 
                     stats    = ranks, 
                     minSize  = min_size, 
                     maxSize  = max_size)
  
  if (nrow(fgsea_res) > 0) {
    fgsea_res$collection <- coll
    all_res_dt[[coll]] <- fgsea_res
    
    # Format for CSV output
    res_df <- as.data.frame(fgsea_res)
    # Reorder columns to put collection first
    res_df <- res_df[, c("collection", setdiff(colnames(res_df), "collection"))]
    
    # Convert leadingEdge list to semicolon-separated string
    res_df$leadingEdge <- sapply(res_df$leadingEdge, paste, collapse = ";")
    
    all_res[[coll]] <- res_df
  }
}

if (length(all_res) > 0) {
  combined_df <- bind_rows(all_res)
  combined_dt <- rbindlist(all_res_dt, fill = TRUE)
  
  # Ensure columns exactly as specified: collection, pathway, pval, padj, log2err, ES, NES, size, leadingEdge
  col_order <- c("collection", "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")
  # Keep only columns present in the requested order (some metrics might be absent if fgsea failed to compute them, but usually they are present)
  col_order <- intersect(col_order, colnames(combined_df))
  combined_df <- combined_df[, col_order]
  
  # Ensure target directory exists
  dir.create(dirname(table_out), showWarnings = FALSE, recursive = TRUE)
  
  write.csv(combined_df, table_out, row.names = FALSE)
  saveRDS(combined_dt, rds_out)
} else {
  # Empty result
  message("No significant enrichment found or no pathways evaluated.")
  empty_df <- data.frame(
    collection = character(),
    pathway = character(),
    pval = numeric(),
    padj = numeric(),
    log2err = numeric(),
    ES = numeric(),
    NES = numeric(),
    size = integer(),
    leadingEdge = character(),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(table_out), showWarnings = FALSE, recursive = TRUE)
  write.csv(empty_df, table_out, row.names = FALSE)
  saveRDS(data.table::data.table(), rds_out)
}
