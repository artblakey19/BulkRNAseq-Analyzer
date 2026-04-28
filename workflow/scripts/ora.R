log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(dplyr)
})

set.seed(42)

# --- Inputs / outputs / params --------------------------------------------
de_path      <- snakemake@input[["de"]]
combined_out <- snakemake@output[["table"]]

primary   <- snakemake@params[["primary"]]
secondary <- snakemake@params[["secondary"]]
min_genes <- snakemake@params[["min_input_genes"]]
max_genes <- snakemake@params[["max_input_genes"]]
databases <- snakemake@params[["databases"]]

empty_df <- data.frame(
  database = character(),
  direction = character(),
  cutoff = character(),
  ID = character(),
  Description = character(),
  GeneRatio = character(),
  BgRatio = character(),
  pvalue = numeric(),
  p.adjust = numeric(),
  qvalue = numeric(),
  geneID = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# --- DEG Selection --------------------------------------------------------
de_res <- read.csv(de_path, stringsAsFactors = FALSE)
de_res <- de_res[!is.na(de_res$padj) & !is.na(de_res$log2FoldChange), ]

cutoff_used <- NA
selected_genes <- data.frame()

# Primary
prim_pass <- de_res[de_res$padj < primary$padj & abs(de_res$log2FoldChange) >= primary$abs_lfc, ]
if (nrow(prim_pass) >= min_genes) {
  cutoff_used <- "primary"
  selected_genes <- prim_pass
} else {
  # Secondary
  sec_pass <- de_res[de_res$padj < secondary$padj & abs(de_res$log2FoldChange) >= secondary$abs_lfc, ]
  if (nrow(sec_pass) >= min_genes) {
    cutoff_used <- "secondary"
    selected_genes <- sec_pass
  }
}

if (is.na(cutoff_used)) {
  message(sprintf("Insufficient DEG: neither primary nor secondary cutoff yielded >= %d genes.", min_genes))
  write.csv(empty_df, file = combined_out, row.names = FALSE)
  quit(save = "no", status = 0)
}

if (nrow(selected_genes) > max_genes) {
  message(sprintf("Total DEGs (%d) exceeds max_input_genes (%d). Truncating based on |log2FoldChange|.", nrow(selected_genes), max_genes))
  selected_genes <- selected_genes[order(abs(selected_genes$log2FoldChange), decreasing = TRUE), ]
  selected_genes <- head(selected_genes, max_genes)
}

up_genes <- selected_genes$gene_name[selected_genes$log2FoldChange > 0]
down_genes <- selected_genes$gene_name[selected_genes$log2FoldChange < 0]

up_genes <- up_genes[!is.na(up_genes) & up_genes != ""]
down_genes <- down_genes[!is.na(down_genes) & down_genes != ""]
universe_genes <- de_res$gene_name[!is.na(de_res$gene_name) & de_res$gene_name != ""]

message(sprintf("Using cutoff: %s. Up genes: %d, Down genes: %d", cutoff_used, length(up_genes), length(down_genes)))

# --- ORA Functions --------------------------------------------------------

run_msigdbr_ora <- function(genes, universe, collection, subcollection = NULL) {
  if (length(genes) == 0) return(NULL)
  m_df <- msigdbr(species = "Homo sapiens",
                  collection = collection,
                  subcollection = subcollection)
  m_t2g <- dplyr::select(m_df, gs_name, gene_symbol)
  enricher(gene = genes, TERM2GENE = m_t2g, universe = universe,
           pvalueCutoff = 1, qvalueCutoff = 1)
}

run_kegg_ora <- function(genes_symbol, universe_symbol) {
  if (length(genes_symbol) == 0) return(NULL)
  
  # Suppress messages from bitr
  gene_map <- suppressMessages(bitr(genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  dropped <- length(genes_symbol) - nrow(gene_map)
  if (dropped > 0) message("KEGG: dropped ", dropped, " genes during SYMBOL->ENTREZID mapping.")
  
  univ_map <- suppressMessages(bitr(universe_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  
  res <- enrichKEGG(gene = gene_map$ENTREZID, organism = "hsa", keyType = "kegg", universe = univ_map$ENTREZID, pvalueCutoff = 1, qvalueCutoff = 1)
  
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }
  return(res)
}

# --- Execute ORA ----------------------------------------------------------

all_res_dfs <- list()

for (db in databases) {
  for (dir in c("up", "down")) {
    genes <- if (dir == "up") up_genes else down_genes
    if (length(genes) == 0) next

    message(sprintf("Running ORA for %s (%s, n=%d)", db, dir, length(genes)))

    res <- NULL
    if (db == "GO_BP") {
      res <- run_msigdbr_ora(genes, universe_genes, collection = "C5", subcollection = "GO:BP")
    } else if (db == "Reactome") {
      res <- run_msigdbr_ora(genes, universe_genes, collection = "C2", subcollection = "CP:REACTOME")
    } else if (db == "Hallmark") {
      res <- run_msigdbr_ora(genes, universe_genes, collection = "H", subcollection = NULL)
    } else if (db == "KEGG") {
      res <- run_kegg_ora(genes, universe_genes)
    } else {
      message("Unknown database: ", db)
    }

    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      df <- as.data.frame(res)
      df$database <- db
      df$direction <- dir
      df$cutoff <- cutoff_used
      all_res_dfs[[paste(db, dir, sep = "_")]] <- df
    }
  }
}

if (length(all_res_dfs) > 0) {
  # bind_rows tolerates differing column sets — enrichKEGG / enricher can
  # return slightly different columns (e.g. geneID absent on zero-row results).
  final_df <- dplyr::bind_rows(all_res_dfs)
  rownames(final_df) <- NULL

  expected_cols <- c("database", "direction", "cutoff", "ID", "Description",
                     "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue",
                     "geneID", "Count")
  for (col in expected_cols) {
    if (!col %in% colnames(final_df)) final_df[[col]] <- NA
  }
  final_df <- final_df[, expected_cols]
} else {
  final_df <- empty_df
}

write.csv(final_df, file = combined_out, row.names = FALSE)

message("ORA analysis complete.")
