log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(dplyr)
})

# --- Inputs / outputs / params --------------------------------------------
de_path <- snakemake@input[["de"]]
combined_out <- snakemake@output[["table"]]
input_genes_out <- snakemake@output[["input_genes"]]

primary <- snakemake@params[["primary"]]
secondary <- snakemake@params[["secondary"]]
min_genes <- snakemake@params[["min_input_genes"]]
max_genes <- snakemake@params[["max_input_genes"]]
databases <- snakemake@params[["databases"]]

result_cols <- c(
  "database", "direction", "cutoff", "ID", "Description",
  "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue",
  "geneID", "Count"
)
input_cols <- c(
  "direction", "input_rank", "cutoff",
  "gene_id", "gene_name", "log2FoldChange", "padj"
)

empty_df <- function(cols) {
  setNames(data.frame(matrix(NA, nrow = 0, ncol = length(cols))), cols)
}

# Filter out Ensembl IDs in gene_name
is_valid_symbol <- function(x) {
  !is.na(x) & nzchar(x) & !grepl("^ENSG\\d+", x)
}

write_csv_safe <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, file = path, row.names = FALSE)
}

# --- DEG selection --------------------------------------------------------

apply_cutoff <- function(de, cut) {
  de[de$padj < cut$padj & abs(de$log2FoldChange) >= cut$abs_lfc, ]
}

# Returns list(genes, cutoff) or NULL when neither tier yields enough DEGs.
# Tries primary first, falls back to secondary; truncates by padj if oversized.
select_degs <- function(de, primary, secondary, min_genes, max_genes) {
  hits <- apply_cutoff(de, primary)
  cutoff <- "primary"
  if (nrow(hits) < min_genes) {
    hits <- apply_cutoff(de, secondary)
    cutoff <- "secondary"
  }
  if (nrow(hits) < min_genes) {
    return(NULL)
  }

  if (nrow(hits) > max_genes) {
    message(sprintf(
      "Total DEGs (%d) exceeds max_input_genes (%d). Truncating by padj.",
      nrow(hits), max_genes
    ))
    hits <- head(hits[order(hits$padj), ], max_genes)
  }
  hits <- hits[is_valid_symbol(hits$gene_name), ]
  list(genes = hits, cutoff = cutoff)
}

de_full <- read.csv(de_path, stringsAsFactors = FALSE)

# Universe must include every tested gene — DESeq2's independent filtering
# leaves padj=NA on low-count genes, but they were still part of the test set
# and belong in the hypergeometric background.
universe_genes <- de_full$gene_name[is_valid_symbol(de_full$gene_name)]

de_res <- de_full[!is.na(de_full$padj) & !is.na(de_full$log2FoldChange), ]
degs <- select_degs(de_res, primary, secondary, min_genes, max_genes)

if (is.null(degs)) {
  message(sprintf("Insufficient DEG: neither primary nor secondary cutoff yielded >= %d genes.", min_genes))
  write_csv_safe(empty_df(result_cols), combined_out)
  write_csv_safe(empty_df(input_cols), input_genes_out)
  quit(save = "no", status = 0)
}

input_genes_df <- degs$genes |>
  dplyr::mutate(
    direction = dplyr::case_when(
      log2FoldChange > 0 ~ "up",
      log2FoldChange < 0 ~ "down",
      TRUE ~ NA_character_
    ),
    cutoff = degs$cutoff
  ) |>
  dplyr::filter(direction %in% c("up", "down")) |>
  dplyr::group_by(direction) |>
  dplyr::arrange(dplyr::desc(abs(log2FoldChange)), .by_group = TRUE) |>
  dplyr::mutate(input_rank = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::select(dplyr::any_of(input_cols))

write_csv_safe(input_genes_df, input_genes_out)

up_genes <- input_genes_df$gene_name[input_genes_df$direction == "up"]
down_genes <- input_genes_df$gene_name[input_genes_df$direction == "down"]

message(sprintf(
  "Using cutoff: %s. Up genes: %d, Down genes: %d",
  degs$cutoff, length(up_genes), length(down_genes)
))

# --- ORA functions --------------------------------------------------------

# Parse a colon-separated MSigDB id like "C5:GO:BP" or "C2:CP:REACTOME" — first
# segment is the collection, the rest is the subcollection.
run_msigdbr_ora <- function(genes, universe, coll_str) {
  tokens <- strsplit(coll_str, ":", fixed = TRUE)[[1]]
  gene_sets <- msigdbr(
    species = "Homo sapiens",
    collection = tokens[1],
    subcollection = if (length(tokens) > 1) paste(tokens[-1], collapse = ":") else NULL
  )
  enricher(
    gene = genes, TERM2GENE = dplyr::select(gene_sets, gs_name, gene_symbol),
    universe = universe, pvalueCutoff = 1, qvalueCutoff = 1
  )
}

symbol_to_entrez <- function(symbols) {
  suppressMessages(bitr(symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID", OrgDb = org.Hs.eg.db
  ))
}

run_kegg_ora <- function(genes_symbol, universe_symbol) {
  genes_entrez <- symbol_to_entrez(genes_symbol)$ENTREZID
  universe_entrez <- symbol_to_entrez(universe_symbol)$ENTREZID

  dropped <- length(genes_symbol) - length(genes_entrez)
  if (dropped > 0) message("KEGG: dropped ", dropped, " genes during SYMBOL->ENTREZID mapping.")

  res <- enrichKEGG(
    gene = genes_entrez, organism = "hsa", keyType = "kegg",
    universe = universe_entrez, pvalueCutoff = 1, qvalueCutoff = 1
  )
  if (!is.null(res) && nrow(res@result) > 0) {
    res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  res
}

run_ora <- function(db, dir) {
  genes <- if (dir == "up") up_genes else down_genes
  if (length(genes) == 0) {
    return(NULL)
  }

  message(sprintf("Running ORA for %s (%s, n=%d)", db, dir, length(genes)))
  res <- if (db == "KEGG") {
    run_kegg_ora(genes, universe_genes)
  } else {
    run_msigdbr_ora(genes, universe_genes, db)
  }
  if (is.null(res) || nrow(res@result) == 0) {
    return(NULL)
  }

  df <- as.data.frame(res)
  df$database <- db
  df$direction <- dir
  df$cutoff <- degs$cutoff
  df
}

# --- Execute ORA ----------------------------------------------------------

jobs <- expand.grid(
  db = databases, dir = c("up", "down"),
  stringsAsFactors = FALSE
)
all_res_dfs <- Map(run_ora, jobs$db, jobs$dir)
all_res_dfs <- Filter(Negate(is.null), all_res_dfs)

if (length(all_res_dfs) > 0) {
  # bind_rows tolerates differing column sets — enrichKEGG / enricher can
  # return slightly different columns (e.g. geneID absent on zero-row results).
  final_df <- dplyr::bind_rows(all_res_dfs)
  rownames(final_df) <- NULL
  missing <- setdiff(result_cols, colnames(final_df))
  final_df[missing] <- NA
  final_df <- final_df[, result_cols]
} else {
  final_df <- empty_df(result_cols)
}

write_csv_safe(final_df, combined_out)

message("ORA analysis complete.")
