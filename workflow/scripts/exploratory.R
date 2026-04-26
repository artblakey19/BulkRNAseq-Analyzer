# workflow/scripts/exploratory.R
# B2: Exploratory analysis — build VST DESeqTransform, PCA on top-500
# variable genes, sample-sample Spearman/Pearson correlation, hierarchical
# clustering on 1 - Spearman.

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(matrixStats)
})

set.seed(42)

counts_path  <- snakemake@input[["counts"]]
samples_path <- snakemake@input[["samples"]]

samples <- read.table(samples_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)

counts_raw <- read.table(counts_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
gene_ids   <- as.character(counts_raw[[1]])
gene_names <- as.character(counts_raw[[2]])

counts_mat <- round(as.matrix(counts_raw[, -c(1, 2), drop = FALSE]))
mode(counts_mat) <- "integer"
rownames(counts_mat) <- gene_ids
counts_mat <- counts_mat[, samples$sample, drop = FALSE]

col_data <- data.frame(
  sample    = samples$sample,
  condition = factor(samples$condition),
  replicate = samples$replicate,
  batch     = factor(samples$batch),
  row.names = samples$sample,
  stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = col_data,
                              design    = ~ condition)
rowData(dds)$gene_name <- gene_names[match(rownames(dds), gene_ids)]

vst_obj <- tryCatch(
  vst(dds, blind = TRUE),
  error = function(e) {
    message("vst() failed (", conditionMessage(e),
            "); falling back to varianceStabilizingTransformation()")
    varianceStabilizingTransformation(dds, blind = TRUE)
  }
)
vst_mat <- assay(vst_obj)

row_vars <- rowVars(vst_mat)
n_top <- min(500L, length(row_vars))
top_idx <- order(row_vars, decreasing = TRUE)[seq_len(n_top)]
top_genes <- rownames(vst_mat)[top_idx]

pca <- prcomp(t(vst_mat[top_idx, , drop = FALSE]),
              center = TRUE, scale. = FALSE)
var_explained <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
names(var_explained) <- paste0("PC", seq_along(var_explained))

scores_df <- as.data.frame(pca$x)
meta_match <- match(rownames(scores_df), rownames(col_data))
scores_df$sample    <- rownames(scores_df)
scores_df$condition <- as.character(col_data$condition)[meta_match]
scores_df$replicate <- col_data$replicate[meta_match]
scores_df$batch     <- as.character(col_data$batch)[meta_match]
pc_cols <- grep("^PC", colnames(scores_df), value = TRUE)
scores_df <- scores_df[, c("sample", "condition", "replicate", "batch", pc_cols)]
rownames(scores_df) <- scores_df$sample

pca_out <- list(
  scores        = scores_df,
  var_explained = var_explained,
  top_genes     = top_genes
)

# Correlate on the top-variable gene set used by PCA (bulk convention, DESeq2
# default). Using all ~tens-of-thousands of genes makes Spearman dominated by
# the baseline rank structure of housekeeping / low-variance genes, collapsing
# within- vs between-group correlations into a narrow band. Restricting to the
# top variable genes produces correlations that track the biology and match
# what PCA shows.
vst_top <- vst_mat[top_idx, , drop = FALSE]
spearman_mat <- cor(vst_top, method = "spearman")
pearson_mat  <- cor(vst_top, method = "pearson")
hc <- hclust(as.dist(1 - spearman_mat), method = "average")

cor_out <- list(
  spearman = spearman_mat,
  pearson  = pearson_mat,
  hclust   = hc
)

saveRDS(vst_obj, snakemake@output[["dds"]])
saveRDS(pca_out, snakemake@output[["pca"]])
saveRDS(cor_out, snakemake@output[["cor"]])

vst_out_path <- snakemake@output[["vst_matrix"]]
if (!is.null(vst_out_path) && nzchar(vst_out_path)) {
  vst_df <- data.frame(
    gene_id   = rownames(vst_mat),
    gene_name = gene_names[match(rownames(vst_mat), gene_ids)],
    vst_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write.table(vst_df, file = vst_out_path, sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
}

message(sprintf("exploratory.R done: samples=%d genes=%d top_genes=%d",
                ncol(vst_mat), nrow(vst_mat), n_top))
