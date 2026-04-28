# Exploratory analysis — build VST DESeqTransform, PCA on top-500
# variable genes, sample-sample Spearman correlation, hierarchical
# clustering on 1 - Spearman.
# Input: counts from `salmon.merged.gene_counts_length_scaled.tsv`
# Output: vst_deseq.rds, pca_results.rds, correlations.rds, vst_matrix.tsv

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(matrixStats)
})

counts_path <- snakemake@input[["counts"]]
samples_path <- snakemake@input[["samples"]]

samples <- read.table(samples_path,
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE
)

counts_raw <- read.table(counts_path,
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE
)
# Convert count_raw(DataFrame) to count_matrix(matrix).
count_matrix <- as.matrix(counts_raw[, samples$sample, drop = FALSE])
rownames(count_matrix) <- counts_raw$gene_id
# DESeqDataSetFromMatrix() only accepts integer counts.
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

col_data <- data.frame(
  sample = samples$sample,
  condition = factor(samples$condition),
  replicate = samples$replicate,
  batch = factor(samples$batch),
  row.names = samples$sample,
  stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  # 'design' is ignored for VST(blind=TRUE), but DESeqDataSetFromMatrix()
  # requires it, so just add a placeholder.
  design = ~condition
)
# Add gene_name(gene symbol) to rowData.
rowData(dds)$gene_name <- counts_raw$gene_name

# DESeq2 uses raw count data for differential expression estimation.
# Variance stabilizing transformation is only for visualization.
vst_obj <- tryCatch(
  vst(dds, blind = TRUE),
  # vst() is faster version of varianceStabilizingTransformation()
  # But it may fail for small samples like toy dataset.
  error = function(e) {
    message(
      "vst() failed (", conditionMessage(e),
      "); falling back to varianceStabilizingTransformation()"
    )
    varianceStabilizingTransformation(dds, blind = TRUE)
  }
)
vst_mat <- assay(vst_obj)

row_vars <- rowVars(vst_mat)
n_top <- min(500L, length(row_vars))
top_idx <- order(row_vars, decreasing = TRUE)[seq_len(n_top)]

pca <- prcomp(t(vst_mat[top_idx, , drop = FALSE]))
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
names(var_explained) <- paste0("PC", seq_along(var_explained))

# 'pca' from prcomp() is simple
scores_df <- as.data.frame(pca$x)
meta_match <- match(rownames(scores_df), rownames(col_data))
scores_df$sample <- rownames(scores_df)
scores_df$condition <- as.character(col_data$condition)[meta_match]
scores_df$replicate <- col_data$replicate[meta_match]
scores_df$batch <- as.character(col_data$batch)[meta_match]
# Reorder columns: meta first, PC scores after.
pc_cols <- colnames(pca$x)
scores_df <- scores_df[, c("sample", "condition", "replicate", "batch", pc_cols)]
rownames(scores_df) <- scores_df$sample

pca_out <- list(
  scores        = scores_df,
  var_explained = var_explained
)

# Sample-to-sample Euclidean distances on VST-transformed counts.
# This may have a problem with curse of dimensionality, but it's common practice.
# DESeq2 vignette also uses this approach.
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization
sample_dist <- dist(t(vst_mat))
sample_dist_mat <- as.matrix(sample_dist)
hc <- hclust(sample_dist)

dist_out <- list(
  dist   = sample_dist_mat,
  hclust = hc
)

saveRDS(vst_obj, snakemake@output[["dds"]])
saveRDS(pca_out, snakemake@output[["pca"]])
saveRDS(dist_out, snakemake@output[["cor"]])

# Save VST matrix as .tsv
vst_df <- data.frame(
  gene_id = rownames(vst_mat),
  gene_name = rowData(vst_obj)$gene_name,
  vst_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.table(vst_df,
  file = snakemake@output[["vst_matrix"]], sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA"
)

message(sprintf(
  "exploratory.R done: samples=%d genes=%d top_genes=%d",
  ncol(vst_mat), nrow(vst_mat), n_top
))