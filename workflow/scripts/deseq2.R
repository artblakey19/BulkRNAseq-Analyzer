# DESeq2 DE analysis per contrast, with apeglm-shrunken log2FoldChange.
# I/O driven by snakemake@input/@output/@params;
# Reference level is taken from contrasts.tsv `denominator`.
# Input: counts from `salmon.merged.gene_counts_length_scaled.tsv`

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

set.seed(42)

# --- Inputs / outputs / params --------------------------------------------
counts_path <- snakemake@input[["counts"]]
samples_path <- snakemake@input[["samples"]]
contrasts_path <- snakemake@input[["contrasts"]]

results_out <- snakemake@output[["results"]]
summary_out <- snakemake@output[["summary"]]

target_contrast <- snakemake@params[["contrast_id"]]
primary <- snakemake@params[["primary"]]
secondary <- snakemake@params[["secondary"]]

# --- Contrast row ---------------------------------------------------------
contrasts_df <- read.table(contrasts_path,
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE
)
contrast_row <- contrasts_df[contrasts_df$contrast_id == target_contrast, ,
  drop = FALSE
]
numerator <- as.character(contrast_row$numerator[1])
denominator <- as.character(contrast_row$denominator[1])

# --- Samples --------------------------------------------------------------
samples_df <- read.table(samples_path,
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE
)

# --- Counts ---------------------------------------------------------------
counts_raw <- read.table(counts_path,
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE
)
count_matrix <- as.matrix(counts_raw[, samples_df$sample, drop = FALSE])
rownames(count_matrix) <- counts_raw$gene_id
# DESeqDataSetFromMatrix() only accepts integer counts.
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

# --- colData + design -----------------------------------------------------
# ColData is metadata for samples, part of SummarizedExperiment object.
# Must be same order with colnames(countData)
# => rownames(colData) == colnames(countData)
col_data <- data.frame(
  sample = samples_df$sample,
  condition = factor(samples_df$condition),
  # 'replicate' is not used for analysis, but kept in colData as metadata.
  replicate = samples_df$replicate,
  batch = as.character(samples_df$batch),
  row.names = samples_df$sample,
  stringsAsFactors = FALSE
)

batch_unique <- unique(col_data$batch[!is.na(col_data$batch) & nzchar(col_data$batch)])
use_batch <- length(batch_unique) >= 2
if (use_batch) {
  col_data$batch <- factor(col_data$batch)
  design_formula <- ~ batch + condition
  message(sprintf("Design: ~ batch + condition  (batch levels: %d)", length(batch_unique)))
} else {
  design_formula <- ~condition
  message("Design: ~ condition  (batch has < 2 unique non-empty values)")
}

# --- Build DESeqDataSet ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = design_formula
)
# Set reference level for DESeq(). Otherwise, alphabetical order is used.
dds$condition <- relevel(dds$condition, ref = denominator)
# Add gene names to rowData(dds).
# rowData is metadata for rows(genes).
rowData(dds)$gene_name <- counts_raw$gene_name

# Prefilter: keep genes with count >= 10 in at least smallest-group samples.
# 10 is the value recommended for bulk RNA-seq in the DESeq2 vignette.
smallest_group <- min(table(col_data$condition))
keep <- rowSums(counts(dds) >= 10) >= smallest_group
message(sprintf(
  "Prefilter: keeping %d / %d genes (smallest_group=%d)",
  sum(keep), length(keep), smallest_group
))
dds <- dds[keep, ]

# --- Fit + shrink ---------------------------------------------------------
dds <- DESeq(dds)

# Check coefficient name(from contrasts.tsv) exists in dds(from samples.tsv)
# Manual editing of coefficient names in contrasts.tsv can cause this error.
coef_name <- paste0("condition_", numerator, "_vs_", denominator)
res_names <- resultsNames(dds)
if (!(coef_name %in% res_names)) {
  stop(sprintf(
    "Expected coefficient '%s' not in resultsNames(dds): [%s]",
    coef_name, paste(res_names, collapse = ", ")
  ))
}

# lfcShrink() shrinks the log2 fold changes to reduce the effect of noise
# in genes with low expression levels.
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

# apeglm does not return a `stat` column.
# Attach the unshrunk Wald statistic from results() manually.
# Wald stat is used for GSEA and TFEA.
res_shrunk$stat <- results(dds, name = coef_name)$stat

out_df <- as.data.frame(res_shrunk)
out_df <- data.frame(
  # Copy EnsemblID(gene_id) & Gene symbol(gene_name)
  # from rownames to the columns, CSV doesn't include rownames.
  gene_id = rownames(out_df),
  gene_name = rowData(dds)[rownames(out_df), "gene_name"],
  out_df,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
# Reorder columns
out_df <- out_df[, c(
  "gene_id", "gene_name", "baseMean", "log2FoldChange",
  "lfcSE", "stat", "pvalue", "padj"
)]
# Sort by adjusted padj
out_df <- out_df[order(out_df$padj, na.last = TRUE), ]

dir.create(dirname(results_out), showWarnings = FALSE, recursive = TRUE)
write.csv(out_df, file = results_out, row.names = FALSE)

# --- DEG cutoff summary ---------------------------------------------------
cutoff_counts <- function(df, padj_cut, lfc_cut) {
  valid <- !is.na(df$padj) & !is.na(df$log2FoldChange)
  sig <- valid & df$padj < padj_cut & abs(df$log2FoldChange) >= lfc_cut
  n_up <- sum(sig & df$log2FoldChange > 0)
  n_down <- sum(sig & df$log2FoldChange < 0)
  c(n_up = n_up, n_down = n_down, n_total = n_up + n_down)
}
prim <- cutoff_counts(out_df, primary$padj, primary$abs_lfc)
sec <- cutoff_counts(out_df, secondary$padj, secondary$abs_lfc)

summary_df <- data.frame(
  cutoff = c("primary", "secondary"),
  padj_cutoff = c(primary$padj, secondary$padj),
  lfc_cutoff = c(primary$abs_lfc, secondary$abs_lfc),
  n_up = c(prim[["n_up"]], sec[["n_up"]]),
  n_down = c(prim[["n_down"]], sec[["n_down"]]),
  n_total = c(prim[["n_total"]], sec[["n_total"]]),
  stringsAsFactors = FALSE
)
write.table(summary_df,
  file = summary_out, sep = "\t",
  quote = FALSE, row.names = FALSE
)

message("DESeq2 analysis complete for contrast: ", target_contrast)
