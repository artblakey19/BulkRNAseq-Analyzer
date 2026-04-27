# DESeq2 DE analysis per contrast, with apeglm-shrunken log2FoldChange.
# I/O is driven entirely by snakemake@input/@output/@params;
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
counts_path    <- snakemake@input[["counts"]]
samples_path   <- snakemake@input[["samples"]]
contrasts_path <- snakemake@input[["contrasts"]]

results_out <- snakemake@output[["results"]]
rds_out     <- snakemake@output[["rds"]]
summary_out <- snakemake@output[["summary"]]

target_contrast     <- snakemake@params[["contrast_id"]]
prefilter_min_count <- snakemake@params[["prefilter_min_count"]]
primary             <- snakemake@params[["primary"]]
secondary           <- snakemake@params[["secondary"]]

# --- Contrast row ---------------------------------------------------------
contrasts_df <- read.table(contrasts_path, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, check.names = FALSE)
contrast_row <- contrasts_df[contrasts_df$contrast_id == target_contrast,
                             , drop = FALSE]
numerator   <- as.character(contrast_row$numerator[1])
denominator <- as.character(contrast_row$denominator[1])

# --- Samples --------------------------------------------------------------
samples_df <- read.table(samples_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)

# --- Counts ---------------------------------------------------------------
counts_raw <- read.table(counts_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
gene_id_all   <- counts_raw[["gene_id"]]
gene_name_all <- counts_raw[["gene_name"]]
names(gene_name_all) <- gene_id_all

count_matrix <- as.matrix(counts_raw[, samples_df$sample, drop = FALSE])
rownames(count_matrix) <- gene_id_all
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

# --- colData + design -----------------------------------------------------
col_data <- data.frame(
  sample    = samples_df$sample,
  condition = factor(samples_df$condition),
  replicate = samples_df$replicate,
  batch     = as.character(samples_df$batch),
  row.names = samples_df$sample,
  stringsAsFactors = FALSE
)
col_data$condition <- relevel(col_data$condition, ref = denominator)

batch_unique <- unique(col_data$batch[!is.na(col_data$batch) & nzchar(col_data$batch)])
use_batch <- length(batch_unique) >= 2
if (use_batch) {
  col_data$batch <- factor(col_data$batch)
  design_formula <- ~ batch + condition
  message(sprintf("Design: ~ batch + condition  (batch levels: %d)", length(batch_unique)))
} else {
  design_formula <- ~ condition
  message("Design: ~ condition  (batch has < 2 unique non-empty values)")
}

# --- Build DESeqDataSet ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData   = col_data,
                              design    = design_formula)
rowData(dds)$gene_name <- gene_name_all[rownames(dds)]

# Prefilter: keep genes with >= min_count in at least smallest_group samples.
smallest_group <- min(table(col_data$condition))
keep <- rowSums(counts(dds) >= prefilter_min_count) >= smallest_group
message(sprintf("Prefilter: keeping %d / %d genes (smallest_group=%d, min_count=%d)",
                sum(keep), length(keep), smallest_group, prefilter_min_count))
dds <- dds[keep, ]

# --- Fit + shrink ---------------------------------------------------------
# DESeq2's parametric / mean dispersion fits can fail on very small or
# low-variance datasets. DESeq2 docs recommend falling back to gene-wise
# dispersion estimates in that case.
fit_with_fallback <- function(dds) {
  tryCatch(
    DESeq(dds),
    error = function(e) {
      if (!grepl("dispersion", conditionMessage(e), ignore.case = TRUE)) stop(e)
      message("Dispersion curve fit failed; falling back to gene-wise estimates. ",
              "Original error: ", conditionMessage(e))
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersionsGeneEst(dds)
      dispersions(dds) <- mcols(dds)$dispGeneEst
      nbinomWaldTest(dds)
    }
  )
}
dds <- fit_with_fallback(dds)

coef_name <- paste0("condition_", numerator, "_vs_", denominator)
res_names <- resultsNames(dds)
if (!(coef_name %in% res_names)) {
  stop(sprintf("Expected coefficient '%s' not in resultsNames(dds): [%s]",
               coef_name, paste(res_names, collapse = ", ")))
}

res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
res_wald   <- results(dds, name = coef_name)

# apeglm does not return a `stat` column; use the Wald test statistic from
# results() so the frozen output schema is satisfied and downstream GSEA
# (ranking="stat") has a ranking vector.
stat_vec <- res_wald$stat[match(rownames(res_shrunk), rownames(res_wald))]

out_df <- data.frame(
  gene_id        = rownames(res_shrunk),
  gene_name      = unname(gene_name_all[rownames(res_shrunk)]),
  baseMean       = res_shrunk$baseMean,
  log2FoldChange = res_shrunk$log2FoldChange,
  lfcSE          = res_shrunk$lfcSE,
  stat           = stat_vec,
  pvalue         = res_shrunk$pvalue,
  padj           = res_shrunk$padj,
  stringsAsFactors = FALSE,
  check.names    = FALSE
)
out_df <- out_df[order(out_df$padj, na.last = TRUE), ]

dir.create(dirname(results_out), showWarnings = FALSE, recursive = TRUE)
write.csv(out_df, file = results_out, row.names = FALSE)
saveRDS(dds, file = rds_out)

# --- DEG cutoff summary ---------------------------------------------------
cutoff_counts <- function(df, padj_cut, lfc_cut) {
  valid <- !is.na(df$padj) & !is.na(df$log2FoldChange)
  sig   <- valid & df$padj < padj_cut & abs(df$log2FoldChange) >= lfc_cut
  n_up   <- sum(sig & df$log2FoldChange > 0)
  n_down <- sum(sig & df$log2FoldChange < 0)
  c(n_up = n_up, n_down = n_down, n_total = n_up + n_down)
}
prim <- cutoff_counts(out_df, primary$padj,   primary$abs_lfc)
sec  <- cutoff_counts(out_df, secondary$padj, secondary$abs_lfc)

summary_df <- data.frame(
  cutoff      = c("primary", "secondary"),
  padj_cutoff = c(primary$padj,    secondary$padj),
  lfc_cutoff  = c(primary$abs_lfc, secondary$abs_lfc),
  n_up        = c(prim[["n_up"]],    sec[["n_up"]]),
  n_down      = c(prim[["n_down"]],  sec[["n_down"]]),
  n_total     = c(prim[["n_total"]], sec[["n_total"]]),
  stringsAsFactors = FALSE
)
write.table(summary_df, file = summary_out, sep = "\t",
            quote = FALSE, row.names = FALSE)

message("DESeq2 analysis complete for contrast: ", target_contrast)
