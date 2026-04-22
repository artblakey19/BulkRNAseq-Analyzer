# Per-sample QC metric table.
# Joins samples.tsv metadata, MultiQC general stats, and counts-derived
# assigned_reads. Metrics-only: no pass/fail columns, no thresholds.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

counts_path  <- snakemake@input[["counts"]]
multiqc_dir  <- snakemake@params[["multiqc_dir"]]
samples_path <- snakemake@input[["samples"]]
out_path     <- snakemake@output[["summary"]]

stopifnot(file.exists(counts_path))
stopifnot(file.exists(samples_path))

REQUIRED_COLS <- c(
  "sample", "condition", "replicate", "batch",
  "total_reads", "mapping_rate", "rrna_pct", "duplication", "assigned_reads"
)

samples <- read_tsv(samples_path, show_col_types = FALSE)
stopifnot("sample" %in% names(samples))
for (col in c("condition", "replicate", "batch")) {
  if (!col %in% names(samples)) samples[[col]] <- NA
}

counts <- read_tsv(counts_path, show_col_types = FALSE)
sample_cols <- setdiff(names(counts), c("gene_id", "gene_name"))

missing_in_counts <- setdiff(samples$sample, sample_cols)
if (length(missing_in_counts) > 0) {
  stop(sprintf(
    "qc_summary: samples in samples.tsv missing from counts matrix: %s",
    paste(missing_in_counts, collapse = ", ")
  ))
}

assigned_tbl <- tibble(
  sample = samples$sample,
  assigned_reads = vapply(
    samples$sample,
    function(s) sum(counts[[s]], na.rm = TRUE),
    numeric(1)
  )
)

# Pull a metric column from a general_stats-style data frame by trying a set
# of alias names. Matches exact or suffix (for module-prefixed MultiQC keys
# such as "STAR_mqc-generalstats-star-uniquely_mapped_percent").
metric_aliases <- list(
  total_reads  = c("total_reads", "total_sequences",
                   "m_reads", "m_seqs", "m_reads_mapped"),
  mapping_rate = c("mapping_rate", "pct_mapped", "percent_mapped",
                   "uniquely_mapped_percent"),
  rrna_pct     = c("rrna_pct", "rrna", "percent_rrna", "pct_rrna"),
  duplication  = c("duplication", "percent_duplication",
                   "pct_duplication", "percent_duplicates")
)

normalise_name <- function(x) gsub("[^a-z0-9]", "_", tolower(x))

pick_metric <- function(df, aliases) {
  norm_names <- normalise_name(names(df))
  for (a in aliases) {
    key <- normalise_name(a)
    hit <- which(norm_names == key)
    if (length(hit) > 0) return(df[[hit[1]]])
    hit <- grep(paste0("(^|_)", key, "$"), norm_names)
    if (length(hit) > 0) return(df[[hit[1]]])
  }
  rep(NA_real_, nrow(df))
}

load_multiqc_general_stats <- function(dir) {
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(NULL)
  candidates <- c(
    file.path(dir, "multiqc_general_stats.txt"),
    file.path(dir, "multiqc_general_stats.tsv"),
    file.path(dir, "general_stats.tsv"),
    file.path(dir, "general_stats.txt")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) return(NULL)
  df <- read_tsv(found[1], show_col_types = FALSE)
  if (ncol(df) < 2) return(NULL)
  names(df)[1] <- "sample"
  df$sample <- as.character(df$sample)
  df
}

mqc <- load_multiqc_general_stats(multiqc_dir)

if (!is.null(mqc)) {
  mqc_tbl <- tibble(
    sample       = mqc$sample,
    total_reads  = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$total_reads))),
    mapping_rate = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$mapping_rate))),
    rrna_pct     = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$rrna_pct))),
    duplication  = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$duplication)))
  )
  message(sprintf("qc_summary: loaded MultiQC general stats (%d samples)", nrow(mqc_tbl)))
} else {
  message("qc_summary: no MultiQC general_stats file found; QC metrics set to NA")
  mqc_tbl <- tibble(
    sample       = samples$sample,
    total_reads  = NA_real_,
    mapping_rate = NA_real_,
    rrna_pct     = NA_real_,
    duplication  = NA_real_
  )
}

summary_tbl <- samples %>%
  select(sample, condition, replicate, batch) %>%
  left_join(mqc_tbl, by = "sample") %>%
  left_join(assigned_tbl, by = "sample")

for (col in REQUIRED_COLS) {
  if (!col %in% names(summary_tbl)) summary_tbl[[col]] <- NA
}
summary_tbl <- summary_tbl[, REQUIRED_COLS]

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write_tsv(summary_tbl, out_path)

message(sprintf("qc_summary: wrote %d rows x %d cols -> %s",
                nrow(summary_tbl), ncol(summary_tbl), out_path))
