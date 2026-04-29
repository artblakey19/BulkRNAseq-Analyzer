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

REQUIRED_COLS <- c(
  "sample", "condition", "replicate", "batch",
  "total_reads", "mapping_rate", "uniquely_mapped_pct",
  "rrna_pct", "duplication", "dupradar_slope",
  "exonic_pct", "intronic_pct", "intergenic_pct",
  "gene_body_5_3_bias", "insert_size_median", "gc_content_pct",
  "strandedness",
  "assigned_reads"
)

samples <- read_tsv(samples_path, show_col_types = FALSE)

counts <- read_tsv(counts_path, show_col_types = FALSE)

missing_counts <- setdiff(samples$sample, names(counts))
if (length(missing_counts) > 0) {
  stop(sprintf("samples.tsv has samples not present as columns in counts: %s",
               paste(missing_counts, collapse = ", ")))
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
  total_reads         = c("total_reads", "total_sequences"),
  # Generic "overall % aligned" — from Salmon/STAR aggregate columns. Includes
  # multi-mappers, so it can look healthy even when unique alignments are low.
  mapping_rate        = c("mapping_rate", "pct_mapped", "percent_mapped"),
  # STAR unique alignment rate — the DE-relevant fraction. Big gap between
  # this and mapping_rate flags repeat/pseudo-gene multi-mapping issues.
  uniquely_mapped_pct = c("uniquely_mapped_percent", "uniquely_mapped_pct",
                          "pct_uniquely_mapped"),
  rrna_pct            = c("rrna_pct", "rrna", "percent_rrna", "pct_rrna"),
  duplication         = c("duplication", "percent_duplication",
                          "pct_duplication", "percent_duplicates"),
  # FastQC %GC. Outliers signal foreign-organism contamination, adapter dimers,
  # or unbalanced rRNA depletion. Tissue-typical GC clusters are tight, so a
  # single sample sitting >5pp from the cohort mean is a red flag.
  gc_content_pct      = c("percent_gc", "gc_pct", "gc_content", "avg_gc"),
  # samtools stats / Picard CollectInsertSizeMetrics. Short median insert size
  # = over-fragmentation or adapter read-through; very long = under-fragmented
  # library or PCR chimera.
  insert_size_median  = c("insert_size_median", "median_insert_size",
                          "insert_size_average", "avg_insert_size")
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

load_mqc_module <- function(dir, basename) {
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(NULL)
  path <- file.path(dir, basename)
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || ncol(df) < 2) return(NULL)
  names(df)[1] <- "sample"
  df$sample <- as.character(df$sample)
  df
}

load_multiqc_general_stats <- function(dir) {
  for (bn in c("multiqc_general_stats.txt", "multiqc_general_stats.tsv")) {
    df <- load_mqc_module(dir, bn)
    if (!is.null(df)) return(df)
  }
  NULL
}

numeric_col <- function(df, pattern) {
  hit <- grep(paste0("^", pattern, "$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) return(rep(NA_real_, nrow(df)))
  suppressWarnings(as.numeric(df[[hit[1]]]))
}

# Strandedness call from RSeQC infer_experiment. Used as a sanity-check column
# only — dUTP chemistry is deterministic, so this is uniform across a properly
# prepared cohort. Unexpected forward/unstranded calls flag a kit mix-up, not
# a biology effect.
parse_infer_experiment <- function(dir) {
  df <- load_mqc_module(dir, "multiqc_rseqc_infer_experiment.txt")
  if (is.null(df)) return(NULL)
  merge_cols <- function(df, suffix) {
    pe <- paste0("pe_", suffix)
    se <- paste0("se_", suffix)
    if (pe %in% names(df) && se %in% names(df)) {
      return(ifelse(is.na(df[[pe]]) | df[[pe]] == 0, df[[se]], df[[pe]]))
    }
    if (pe %in% names(df)) return(df[[pe]])
    if (se %in% names(df)) return(df[[se]])
    rep(NA_real_, nrow(df))
  }
  sense <- suppressWarnings(as.numeric(merge_cols(df, "sense")))
  anti  <- suppressWarnings(as.numeric(merge_cols(df, "antisense")))
  und   <- suppressWarnings(as.numeric(merge_cols(df, "undetermined")))
  call <- mapply(function(s, a, u) {
    vals <- c(forward = s, reverse = a, unstranded = u)
    if (all(is.na(vals))) return(NA_character_)
    names(which.max(vals))
  }, sense, anti, und, USE.NAMES = FALSE)
  tibble(sample = df$sample, strandedness = call)
}

# RSeQC read_distribution. exonic = CDS + 5'UTR + 3'UTR. Intronic and
# intergenic percentages flag distinct failures: high intronic = pre-mRNA
# (poly-A enrichment failure or nuclear RNA dominance); high intergenic =
# genomic DNA contamination (DNase miss).
parse_read_distribution <- function(dir) {
  df <- load_mqc_module(dir, "multiqc_rseqc_read_distribution.txt")
  if (is.null(df)) return(NULL)
  cds  <- numeric_col(df, "cds_exons_tag_pct")
  utr5 <- numeric_col(df, "5_utr_exons_tag_pct")
  utr3 <- numeric_col(df, "3_utr_exons_tag_pct")
  exonic <- rowSums(cbind(cds, utr5, utr3), na.rm = TRUE)
  exonic[is.na(cds) & is.na(utr5) & is.na(utr3)] <- NA_real_
  introns <- numeric_col(df, "introns_tag_pct")
  inter_mat <- cbind(
    numeric_col(df, "tss_up_1kb_tag_pct"),
    numeric_col(df, "tss_up_5kb_tag_pct"),
    numeric_col(df, "tss_up_10kb_tag_pct"),
    numeric_col(df, "tes_down_1kb_tag_pct"),
    numeric_col(df, "tes_down_5kb_tag_pct"),
    numeric_col(df, "tes_down_10kb_tag_pct"),
    numeric_col(df, "other_intergenic_tag_pct")
  )
  intergenic <- rowSums(inter_mat, na.rm = TRUE)
  intergenic[rowSums(!is.na(inter_mat)) == 0] <- NA_real_
  tibble(
    sample        = df$sample,
    exonic_pct    = exonic,
    intronic_pct  = introns,
    intergenic_pct = intergenic
  )
}

# RSeQC gene-body coverage 5'->3' bias. Coverage is binned into 100 percentiles
# of transcript length; the ratio of mean coverage at the 5' flank to the 3'
# flank summarises RNA integrity. ~1.0 = even coverage; <0.5 = strong 3' bias
# (typical of degraded mRNA, since poly-A enrichment retains 3' ends while
# truncated transcripts lose 5' coverage).
parse_gene_body_coverage <- function(dir) {
  df <- load_mqc_module(dir, "multiqc_rseqc_gene_body_coverage.txt")
  if (is.null(df)) return(NULL)
  pct_cols <- names(df)[grepl("^(i_)?\\d+$", names(df))]
  if (length(pct_cols) < 50) return(NULL)
  idx <- suppressWarnings(as.integer(sub("^i_", "", pct_cols)))
  ord <- order(idx)
  pct_cols <- pct_cols[ord]; idx <- idx[ord]
  mat <- as.matrix(df[, pct_cols])
  five  <- rowMeans(mat[, idx >=  5 & idx <= 15, drop = FALSE], na.rm = TRUE)
  three <- rowMeans(mat[, idx >= 85 & idx <= 95, drop = FALSE], na.rm = TRUE)
  bias <- ifelse(is.finite(three) & three > 0, five / three, NA_real_)
  tibble(sample = df$sample, gene_body_5_3_bias = bias)
}

# dupRadar slope. Low slope => PCR-driven duplication (real complexity loss);
# high slope (~1) => duplication tracks expression (benign for high-expr genes).
parse_dupradar <- function(dir) {
  df <- load_mqc_module(dir, "multiqc_dupradar.txt")
  if (is.null(df)) return(NULL)
  slope_col <- grep("slope", names(df), ignore.case = TRUE, value = TRUE)
  if (length(slope_col) == 0) return(NULL)
  tibble(
    sample = df$sample,
    dupradar_slope = suppressWarnings(as.numeric(df[[slope_col[1]]]))
  )
}

mqc <- load_multiqc_general_stats(multiqc_dir)

if (!is.null(mqc)) {
  missing_in_mqc <- setdiff(samples$sample, mqc$sample)
  extra_in_mqc  <- setdiff(mqc$sample, samples$sample)
  if (length(missing_in_mqc) > 0) {
    message(sprintf("qc_summary: samples missing from MultiQC general_stats: %s",
                    paste(missing_in_mqc, collapse = ", ")))
  }
  if (length(extra_in_mqc) > 0) {
    message(sprintf("qc_summary: MultiQC general_stats has rows not in samples.tsv: %s",
                    paste(extra_in_mqc, collapse = ", ")))
  }

  mqc_tbl <- tibble(
    sample              = mqc$sample,
    total_reads         = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$total_reads))),
    mapping_rate        = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$mapping_rate))),
    uniquely_mapped_pct = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$uniquely_mapped_pct))),
    rrna_pct            = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$rrna_pct))),
    duplication         = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$duplication))),
    gc_content_pct      = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$gc_content_pct))),
    insert_size_median  = suppressWarnings(as.numeric(pick_metric(mqc, metric_aliases$insert_size_median)))
  )
  message(sprintf("qc_summary: loaded MultiQC general stats (%d samples)", nrow(mqc_tbl)))
} else {
  message("qc_summary: no MultiQC general_stats file found; QC metrics set to NA")
  mqc_tbl <- tibble(
    sample              = samples$sample,
    total_reads         = NA_real_,
    mapping_rate        = NA_real_,
    uniquely_mapped_pct = NA_real_,
    rrna_pct            = NA_real_,
    duplication         = NA_real_,
    gc_content_pct      = NA_real_,
    insert_size_median  = NA_real_
  )
}

# Per-tool MultiQC modules (optional — absence does not fail the rule).
extra_tbls <- list(
  strand    = parse_infer_experiment(multiqc_dir),
  readdist  = parse_read_distribution(multiqc_dir),
  dupradar  = parse_dupradar(multiqc_dir),
  genebody  = parse_gene_body_coverage(multiqc_dir)
)

summary_tbl <- samples %>%
  select(sample, condition, replicate, batch) %>%
  left_join(mqc_tbl, by = "sample")
for (nm in names(extra_tbls)) {
  tbl <- extra_tbls[[nm]]
  if (is.null(tbl)) next
  message(sprintf("qc_summary: loaded %s (%d rows)", nm, nrow(tbl)))
  summary_tbl <- left_join(summary_tbl, tbl, by = "sample")
}
summary_tbl <- left_join(summary_tbl, assigned_tbl, by = "sample")

for (col in REQUIRED_COLS) {
  if (!col %in% names(summary_tbl)) summary_tbl[[col]] <- NA
}
summary_tbl <- summary_tbl[, REQUIRED_COLS]

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write_tsv(summary_tbl, out_path)

message(sprintf("qc_summary: wrote %d rows x %d cols -> %s",
                nrow(summary_tbl), ncol(summary_tbl), out_path))
