#!/usr/bin/env R

# Sample_QC.R script filters poor quality samples based on  QC metrics for downstream germline SV analysis
# Input: sample metadata, nextflow working directory, number of samples to process
# Output:
# Algorithm: filters outlier samples based on:
# I) Extremely high/low whole-genome dosage score (<1st quartile - 7*MAD/>3rd quartile + 7*MAD)
# II) Extremely high/low median coverage (<15X or >60X)
# III) Extremely high/low median insert size (<1st quartile - 7*MAD/>3rd quartile + 7*MAD)
# IV) Autosomal aneuploidy (MAD of autosomal ploidy > 0.1)
# V) Sex discordance (discordance between actual and assigned sex)
# VI) Ambiguous sex assignment (probability of mosaicism > 0.5)
# VII) Overall outlier status (overall_high_outlier = 1 or overall_low_outlier = 1)
# Insert size is derived from Picard CollectInsertSizeMetrics output
# All other metrics are derived from Evidence_QC output

## Load required libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
metadata_path <- args[1]
nf_work_dir <- args[2]
num_samples  <- as.numeric(args[3])

## Set output directory
outdir <- "Sample_QC_out" 
if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}

## Load input QC data
# If you change the input to 'path evid_results', you can just do:
evid_qc_path <- list.files(path = "evid_qc_results", pattern = "evidence_qc_table.tsv", recursive = TRUE, full.names = TRUE)
evid_qc_file <- read.table(evid_qc_path, sep = "\t", header = TRUE)[1:num_samples,]

### Load insertsize data
insertsize_paths <- list.files(path = "insert_size_files", 
                               pattern = "\\.insert_size_metrics\\.txt$", 
                               full.names = TRUE)

## Load sex assignment files
evid_sex_path <- list.files(path = "evid_qc_results", pattern = "sample_sex_assignments.txt.gz", recursive = TRUE, full.names = TRUE)
evid_sex_file <- read.table(gzfile(evid_sex_path), sep = "\t", header = TRUE)


# -------------------- Filters --------------------
excluded_samples <- c()
reason_excluded <- c()
## I. Filter by WGD, whole-genome dosage score = log-ratio of observed coverage versus expected coverage in specific, stable regions of the genome defined by wgd_mask file
### Remove samples with outlier dosage scores (based on 1st quartile - 7 * MAD and 3rd quartile + 7 * MAD)
### In general, distribution of WGD for PCR- samples is expected to be slightly lower than 0, and the distribution of WGD for PCR+ samples is 
### expected to be slightly greater than 0.

first_quartile <- quantile(evid_qc_file$wgd_score, 0.25)
third_quartile <- quantile(evid_qc_file$wgd_score, 0.75)
mad_wgd <- mad(evid_qc_file$wgd_score)

if (min(evid_qc_file$wgd_score) < (first_quartile - 7 * mad_wgd) | max(evid_qc_file$wgd_score) > (third_quartile + 7 * mad_wgd)) {
    low_outliers <- evid_qc_file$sample_id[evid_qc_file$wgd_score < (first_quartile - 7 * mad_wgd)]
    excluded_samples <- c(excluded_samples, as.character(low_outliers))
    reason_excluded <- c(reason_excluded, rep("low_WGD_score", length(low_outliers)))

    high_outliers <- evid_qc_file$sample_id[evid_qc_file$wgd_score > (third_quartile + 7 * mad_wgd)]
    excluded_samples <- c(excluded_samples, as.character(high_outliers))
    reason_excluded <- c(reason_excluded, rep("high_WGD_score", length(high_outliers)))
}

## II. Median Coverage 
### Remove samples with median coverage < 15X and > 60X
low_cov_samples <- evid_qc_file$sample_id[evid_qc_file$median_coverage < 15]
excluded_samples <- c(excluded_samples, as.character(low_cov_samples))
reason_excluded <- c(reason_excluded, rep("low_median_coverage", length(low_cov_samples)))

high_cov_samples <- evid_qc_file$sample_id[evid_qc_file$median_coverage > 60]
excluded_samples <- c(excluded_samples, as.character(high_cov_samples))
reason_excluded <- c(reason_excluded, rep("high_median_coverage", length(high_cov_samples)))

### III. Median Insert Size
### Collect median insert size for each sample
for (insertsize_path in insertsize_paths) {
    insertsize_file <- read.table(insertsize_path, sep = "\t", header = TRUE, skip = 6, nrows = 1, comment.char = "")
    # Use basename to get the sample ID from the filename
    sample_id <- gsub(".insert_size_metrics.txt", "", basename(insertsize_path))
    evid_qc_file$median_insert_size[evid_qc_file$sample_id == sample_id] <- insertsize_file$MEDIAN_INSERT_SIZE
}

## Remove samples with outlier insert size (based on 1st quartile - 7 * MAD and 3rd quartile + 7 * MAD)
first_quartile <- quantile(evid_qc_file$median_insert_size, 0.25)
third_quartile <- quantile(evid_qc_file$median_insert_size, 0.75)
mad_insert_size <- mad(evid_qc_file$median_insert_size)

if (min(evid_qc_file$median_insert_size) < (first_quartile - 7 * mad_insert_size) | max(evid_qc_file$median_insert_size) > (third_quartile + 7 * mad_insert_size)) {
    low_outliers <- evid_qc_file$sample_id[evid_qc_file$median_insert_size < (first_quartile - 7 * mad_insert_size)]
    excluded_samples <- c(excluded_samples, as.character(low_outliers))
    reason_excluded <- c(reason_excluded, rep("low_median_insert_size", length(low_outliers)))

    high_outliers <- evid_qc_file$sample_id[evid_qc_file$median_insert_size > (third_quartile + 7 * mad_insert_size)]
    excluded_samples <- c(excluded_samples, as.character(high_outliers))
    reason_excluded <- c(reason_excluded, rep("high_median_insert_size", length(high_outliers)))
}

### IV. Autosomal Ploidy
## Calculate median absolute deviation (MAD) of autosomal ploidy (chromosomes 1-22) for each sample
mad_autosomal_ploidy <- c()
for (i in 1:nrow(evid_qc_file)) {
    cn_cols <- grep("^chr[1-9][0-9]*_CopyNumber", colnames(evid_qc_file))
    mad_autosomal_ploidy[[i]] <- mad(as.numeric(evid_qc_file[i, cn_cols]))
    names(mad_autosomal_ploidy)[i] <- evid_qc_file$sample_id[i]
}

## Remove samples with autosomal ploidy MAD > 0.1
for (i in 1:length(mad_autosomal_ploidy)) {
    if (mad_autosomal_ploidy[[i]] > 0.1) {
        excluded_samples <- c(excluded_samples, names(mad_autosomal_ploidy)[i])
        reason_excluded <- c(reason_excluded, "autosomal_aneuploidy")
    }
}

### V. Sex Discordance
## Collect actual sex of samples
actual_sex_assign <- read.table(metadata_path, sep = "\t", header = FALSE)[, c(1,2)]
colnames(actual_sex_assign) <- c("sample_id", "Assignment")
actual_sex_assign$actual_sex <- ifelse(actual_sex_assign$Assignment %in% c("XY"), "MALE", "FEMALE")
actual_sex_assign$inferred_sex <- evid_sex_file$Assignment[match(actual_sex_assign$sample_id, evid_sex_file$sample_id)]

## Remove samples with discordance between actual and inferred sex
discordant_samples <- actual_sex_assign$sample_id[actual_sex_assign$actual_sex != actual_sex_assign$inferred_sex]
excluded_samples <- c(excluded_samples, as.character(discordant_samples))
reason_excluded <- c(reason_excluded, rep("sex_discordance", length(discordant_samples)))

### VI. Overall Outlier Status
## Remove samples with overall_high_outlier = 1 or overall_low_outlier = 1
overall_high_outliers <- evid_qc_file$sample_id[evid_qc_file$overall_high_outlier == 1]
excluded_samples <- c(excluded_samples, as.character(overall_high_outliers))
reason_excluded <- c(reason_excluded, rep("overall_high_outlier", length(overall_high_outliers)))

overall_low_outliers <- evid_qc_file$sample_id[evid_qc_file$overall_low_outlier == 1]
excluded_samples <- c(excluded_samples, as.character(overall_low_outliers))
reason_excluded <- c(reason_excluded, rep("overall_low_outlier", length(overall_low_outliers)))



# -------------------- Plots --------------------
qc_colors <- c("PASS" = "#377eb8", "FAIL" = "#e41a1c")

### I. WGD Score Distribution
wgd_fail_ids <- excluded_samples[reason_excluded %in% c("low_WGD_score", "high_WGD_score")]
evid_qc_file$wgd_qc_status <- factor(ifelse(evid_qc_file$sample_id %in% wgd_fail_ids, "FAIL", "PASS"), levels = c("PASS", "FAIL"))

wgd_low <- quantile(evid_qc_file$wgd_score, 0.25) - (7 * mad(evid_qc_file$wgd_score))
wgd_high <- quantile(evid_qc_file$wgd_score, 0.75) + (7 * mad(evid_qc_file$wgd_score))

p1 <- ggplot(evid_qc_file, aes(x = wgd_score, fill = wgd_qc_status)) +
        geom_histogram(bins = 30, color = "black", alpha = 0.7, show.legend = TRUE) +
        geom_vline(xintercept = c(wgd_low, wgd_high), linetype = "dashed", color = "red", linewidth=1) +
        scale_fill_manual(values = qc_colors, drop = FALSE) +
        labs(x = "WGD Score", y = "Count", fill = "WGD QC Status") +
        theme_minimal()
ggsave(filename = file.path(outdir, "WGD_Score_Distribution.pdf"), plot = p1, width = 10, height = 8)

### II. Median Coverage Distribution
cov_fail_ids <- excluded_samples[reason_excluded %in% c("low_median_coverage", "high_median_coverage")]
evid_qc_file$cov_qc_status <- factor(ifelse(evid_qc_file$sample_id %in% cov_fail_ids, "FAIL", "PASS"), levels = c("PASS", "FAIL"))

p2 <- ggplot(evid_qc_file, aes(x = median_coverage, fill = cov_qc_status)) +
        geom_histogram(bins = 30, color = "black", alpha = 0.7, show.legend = TRUE) +
        geom_vline(xintercept = c(15, 60), linetype = "dashed", color = "red", linewidth=1) +
        scale_fill_manual(values = qc_colors, drop = FALSE) +
        labs(x = "Median Coverage", y = "Count", fill = "Coverage QC Status") +
        theme_minimal()
ggsave(filename = file.path(outdir, "Median_Coverage_Distribution.pdf"), plot = p2, width = 10, height = 8)

### III. Median Insert Size Distribution
insert_fail_ids <- excluded_samples[reason_excluded %in% c("low_median_insert_size", "high_median_insert_size")]
evid_qc_file$insert_qc_status <- factor(ifelse(evid_qc_file$sample_id %in% insert_fail_ids, "FAIL", "PASS"), levels = c("PASS", "FAIL"))

insert_low <- quantile(evid_qc_file$median_insert_size, 0.25) - (7 * mad(evid_qc_file$median_insert_size))
insert_high <- quantile(evid_qc_file$median_insert_size, 0.75) + (7 * mad(evid_qc_file$median_insert_size))

p3 <- ggplot(evid_qc_file, aes(x = median_insert_size, fill = insert_qc_status)) +
        geom_histogram(bins = 30, color = "black", alpha = 0.7, show.legend = TRUE) +
        geom_vline(xintercept = c(insert_low, insert_high), linetype = "dashed", color = "red", linewidth=1) +
        scale_fill_manual(values = qc_colors, drop = FALSE) +
        labs(x = "Median Insert Size", y = "Count", fill = "Insert Size QC Status") +
        theme_minimal()
ggsave(filename = file.path(outdir, "Median_Insert_Size_Distribution.pdf"), plot = p3, width = 10, height = 8)

### IV. Autosomal Ploidy MAD Distribution
autos_ploidy_df <- data.frame(sample_id = names(mad_autosomal_ploidy),
                             mad_autosomal_ploidy = unlist(mad_autosomal_ploidy))

autos_ploidy_fail_ids <- excluded_samples[reason_excluded %in% c("autosomal_aneuploidy")]
autos_ploidy_df$autos_ploidy_qc_status <- factor(ifelse(autos_ploidy_df$sample_id %in% autos_ploidy_fail_ids, "FAIL", "PASS"), 
levels = c("PASS", "FAIL"))

p4 <- ggplot(autos_ploidy_df, aes(x = mad_autosomal_ploidy, fill = autos_ploidy_qc_status)) +
        geom_histogram(bins = 30, color = "black", alpha = 0.7, show.legend = TRUE) +
        geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", linewidth=1) +
        scale_fill_manual(values = qc_colors, drop = FALSE) +
        labs(x = "Autosomal Ploidy MAD", y = "Count", fill = "Autosomal Ploidy QC Status") +
        theme_minimal()
ggsave(filename = file.path(outdir, "Autosomal_Ploidy_MAD_Distribution.pdf"), plot = p4, width = 10, height = 8)

### V. Sex Discordance Plot
evid_qc_file$actual_sex <- factor(actual_sex_assign$actual_sex[match(evid_qc_file$sample_id, actual_sex_assign$sample_id)],
    levels = c("MALE", "FEMALE"))
evid_qc_file$inferred_sex <- factor(actual_sex_assign$inferred_sex[match(evid_qc_file$sample_id, actual_sex_assign$sample_id)])
evid_qc_file$sex_discordance <- factor(ifelse(evid_qc_file$sample_id %in% discordant_samples, "DISCORDANT", "CONCORDANT"), 
    levels = c("CONCORDANT", "DISCORDANT"))

p5 <- ggplot(evid_qc_file, aes(x = chrX_CopyNumber, y = chrY_CopyNumber)) +
    geom_point(aes(color = actual_sex, shape = sex_discordance), size = 3, alpha = 0.8, show.legend = TRUE) +
    labs(x = "ChrX Copy Number", y = "ChrY Copy Number",
            color = "Sex Known", shape = "Sex Discordance") +
    scale_color_manual(
        values = c("MALE" = "#1f78b4", "FEMALE" = "#e31a1c"),
        drop = FALSE) +
    scale_shape_manual(
        values = c("CONCORDANT" = 16, "DISCORDANT" = 17),
        drop = FALSE) +
    theme_bw() +
    theme(legend.position = "right")
ggsave(filename = file.path(outdir, "Sex_Discordance.pdf"), plot = p5, width = 10, height = 8)

### VI. Structural Variant Counts by Caller
evid_qc_long <- evid_qc_file %>%
  pivot_longer(
    cols = starts_with(c("Manta_", "Wham_", "Scramble_", "Melt_"), ignore.case = FALSE), 
    names_to = "variable", 
    values_to = "count"
  ) %>%
  # Split the column names into Caller, Type, and Chrom
  tidyr::separate(variable, into = c("Caller", "SV_Type", "Chrom"), sep = "_") %>%
  # Ensure chromosomes are ordered
  mutate(Chrom = factor(Chrom, levels = paste0("chr", c(1:22, "X")))) %>%
  filter(!is.na(Chrom))

callers <- unique(evid_qc_long$Caller)

for (cmd_caller in callers) {
  
  # Filter for the specific caller
  plot_data <- evid_qc_long %>% filter(Caller == cmd_caller)
  
  p <- ggplot(plot_data, aes(x = SV_Type, y = count, fill = SV_Type)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8, aes(color = SV_Type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
    facet_wrap(~Chrom, scales = "free_y") + 
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90"),
      panel.border = element_rect(color = "grey", fill = NA)
    ) +
    labs(x = "SV Type",
        y = "SV Count")

  ggsave(filename = file.path(outdir, paste0("SV_Counts_by_", cmd_caller, ".pdf")), plot = p, width = 12, height = 8)
}


# -------------------- Outputs --------------------
if (length(excluded_samples) > 0) {
    raw_df <- data.frame(sample_id = excluded_samples, reason = reason_excluded, stringsAsFactors = FALSE)
    final_df <- aggregate(reason ~ sample_id, data = raw_df, 
                               FUN = function(x) paste(unique(x), collapse = ";"))
    write.table(final_df, file.path(outdir, "Excluded_Samples_Report.tsv"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
} else {
    final_df <- data.frame(sample_id = character(), reason = character()) 
    write.table(final_df, file.path(outdir, "Excluded_Samples_Report.tsv"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
}