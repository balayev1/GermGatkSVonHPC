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


args <- commandArgs(trailingOnly = TRUE)
metadata_path <- args[1]
nf_work_dir <- args[2]
num_samples  <- as.numeric(args[3])

## Set output directory
outdir <- file.path(nf_work_dir, "Sample_QC_out")
if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}

## Load input QC data
evid_qc_path <- system(paste0("find ", nf_work_dir, " -path '*/*/evidence_qc_results/call-MakeQcTable/execution/*.evidence_qc_table.tsv'"), intern = TRUE)
evid_qc_file <- read.table(evid_qc_path, sep = "\t", header = TRUE)[1:num_samples,]

### Load insertsize data
insertsize_paths <- system(paste0("find ", nf_work_dir, " -name '*.insert_size_metrics.txt'"), intern = TRUE)

## Load sex assignment files
evid_sex_path <- system(paste0("find ", nf_work_dir, " -path '*/*/evidence_qc_results/call-MakeQcTable/execution/ploidy_est/*'"), intern = TRUE)

# -------------------- Filters --------------------
excluded_samples <- c()
reason_excluded <- c()
## I. Filter by WGD, whole-genome dosage score = log-ratio of observed coverage versus expected coverage in specific, stable regions of the genome defined by wgd_mask file
### Remove samples with outlier dosage scores (based on 1st quartile - 7 * MAD and 3rd quartile + 7 * MAD)
### In general, distribution of WGD for PCR- samples is expected to be slightly lower than 0, and the distribution of WGD for PCR+ samples is 
### expected to be slightly greater than 0.

# summary(evid_qc_file$wgd_score)

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
    mad_autosomal_ploidy[[i]] <- mad(evid_qc_file[i, cn_cols])
    names(mad_autosomal_ploidy)[i] <- evid_qc_file$sample_id[i]
}

## Remove samples with autosomal ploidy MAD > 0.1
for (i in 1:length(mad_autosomal_ploidy)) {
    if (mad_autosomal_ploidy[[i]] > 0.1) {
        excluded_samples <- c(excluded_samples, names(mad_autosomal_ploidy)[i])
        reason_excluded <- c(reason_excluded, "autosomal_aneuploidy")
    }
}

### V. Sex Concordance
## Collect actual sex of samples
actual_sex_assign <- read.table(metadata_path, sep = "\t", header = FALSE)[, c(1,2)]
colnames(actual_sex_assign) <- c("sample_id", "sex")