#!/usr/bin/env R

# Sample_QC.R script to filter poor quality samples based on  QC metrics for downstream germline SV analysis
# /scratch.global/balay011/GermlineSV_outs/EvidenceQC_out/results/cromwell-executions/EvidenceQC/a54207f1-a499-4434-b02f-162fd10b8f77/call-MakeQcTable/execution/CHORDBAI.evidence_qc_table.tsv:
# - includes ploidy per chromosome, sex assignment, outlier status, median coverage, whole-genome dosage score, mean_insert_size, 
# manta_high_outlier, wham_high_outlier, scramble_high_outlier, dragen_high_outlier, overall_high_outlier, manta_low_outlier, melt_low_outlier, wham_low_outlier, scramble_low_outlier, dragen_low_outlier, overall_low_outlier

# mean insert size must be taken from multiqc output
# put side-by-side sex assignment from plink and this
# sample_sex_assignments.txt.gz contains probability of mosaicism

## Load required libraries
# suppressPackageStartupMessages({})

## Set output directory
outdir <- "/scratch.global/balay011/GermlineSV_outs/Sample_QC_out/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

num_samples <- 80 

## Load input QC data
evid_qc_path <- "/scratch.global/balay011/GermlineSV_outs/EvidenceQC_out/results/cromwell-executions/EvidenceQC/a54207f1-a499-4434-b02f-162fd10b8f77/call-MakeQcTable/execution/CHORDBAI.evidence_qc_table.tsv"
evid_qc_file <- read.table(evid_qc_path, sep = "\t", header = TRUE)[1:num_samples,]

### Load insertsize data
insertsize_path <- "/users/1/balay011/helper_scripts/multiqc_qualimap_bamqc_genome_results.txt"
insertsize_file <- read.table(insertsize_path, sep = "\t", header = TRUE)

## Load sex assignment files
plink.sex.path <- ""
evid.sex.path <- ""

# -------------------- Filters --------------------
excluded_samples <- c()
reason_excluded <- c()
## I. Filter by WGD, whole-genome dosage score = log-ratio of observed coverage versus expected coverage in specific, stable regions of the genome defined by wgd_mask file
### Remove samples with outlier dosage scores (based on 1st quartile - 8 * MAD and 3rd quartile + 8 * MAD)
### In general, distribution of WGD for PCR- samples is expected to be slightly lower than 0, and the distribution of WGD for PCR+ samples is 
### expected to be slightly greater than 0.

# summary(evid_qc_file$wgd_score)

first_quartile <- quantile(evid_qc_file$wgd_score, 0.25)
third_quartile <- quantile(evid_qc_file$wgd_score, 0.75)
mad <- mad(evid_qc_file$wgd_score)

if (min(evid_qc_file$wgd_score) < (first_quartile - 8 * mad) | max(evid_qc_file$wgd_score) > (third_quartile + 8 * mad)) {
    low_outliers <- evid_qc_file$sample_id[evid_qc_file$wgd_score < (first_quartile - 8 * mad)]
    excluded_samples <- c(excluded_samples, as.character(low_outliers))
    reason_excluded <- c(reason_excluded, rep("low_WGD_score", length(low_outliers)))

    high_outliers <- evid_qc_file$sample_id[evid_qc_file$wgd_score > (third_quartile + 8 * mad)]
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
## Remove samples with outlier insert size (based on 1st quartile - 8 * MAD and 3rd quartile + 8 * MAD)
if (length(insertsize_path) > 0 && file.exists(insertsize_path)) {
    insertsize_file <- read.table(insertsize_path, sep = "\t", header = TRUE)
    evid_qc_file$median_insert_size <- insertsize_file$median_insert_size[match(evid_qc_file$sample_id, insertsize_file$sample_id)]
}
