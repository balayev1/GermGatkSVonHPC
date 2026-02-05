#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SAMPLE_QC {
    tag "Sample_QC"
    
    container 'docker://rocker/tidyverse:latest'

    input:
    path metadata_tsv
    val num_samples
    path "insert_size_files/*"
    path "evid_qc_results"

    output:
    path "*.pdf", emit: sample_qc_plots
    path "*.tsv", emit: sample_qc_reports

    script:
    """
    Rscript ${baseDir}/R/Sample_QC.R \\
        ${metadata_tsv} \\
        "evid_qc_results" \\
        "insert_size_files" \\
        ${num_samples}
    """
}