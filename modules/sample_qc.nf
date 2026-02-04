#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SAMPLE_QC {
    tag "Sample_QC"
    
    input:
    path metadata_tsv
    val work_dir
    val num_samples

    output:
    path "Sample_QC_out/*", emit: sample_qc_results

    script:
    """
    Rscript ${baseDir}/R/Sample_QC.R ${metadata_tsv} "${work_dir}" ${num_samples}
    """
}