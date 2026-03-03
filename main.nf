#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GermGatkSVonHPC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/balayev1/GermGatkSVonHPC

----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.reference_fasta = WorkflowMain.getGenomeAttribute(params, 'reference_fasta')
params.reference_index = WorkflowMain.getGenomeAttribute(params, 'reference_index')
params.reference_dict = WorkflowMain.getGenomeAttribute(params, 'reference_dict')
params.reference_version = WorkflowMain.getGenomeAttribute(params, 'reference_version')
params.primary_contigs_list = WorkflowMain.getGenomeAttribute(params, 'primary_contigs_list')
params.primary_contigs_fai = WorkflowMain.getGenomeAttribute(params, 'primary_contigs_fai')
params.preprocessed_intervals = WorkflowMain.getGenomeAttribute(params, 'preprocessed_intervals')
params.manta_region_bed = WorkflowMain.getGenomeAttribute(params, 'manta_region_bed')
params.manta_region_bed_index = WorkflowMain.getGenomeAttribute(params, 'manta_region_bed_index')
params.mei_bed = WorkflowMain.getGenomeAttribute(params, 'mei_bed')
params.sd_locs_vcf = WorkflowMain.getGenomeAttribute(params, 'sd_locs_vcf')
params.sd_locs_vcf_index = WorkflowMain.getGenomeAttribute(params, 'sd_locs_vcf_index')
params.wham_include_list_bed_file = WorkflowMain.getGenomeAttribute(params, 'wham_include_list_bed_file')
params.genome_file = WorkflowMain.getGenomeAttribute(params, 'genome_file')
params.wgd_scoring_mask = WorkflowMain.getGenomeAttribute(params, 'wgd_scoring_mask')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATKSVCALLER } from './workflows/gatksvcaller.nf'

//
// WORKFLOW: Run main GERMGATKSVONHPC analysis pipeline
//
workflow NFCORE_GATKSV {
    GATKSVCALLER ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {
    NFCORE_GATKSV ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
