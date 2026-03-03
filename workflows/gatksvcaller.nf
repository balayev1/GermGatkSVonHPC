/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
if (!params.input) {
    exit 1, 'Input samplesheet not specified!'
}
if (!params.sd_locs_vcf) {
    exit 1, "Genome resource 'sd_locs_vcf' is not configured for genome '${params.genome}'."
}
if (!params.sd_locs_vcf_index) {
    exit 1, "Genome resource 'sd_locs_vcf_index' is not configured for genome '${params.genome}'."
}
if (!params.reference_dict) {
    exit 1, "Genome resource 'reference_dict' is not configured for genome '${params.genome}'."
}

ch_input = file(params.input)
ch_sd_locs_vcf = Channel.value(file(params.sd_locs_vcf, checkIfExists: true))
ch_sd_locs_vcf_index = Channel.value(file(params.sd_locs_vcf_index, checkIfExists: true))
ch_reference_dict = Channel.value(file(params.reference_dict, checkIfExists: true))
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK              } from '../subworkflows/local/input_check'
include { SAMPLE_PROCESSING        } from '../subworkflows/local/sample_processing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow GATKSVCALLER {
    main:
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_sample = INPUT_CHECK.out.reads

    //
    // SUBWORKFLOW: Collect per-sample SV evidence, perform sample QC and create batches if necessary
    //
    SAMPLE_PROCESSING (
        ch_sample,
        ch_sd_locs_vcf,
        ch_sd_locs_vcf_index,
        ch_reference_dict
    )
    ch_versions = ch_versions.mix(SAMPLE_PROCESSING.out.versions)
}
