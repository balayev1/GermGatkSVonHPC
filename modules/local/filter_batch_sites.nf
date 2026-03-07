#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_FILTERBATCHSITES {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(depth_vcf), path(manta_vcf), path(wham_vcf), path(scramble_vcf), path(evidence_metrics)

    output:
    tuple val(batch_name), path("**/*.depth.with_evidence.vcf.gz"), emit: sites_filtered_depth_vcf
    tuple val(batch_name), path("**/*.manta.with_evidence.vcf.gz"), emit: sites_filtered_manta_vcf
    tuple val(batch_name), path("**/*.scramble.with_evidence.vcf.gz"), emit: sites_filtered_scramble_vcf
    tuple val(batch_name), path("**/*.wham.with_evidence.vcf.gz"), emit: sites_filtered_wham_vcf
    tuple val(batch_name), path("**/${batch_name}.cutoffs"), emit: cutoffs
    tuple val(batch_name), path("**/${batch_name}.scores"), emit: scores
    tuple val(batch_name), path("**/${batch_name}.RF_intermediate_files.tar.gz"), emit: rf_intermediates
    tuple val(batch_name), path("**/*.svcounts.txt"), emit: sv_counts
    tuple val(batch_name), path("**/*.all_SVTYPEs.counts_per_sample.png"), emit: sv_count_plots
    tuple val(batch_name), path("**/*.outliers_preview.samples.txt"), emit: outlier_samples_preview
    tuple val(batch_name), path("**/*.outliers_preview_with_reason.samples.tsv"), emit: outlier_samples_with_reason
    tuple val(batch_name), path("num_outliers.txt"), emit: num_outlier_samples
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.filterbatchsites_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.filter_batch_sites ?: [:])

    def batch_id = batch_name?.toString()?.trim()
    if (!batch_id) {
        throw new IllegalArgumentException("FilterBatchSites batch name is required")
    }
    def cohort_id = cohort?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("FilterBatchSites cohort is required")
    }

    def dynamic = [
        "FilterBatchSites.batch"           : batch_id,
        "FilterBatchSites.depth_vcf"       : depth_vcf.toRealPath().toString(),
        "FilterBatchSites.manta_vcf"       : manta_vcf.toRealPath().toString(),
        "FilterBatchSites.wham_vcf"        : wham_vcf.toRealPath().toString(),
        "FilterBatchSites.scramble_vcf"    : scramble_vcf.toRealPath().toString(),
        "FilterBatchSites.evidence_metrics": evidence_metrics.toRealPath().toString()
    ]
    file("filter_batch_sites_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV FILTERBATCHSITES] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out filter_batch_sites_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file filter_batch_sites_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.filterbatchsites_wdl} \\
        -i filter_batch_sites_inputs.json \\
        -p ${params.deps_zip}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p call-stub
    touch call-stub/.stub
    touch ${batch_name}.depth.with_evidence.vcf.gz
    touch ${batch_name}.manta.with_evidence.vcf.gz
    touch ${batch_name}.scramble.with_evidence.vcf.gz
    touch ${batch_name}.wham.with_evidence.vcf.gz
    touch ${batch_name}.cutoffs
    touch ${batch_name}.scores
    touch ${batch_name}.RF_intermediate_files.tar.gz
    touch ${batch_name}.svcounts.txt
    touch ${batch_name}.all_SVTYPEs.counts_per_sample.png
    touch ${batch_name}.outliers_preview.samples.txt
    touch ${batch_name}.outliers_preview_with_reason.samples.tsv
    touch num_outliers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
