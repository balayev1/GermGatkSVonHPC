#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_FILTERBATCHSAMPLES {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(depth_vcf), path(manta_vcf), path(wham_vcf), path(scramble_vcf), path(outlier_cutoff_table)

    output:
    tuple val(batch_name), path("**/*.depth.outliers_removed.vcf*"), emit: filtered_depth_vcf
    tuple val(batch_name), path("**/*.filtered_pesr_merged.vcf*"), emit: filtered_pesr_vcf
    tuple val(batch_name), path("**/*.outliers.samples.list"), emit: outlier_samples_excluded_file
    tuple val(batch_name), path("**/*.outliers_excluded.samples.list"), emit: filtered_batch_samples_file
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.filterbatchsamples_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.filter_batch_samples ?: [:])

    def batch_id = batch_name?.toString()?.trim()
    if (!batch_id) {
        throw new IllegalArgumentException("FilterBatchSamples batch name is required")
    }
    def cohort_id = cohort?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("FilterBatchSamples cohort is required")
    }

    def dynamic = [
        "FilterBatchSamples.batch"              : batch_id,
        "FilterBatchSamples.depth_vcf"          : depth_vcf.toRealPath().toString(),
        "FilterBatchSamples.manta_vcf"          : manta_vcf.toRealPath().toString(),
        "FilterBatchSamples.wham_vcf"           : wham_vcf.toRealPath().toString(),
        "FilterBatchSamples.scramble_vcf"       : scramble_vcf.toRealPath().toString(),
        "FilterBatchSamples.outlier_cutoff_table": outlier_cutoff_table.toRealPath().toString()
    ]
    file("filter_batch_samples_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV FILTERBATCHSAMPLES] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out filter_batch_samples_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file filter_batch_samples_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.filterbatchsamples_wdl} \\
        -i filter_batch_samples_inputs.json \\
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
    touch ${batch_name}.depth.outliers_removed.vcf.gz
    touch ${batch_name}.filtered_pesr_merged.vcf.gz
    touch ${batch_name}.outliers.samples.list
    touch ${batch_name}.outliers_excluded.samples.list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
