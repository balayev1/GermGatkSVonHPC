#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_FILTERBATCHSAMPLES {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(depth_vcf), path(manta_vcf), path(wham_vcf), path(scramble_vcf), path(outlier_cutoff_table)

    output:
    tuple val(batch_name), path("filter_batch_samples_results"), emit: filter_batch_samples_results
    tuple val(batch_name), path("filter_batch_samples_results/exposed/filtered_depth.vcf*"), emit: filtered_depth_vcf
    tuple val(batch_name), path("filter_batch_samples_results/exposed/filtered_pesr.vcf*"), emit: filtered_pesr_vcf
    tuple val(batch_name), path("filter_batch_samples_results/exposed/outlier_samples_excluded.txt"), emit: outlier_samples_excluded_file
    tuple val(batch_name), path("filter_batch_samples_results/exposed/filtered_batch_samples.txt"), emit: filtered_batch_samples_file
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

    mkdir -p filter_batch_samples_results
    cp filter_batch_samples_inputs.json filter_batch_samples_results/
    find cromwell-executions/FilterBatchSamples/ -name "call-*" -type d -exec cp -r {} filter_batch_samples_results/ \\;

    mkdir -p filter_batch_samples_results/exposed

    depth_candidate=\$(find filter_batch_samples_results -type f \\( -name "*.vcf.gz" -o -name "*.vcf" \\) | grep -E -i "depth" | grep -E -i "filtered" | head -n 1 || true)
    pesr_candidate=\$(find filter_batch_samples_results -type f \\( -name "*.vcf.gz" -o -name "*.vcf" \\) | grep -E -i "pesr" | grep -E -i "filtered" | head -n 1 || true)
    if [[ -z "\${pesr_candidate}" ]]; then
        pesr_candidate=\$(find filter_batch_samples_results -type f \\( -name "*.vcf.gz" -o -name "*.vcf" \\) | grep -E -i "manta|wham|scramble" | grep -E -i "filtered" | head -n 1 || true)
    fi

    outlier_file=\$(find filter_batch_samples_results -type f \\( -name "*.txt" -o -name "*.tsv" \\) | grep -E -i "outlier.*excluded|excluded.*outlier" | head -n 1 || true)
    filtered_batch_file=\$(find filter_batch_samples_results -type f \\( -name "*.txt" -o -name "*.tsv" \\) | grep -E -i "filtered.*batch.*samples|batch.*samples.*filtered|remaining.*samples" | head -n 1 || true)

    if [[ -z "\${depth_candidate}" || -z "\${pesr_candidate}" || -z "\${outlier_file}" || -z "\${filtered_batch_file}" ]]; then
        echo "ERROR: FilterBatchSamples could not resolve one or more required outputs" >&2
        exit 1
    fi

    depth_suffix=".vcf"
    [[ "\${depth_candidate}" == *.vcf.gz ]] && depth_suffix=".vcf.gz"
    pesr_suffix=".vcf"
    [[ "\${pesr_candidate}" == *.vcf.gz ]] && pesr_suffix=".vcf.gz"

    cp "\${depth_candidate}" "filter_batch_samples_results/exposed/filtered_depth\${depth_suffix}"
    cp "\${pesr_candidate}" "filter_batch_samples_results/exposed/filtered_pesr\${pesr_suffix}"
    cp "\${outlier_file}" filter_batch_samples_results/exposed/outlier_samples_excluded.txt
    cp "\${filtered_batch_file}" filter_batch_samples_results/exposed/filtered_batch_samples.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p filter_batch_samples_results/call-stub
    mkdir -p filter_batch_samples_results/exposed
    touch filter_batch_samples_results/call-stub/.stub
    touch filter_batch_samples_results/exposed/filtered_depth.vcf.gz
    touch filter_batch_samples_results/exposed/filtered_pesr.vcf.gz
    touch filter_batch_samples_results/exposed/outlier_samples_excluded.txt
    touch filter_batch_samples_results/exposed/filtered_batch_samples.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
