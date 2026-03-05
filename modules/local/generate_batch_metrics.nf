#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GENERATEBATCHMETRICS {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(ped_file), path(pe_file), path(baf_file), path(rd_file), path(sr_file), path(median_file), path(clustered_depth_vcf), path(clustered_manta_vcf), path(clustered_wham_vcf), path(clustered_scramble_vcf)

    output:
    tuple val(batch_name), path("generate_batch_metrics_results"), emit: generate_batch_metrics_results
    tuple val(batch_name), path("generate_batch_metrics_results/exposed/metrics.tsv"), emit: metrics
    tuple val(batch_name), path("generate_batch_metrics_results/exposed/ploidy_table.tsv"), emit: ploidy_table
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.generatebatchmetrics_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.generate_batch_metrics ?: [:])

    def batch_id = batch_name?.toString()?.trim()
    if (!batch_id) {
        throw new IllegalArgumentException("GenerateBatchMetrics batch name is required")
    }
    def cohort_id = cohort?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("GenerateBatchMetrics cohort is required")
    }

    def dynamic = [
        "GenerateBatchMetrics.batch"       : batch_id,
        "GenerateBatchMetrics.ped_file"    : ped_file.toRealPath().toString(),
        "GenerateBatchMetrics.pe_file"     : pe_file.toRealPath().toString(),
        "GenerateBatchMetrics.baf_file"    : baf_file.toRealPath().toString(),
        "GenerateBatchMetrics.rd_file"     : rd_file.toRealPath().toString(),
        "GenerateBatchMetrics.sr_file"     : sr_file.toRealPath().toString(),
        "GenerateBatchMetrics.median_file" : median_file.toRealPath().toString(),
        "GenerateBatchMetrics.depth_vcf"   : clustered_depth_vcf.toRealPath().toString(),
        "GenerateBatchMetrics.manta_vcf"   : clustered_manta_vcf.toRealPath().toString(),
        "GenerateBatchMetrics.wham_vcf"    : clustered_wham_vcf.toRealPath().toString(),
        "GenerateBatchMetrics.scramble_vcf": clustered_scramble_vcf.toRealPath().toString()
    ]

    file("generate_batch_metrics_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV GENERATEBATCHMETRICS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out generate_batch_metrics_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file generate_batch_metrics_dynamic.json

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.generatebatchmetrics_wdl} \\
        -i generate_batch_metrics_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p generate_batch_metrics_results
    cp generate_batch_metrics_inputs.json generate_batch_metrics_results/
    find cromwell-executions/GenerateBatchMetrics/ -name "call-*" -type d -exec cp -r {} generate_batch_metrics_results/ \\;
    mkdir -p generate_batch_metrics_results/exposed

    metrics_file=\$(find generate_batch_metrics_results -type f -name "*.tsv" | grep -E -i "metrics" | head -n 1 || true)
    ploidy_file=\$(find generate_batch_metrics_results -type f -name "*.tsv" | grep -E -i "ploidy.*table|table.*ploidy" | head -n 1 || true)

    if [[ -z "\${metrics_file}" || -z "\${ploidy_file}" ]]; then
        echo "ERROR: GenerateBatchMetrics could not resolve metrics/ploidy_table outputs" >&2
        exit 1
    fi

    cp "\${metrics_file}" generate_batch_metrics_results/exposed/metrics.tsv
    cp "\${ploidy_file}" generate_batch_metrics_results/exposed/ploidy_table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p generate_batch_metrics_results/call-stub
    mkdir -p generate_batch_metrics_results/exposed
    touch generate_batch_metrics_results/call-stub/.stub
    touch generate_batch_metrics_results/exposed/metrics.tsv
    touch generate_batch_metrics_results/exposed/ploidy_table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
