#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GENERATEBATCHMETRICS {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), path(ped_file), path(pe_file), path(baf_file), path(rd_file), path(sr_file), path(median_file), path(clustered_depth_vcf), path(clustered_manta_vcf), path(clustered_wham_vcf), path(clustered_scramble_vcf), val(outlier_sample_ids)

    output:
    tuple val(batch_name), path("${batch_name}/${batch_name}.batch_metrics.aggregate_tests.metrics.tsv"), emit: metrics
    tuple val(batch_name), path("${batch_name}/GenerateBatchMetrics.${batch_name}.metrics.tsv"), emit: metrics_file_batchmetrics
    tuple val(batch_name), path("${batch_name}/${batch_name}.ploidy.tsv"), emit: ploidy_table
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.generatebatchmetrics_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.generate_batch_metrics ?: [:])

    def batch_id = batch_name?.toString()?.trim()
    if (!batch_id) {
        throw new IllegalArgumentException("GenerateBatchMetrics batch name is required")
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
    if (outlier_sample_ids != null && outlier_sample_ids.toString().trim()) {
        dynamic["GenerateBatchMetrics.outlier_sample_ids"] = file(outlier_sample_ids.toString()).toRealPath().toString()
    }

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

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.generatebatchmetrics_wdl} \\
        -i generate_batch_metrics_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p "${batch_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${batch_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/GenerateBatchMetrics -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected GenerateBatchMetrics output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${batch_id}.batch_metrics.aggregate_tests.metrics.tsv" 1
    copy_outputs "GenerateBatchMetrics.${batch_id}.metrics.tsv" 1
    copy_outputs "${batch_id}.ploidy.tsv" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${batch_name}
    touch ${batch_name}/${batch_name}.batch_metrics.aggregate_tests.metrics.tsv
    touch ${batch_name}/GenerateBatchMetrics.${batch_name}.metrics.tsv
    touch ${batch_name}/${batch_name}.ploidy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
