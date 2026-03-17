#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_FILTERBATCHSAMPLES {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(depth_vcf), path(manta_vcf), path(wham_vcf), path(scramble_vcf), val(outlier_cutoff_table)

    output:
    tuple val(batch_name), path("${batch_name}/${batch_name}.depth.outliers_removed.vcf.gz"), emit: outlier_filtered_depth_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.depth.outliers_removed.vcf.gz.tbi"), emit: outlier_filtered_depth_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.manta.outliers_removed.vcf.gz"), emit: outlier_filtered_manta_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.manta.outliers_removed.vcf.gz.tbi"), emit: outlier_filtered_manta_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.scramble.outliers_removed.vcf.gz"), emit: outlier_filtered_scramble_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.scramble.outliers_removed.vcf.gz.tbi"), emit: outlier_filtered_scramble_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.wham.outliers_removed.vcf.gz"), emit: outlier_filtered_wham_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.wham.outliers_removed.vcf.gz.tbi"), emit: outlier_filtered_wham_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.filtered_pesr_merged.vcf.gz"), emit: outlier_filtered_pesr_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.filtered_pesr_merged.vcf.gz.tbi"), emit: outlier_filtered_pesr_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.outliers_excluded.samples.list"), emit: filtered_batch_samples_file
    tuple val(batch_name), path("${batch_name}/${batch_name}.outliers.samples.list"), emit: outlier_samples_excluded_file
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
        "FilterBatchSamples.scramble_vcf"       : scramble_vcf.toRealPath().toString()
    ]
    if (outlier_cutoff_table != null && outlier_cutoff_table.toString().trim()) {
        dynamic["FilterBatchSamples.outlier_cutoff_table"] = file(outlier_cutoff_table.toString()).toRealPath().toString()
    }
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

    mkdir -p "${batch_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${batch_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/FilterBatchSamples -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected FilterBatchSamples output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${batch_id}.depth.outliers_removed.vcf.gz" 1
    copy_outputs "${batch_id}.depth.outliers_removed.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.manta.outliers_removed.vcf.gz" 1
    copy_outputs "${batch_id}.manta.outliers_removed.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.scramble.outliers_removed.vcf.gz" 1
    copy_outputs "${batch_id}.scramble.outliers_removed.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.wham.outliers_removed.vcf.gz" 1
    copy_outputs "${batch_id}.wham.outliers_removed.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.filtered_pesr_merged.vcf.gz" 1
    copy_outputs "${batch_id}.filtered_pesr_merged.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.outliers_excluded.samples.list" 1
    copy_outputs "${batch_id}.outliers.samples.list" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${batch_name}
    touch ${batch_name}/${batch_name}.depth.outliers_removed.vcf.gz
    touch ${batch_name}/${batch_name}.depth.outliers_removed.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.manta.outliers_removed.vcf.gz
    touch ${batch_name}/${batch_name}.manta.outliers_removed.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.scramble.outliers_removed.vcf.gz
    touch ${batch_name}/${batch_name}.scramble.outliers_removed.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.wham.outliers_removed.vcf.gz
    touch ${batch_name}/${batch_name}.wham.outliers_removed.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.filtered_pesr_merged.vcf.gz
    touch ${batch_name}/${batch_name}.filtered_pesr_merged.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.outliers_excluded.samples.list
    touch ${batch_name}/${batch_name}.outliers.samples.list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
