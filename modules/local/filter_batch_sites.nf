#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_FILTERBATCHSITES {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(depth_vcf), path(manta_vcf), path(wham_vcf), path(scramble_vcf), path(evidence_metrics)

    output:
    tuple val(batch_name), path("${batch_name}/${batch_name}.depth.with_evidence.vcf.gz"), emit: sites_filtered_depth_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.manta.with_evidence.vcf.gz"), emit: sites_filtered_manta_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.scramble.with_evidence.vcf.gz"), emit: sites_filtered_scramble_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.wham.with_evidence.vcf.gz"), emit: sites_filtered_wham_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.cutoffs"), emit: cutoffs
    tuple val(batch_name), path("${batch_name}/${batch_name}.scores"), emit: scores
    tuple val(batch_name), path("${batch_name}/${batch_name}.RF_intermediate_files.tar.gz"), emit: RF_intermediate_files
    tuple val(batch_name), path("${batch_name}/*.with_evidence.svcounts.txt"), emit: sv_counts
    tuple val(batch_name), path("${batch_name}/*.with_evidence.all_SVTYPEs.counts_per_sample.png"), emit: sv_count_plots
    tuple val(batch_name), path("${batch_name}/${batch_name}.outliers_preview.samples.txt"), emit: sites_filtered_outlier_samples_preview
    tuple val(batch_name), path("${batch_name}/${batch_name}.outliers_preview_with_reason.samples.tsv"), emit: sites_filtered_outlier_samples_with_reason
    tuple val(batch_name), path("${batch_name}/num_outliers.txt"), emit: num_outlier_samples
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

    mkdir -p "${batch_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${batch_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/FilterBatchSites -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected FilterBatchSites output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${batch_id}.depth.with_evidence.vcf.gz" 1
    copy_outputs "${batch_id}.manta.with_evidence.vcf.gz" 1
    copy_outputs "${batch_id}.scramble.with_evidence.vcf.gz" 1
    copy_outputs "${batch_id}.wham.with_evidence.vcf.gz" 1
    copy_outputs "${batch_id}.cutoffs" 1
    copy_outputs "${batch_id}.scores" 1
    copy_outputs "${batch_id}.RF_intermediate_files.tar.gz" 1
    copy_outputs "*.with_evidence.svcounts.txt" 1
    copy_outputs "*.with_evidence.all_SVTYPEs.counts_per_sample.png" 1
    copy_outputs "${batch_id}.outliers_preview.samples.txt" 1
    copy_outputs "${batch_id}.outliers_preview_with_reason.samples.tsv" 1
    copy_outputs "num_outliers.txt" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${batch_name}
    touch ${batch_name}/${batch_name}.depth.with_evidence.vcf.gz
    touch ${batch_name}/${batch_name}.manta.with_evidence.vcf.gz
    touch ${batch_name}/${batch_name}.scramble.with_evidence.vcf.gz
    touch ${batch_name}/${batch_name}.wham.with_evidence.vcf.gz
    touch ${batch_name}/${batch_name}.cutoffs
    touch ${batch_name}/${batch_name}.scores
    touch ${batch_name}/${batch_name}.RF_intermediate_files.tar.gz
    touch ${batch_name}/${batch_name}.with_evidence.svcounts.txt
    touch ${batch_name}/${batch_name}.with_evidence.all_SVTYPEs.counts_per_sample.png
    touch ${batch_name}/${batch_name}.outliers_preview.samples.txt
    touch ${batch_name}/${batch_name}.outliers_preview_with_reason.samples.tsv
    touch ${batch_name}/num_outliers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
