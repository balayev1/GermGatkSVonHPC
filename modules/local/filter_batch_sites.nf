#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_FILTERBATCHSITES {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(depth_vcf), path(manta_vcf), path(wham_vcf), path(scramble_vcf), path(evidence_metrics)

    output:
    tuple val(batch_name), path("filter_batch_sites_results"), emit: filter_batch_sites_results
    tuple val(batch_name), path("filter_batch_sites_results/exposed/cutoffs.tsv"), emit: cutoffs
    tuple val(batch_name), path("filter_batch_sites_results/exposed/sv_counts/*.tsv"), emit: sv_counts
    tuple val(batch_name), path("filter_batch_sites_results/exposed/sv_count_plots/*"), emit: sv_count_plots
    tuple val(batch_name), path("filter_batch_sites_results/exposed/sites_filtered_depth.vcf*"), emit: sites_filtered_depth_vcf
    tuple val(batch_name), path("filter_batch_sites_results/exposed/sites_filtered_manta.vcf*"), emit: sites_filtered_manta_vcf
    tuple val(batch_name), path("filter_batch_sites_results/exposed/sites_filtered_wham.vcf*"), emit: sites_filtered_wham_vcf
    tuple val(batch_name), path("filter_batch_sites_results/exposed/sites_filtered_scramble.vcf*"), emit: sites_filtered_scramble_vcf
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

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.filterbatchsites_wdl} \\
        -i filter_batch_sites_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p filter_batch_sites_results
    cp filter_batch_sites_inputs.json filter_batch_sites_results/
    find cromwell-executions/FilterBatchSites/ -name "call-*" -type d -exec cp -r {} filter_batch_sites_results/ \\;

    mkdir -p filter_batch_sites_results/exposed/sv_counts
    mkdir -p filter_batch_sites_results/exposed/sv_count_plots

    cutoff_file=\$(find filter_batch_sites_results -type f -name "*.tsv" | grep -E -i "cutoff" | head -n 1 || true)
    if [[ -z "\${cutoff_file}" ]]; then
        echo "ERROR: FilterBatchSites could not resolve cutoffs output" >&2
        exit 1
    fi
    cp "\${cutoff_file}" filter_batch_sites_results/exposed/cutoffs.tsv

    count_files=\$(find filter_batch_sites_results -type f -name "*.tsv" | grep -E -i "count" || true)
    if [[ -z "\${count_files}" ]]; then
        echo "ERROR: FilterBatchSites could not resolve sv_counts outputs" >&2
        exit 1
    fi
    while IFS= read -r count_file; do
        [[ -n "\${count_file}" ]] && cp "\${count_file}" "filter_batch_sites_results/exposed/sv_counts/\$(basename "\${count_file}")"
    done <<< "\${count_files}"

    plot_files=\$(find filter_batch_sites_results -type f \\( -name "*.png" -o -name "*.pdf" \\) | grep -E -i "count|sv" || true)
    if [[ -z "\${plot_files}" ]]; then
        echo "ERROR: FilterBatchSites could not resolve sv_count_plots outputs" >&2
        exit 1
    fi
    while IFS= read -r plot_file; do
        [[ -n "\${plot_file}" ]] && cp "\${plot_file}" "filter_batch_sites_results/exposed/sv_count_plots/\$(basename "\${plot_file}")"
    done <<< "\${plot_files}"

    copy_sites_filtered_vcf() {
        local label="\$1"
        local candidate
        candidate=\$(find filter_batch_sites_results -type f \\( -name "*.vcf.gz" -o -name "*.vcf" \\) | grep -E -i "\${label}" | grep -E -i "sites[_-]?filtered|site[_-]?filtered|filtered" | head -n 1 || true)
        if [[ -z "\${candidate}" ]]; then
            echo "ERROR: FilterBatchSites could not resolve sites_filtered_\${label}.vcf output" >&2
            exit 1
        fi
        local suffix=".vcf"
        [[ "\${candidate}" == *.vcf.gz ]] && suffix=".vcf.gz"
        cp "\${candidate}" "filter_batch_sites_results/exposed/sites_filtered_\${label}\${suffix}"
    }

    copy_sites_filtered_vcf depth
    copy_sites_filtered_vcf manta
    copy_sites_filtered_vcf wham
    copy_sites_filtered_vcf scramble

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p filter_batch_sites_results/call-stub
    mkdir -p filter_batch_sites_results/exposed/sv_counts
    mkdir -p filter_batch_sites_results/exposed/sv_count_plots
    touch filter_batch_sites_results/call-stub/.stub
    touch filter_batch_sites_results/exposed/cutoffs.tsv
    touch filter_batch_sites_results/exposed/sv_counts/depth.counts.tsv
    touch filter_batch_sites_results/exposed/sv_count_plots/depth.counts.png
    touch filter_batch_sites_results/exposed/sites_filtered_depth.vcf.gz
    touch filter_batch_sites_results/exposed/sites_filtered_manta.vcf.gz
    touch filter_batch_sites_results/exposed/sites_filtered_wham.vcf.gz
    touch filter_batch_sites_results/exposed/sites_filtered_scramble.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
