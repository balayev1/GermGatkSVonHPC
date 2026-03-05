#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_CLUSTERBATCH {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(ped_file), path(del_bed), path(dup_bed), path(wham_vcf_tar), path(manta_vcf_tar), path(scramble_vcf_tar)

    output:
    tuple val(batch_name), path("cluster_batch_results"), emit: cluster_batch_results
    tuple val(batch_name), path("cluster_batch_results/exposed_vcfs/clustered_*.vcf*"), emit: clustered_vcfs, optional: true
    tuple val(batch_name), path("cluster_batch_results/exposed_vcfs/clustered_depth.vcf*"), emit: clustered_depth_vcf
    tuple val(batch_name), path("cluster_batch_results/exposed_vcfs/clustered_manta.vcf*"), emit: clustered_manta_vcf
    tuple val(batch_name), path("cluster_batch_results/exposed_vcfs/clustered_wham.vcf*"), emit: clustered_wham_vcf
    tuple val(batch_name), path("cluster_batch_results/exposed_vcfs/clustered_scramble.vcf*"), emit: clustered_scramble_vcf
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.clusterbatch_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.cluster_batch ?: [:])

    def batch_id = batch_name?.toString()?.trim()
    if (!batch_id) {
        throw new IllegalArgumentException("ClusterBatch batch name is required")
    }
    def cohort_id = cohort?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("ClusterBatch cohort is required")
    }

    def dynamic = [
        "ClusterBatch.batch"            : batch_id,
        "ClusterBatch.del_bed"          : del_bed.toRealPath().toString(),
        "ClusterBatch.dup_bed"          : dup_bed.toRealPath().toString(),
        "ClusterBatch.wham_vcf_tar"     : wham_vcf_tar.toRealPath().toString(),
        "ClusterBatch.manta_vcf_tar"    : manta_vcf_tar.toRealPath().toString(),
        "ClusterBatch.scramble_vcf_tar" : scramble_vcf_tar.toRealPath().toString(),
        "ClusterBatch.ped_file"         : ped_file.toRealPath().toString()
    ]
    file("cluster_batch_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV CLUSTERBATCH] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out cluster_batch_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file cluster_batch_dynamic.json

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.clusterbatch_wdl} \\
        -i cluster_batch_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p cluster_batch_results
    cp cluster_batch_inputs.json cluster_batch_results/
    find cromwell-executions/ClusterBatch/ -name "call-*" -type d -exec cp -r {} cluster_batch_results/ \\;
    mkdir -p cluster_batch_results/exposed_vcfs

    choose_and_copy_vcf() {
        local label="\$1"
        local regex="\$2"
        local chosen
        chosen=\$(find cluster_batch_results -type f \\( -name "*.vcf.gz" -o -name "*.vcf" \\) | grep -E -i "\${regex}" | head -n 1 || true)
        if [[ -n "\${chosen}" ]]; then
            local suffix=".vcf"
            if [[ "\${chosen}" == *.vcf.gz ]]; then
                suffix=".vcf.gz"
            fi
            cp "\${chosen}" "cluster_batch_results/exposed_vcfs/clustered_\${label}\${suffix}"
        fi
    }

    choose_and_copy_vcf depth "depth|cnmops|gcnv"
    choose_and_copy_vcf manta "manta"
    choose_and_copy_vcf wham "wham"
    choose_and_copy_vcf scramble "scramble"

    for required_caller in depth manta wham scramble; do
        if ! ls cluster_batch_results/exposed_vcfs/clustered_\${required_caller}.vcf* >/dev/null 2>&1; then
            echo "ERROR: ClusterBatch could not resolve clustered_\${required_caller}.vcf output" >&2
            exit 1
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cluster_batch_results/call-stub
    mkdir -p cluster_batch_results/exposed_vcfs
    touch cluster_batch_results/call-stub/.stub
    touch cluster_batch_results/exposed_vcfs/clustered_depth.vcf.gz
    touch cluster_batch_results/exposed_vcfs/clustered_manta.vcf.gz
    touch cluster_batch_results/exposed_vcfs/clustered_wham.vcf.gz
    touch cluster_batch_results/exposed_vcfs/clustered_scramble.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
