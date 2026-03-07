#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_CLUSTERBATCH {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(ped_file), path(del_bed), path(dup_bed), path(wham_vcf_tar), path(manta_vcf_tar), path(scramble_vcf_tar), val(n_iqr_cutoff_plotting)

    output:
    tuple val(batch_name), path("**/${batch_name}.cluster_batch.depth.vcf.gz"), emit: clustered_depth_vcf
    tuple val(batch_name), path("**/${batch_name}.cluster_batch.manta.vcf.gz"), emit: clustered_manta_vcf
    tuple val(batch_name), path("**/${batch_name}.cluster_batch.wham.vcf.gz"), emit: clustered_wham_vcf
    tuple val(batch_name), path("**/${batch_name}.cluster_batch.scramble.vcf.gz"), emit: clustered_scramble_vcf
    tuple val(batch_name), path("**/*.svcounts.txt"), emit: clustered_sv_counts, optional: true
    tuple val(batch_name), path("**/*.all_SVTYPEs.counts_per_sample.png"), emit: clustered_sv_count_plots, optional: true
    tuple val(batch_name), path("**/*.outliers_preview.samples.txt"), emit: clustered_outlier_samples_preview, optional: true
    tuple val(batch_name), path("**/*.outliers_preview_with_reason.samples.tsv"), emit: clustered_outlier_samples_with_reason, optional: true
    tuple val(batch_name), path("num_outliers.txt"), emit: num_outlier_samples, optional: true
   
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
    if (n_iqr_cutoff_plotting != null && n_iqr_cutoff_plotting.toString().trim()) {
        dynamic["ClusterBatch.N_IQR_cutoff_plotting"] = n_iqr_cutoff_plotting.toString().toInteger()
    }
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

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.clusterbatch_wdl} \\
        -i cluster_batch_inputs.json \\
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
    touch ${batch_name}.cluster_batch.depth.vcf.gz
    touch ${batch_name}.cluster_batch.manta.vcf.gz
    touch ${batch_name}.cluster_batch.wham.vcf.gz
    touch ${batch_name}.cluster_batch.scramble.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
