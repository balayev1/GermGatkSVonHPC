#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_CLUSTERBATCH {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(ped_file), path(del_bed), path(dup_bed), path(wham_vcf_tar), path(manta_vcf_tar), path(scramble_vcf_tar), val(n_iqr_cutoff_plotting)

    output:
    tuple val(batch_name), path("${batch_name}/${batch_name}.cluster_batch.depth.vcf.gz"), emit: clustered_depth_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.cluster_batch.depth.vcf.gz.tbi"), emit: clustered_depth_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.cluster_batch.manta.vcf.gz"), emit: clustered_manta_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.cluster_batch.manta.vcf.gz.tbi"), emit: clustered_manta_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.cluster_batch.wham.vcf.gz"), emit: clustered_wham_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.cluster_batch.wham.vcf.gz.tbi"), emit: clustered_wham_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.cluster_batch.scramble.vcf.gz"), emit: clustered_scramble_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.cluster_batch.scramble.vcf.gz.tbi"), emit: clustered_scramble_vcf_index
    tuple val(batch_name), path("${batch_name}/*.svcounts.txt"), emit: clustered_sv_counts, optional: true
    tuple val(batch_name), path("${batch_name}/*.all_SVTYPEs.counts_per_sample.png"), emit: clustered_sv_count_plots, optional: true
    tuple val(batch_name), path("${batch_name}/${batch_name}.outliers_preview.samples.txt"), emit: clustered_outlier_samples_preview, optional: true
    tuple val(batch_name), path("${batch_name}/${batch_name}.outliers_preview_with_reason.samples.tsv"), emit: clustered_outlier_samples_with_reason, optional: true
    tuple val(batch_name), path("${batch_name}/num_outliers.txt"), emit: num_outlier_samples, optional: true
    tuple val(batch_name), path("${batch_name}/ClusterBatch.${batch_name}.metrics.tsv"), emit: metrics_file_clusterbatch, optional: true   
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

    mkdir -p "${batch_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${batch_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/ClusterBatch -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected ClusterBatch output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${batch_id}.cluster_batch.depth.vcf.gz" 1
    copy_outputs "${batch_id}.cluster_batch.depth.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.cluster_batch.manta.vcf.gz" 1
    copy_outputs "${batch_id}.cluster_batch.manta.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.cluster_batch.wham.vcf.gz" 1
    copy_outputs "${batch_id}.cluster_batch.wham.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.cluster_batch.scramble.vcf.gz" 1
    copy_outputs "${batch_id}.cluster_batch.scramble.vcf.gz.tbi" 1
    copy_outputs "*.svcounts.txt" 0
    copy_outputs "*.all_SVTYPEs.counts_per_sample.png" 0
    copy_outputs "${batch_id}.outliers_preview.samples.txt" 0
    copy_outputs "${batch_id}.outliers_preview_with_reason.samples.tsv" 0
    copy_outputs "num_outliers.txt" 0
    copy_outputs "ClusterBatch.${batch_id}.metrics.tsv" 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${batch_name}
    touch ${batch_name}/${batch_name}.cluster_batch.depth.vcf.gz
    touch ${batch_name}/${batch_name}.cluster_batch.depth.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.cluster_batch.manta.vcf.gz
    touch ${batch_name}/${batch_name}.cluster_batch.manta.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.cluster_batch.wham.vcf.gz
    touch ${batch_name}/${batch_name}.cluster_batch.wham.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.cluster_batch.scramble.vcf.gz
    touch ${batch_name}/${batch_name}.cluster_batch.scramble.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.svcounts.txt
    touch ${batch_name}/${batch_name}.all_SVTYPEs.counts_per_sample.png
    touch ${batch_name}/${batch_name}.outliers_preview.samples.txt
    touch ${batch_name}/${batch_name}.outliers_preview_with_reason.samples.tsv
    touch ${batch_name}/num_outliers.txt
    touch ${batch_name}/ClusterBatch.${batch_name}.metrics.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
