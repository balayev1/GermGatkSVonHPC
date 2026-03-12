#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATK_TRAINGCNV {
    tag "${batch_key}"
    label 'process_medium'

    input:
    tuple val(batch_key), val(sample_ids), path(count_files), val(outlier_sample_ids)  // channel: [batch_key, sample_ids, count_files, outlier_sample_ids?]

    output:
    tuple val(batch_key), path("${batch_key}/${batch_key}-contig-ploidy-model.tar.gz"), emit: cohort_contig_ploidy_model_tar
    tuple val(batch_key), path("${batch_key}/${batch_key}-gcnv-model-shard-*.tar.gz"), emit: cohort_gcnv_model_tars
    tuple val(batch_key), path("${batch_key}/${batch_key}-contig-ploidy-calls.tar.gz"), emit: cohort_contig_ploidy_calls_tar
    tuple val(batch_key), path("${batch_key}/${batch_key}-gcnv-calls-shard-*.tar.gz"), emit: cohort_gcnv_calls_tars
    tuple val(batch_key), path("${batch_key}/genotyped-segments-*.vcf.gz"), emit: cohort_genotyped_segments_vcfs
    tuple val(batch_key), path("${batch_key}/${batch_key}-gcnv-tracking-shard-*.tar.gz"), emit: cohort_gcnv_tracking_tars
    tuple val(batch_key), path("${batch_key}/genotyped-intervals-*.vcf.gz"), emit: cohort_genotyped_intervals_vcfs
    tuple val(batch_key), path("${batch_key}/denoised_copy_ratios-*.tsv"), emit: cohort_denoised_copy_ratios
    tuple val(batch_key), path("${batch_key}/*.annotated.tsv"), emit: annotated_intervals, optional: true
    tuple val(batch_key), path("${batch_key}/*.filtered.interval_list"), emit: filtered_intervals_cnv, optional: true
    tuple val(batch_key), path("${batch_key}/*.filtered.interval_list"), emit: filtered_intervals_ploidy, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.traingcnv_template).toAbsolutePath()
    def static_map = params.tool_inputs?.traingcnv ?: [:]
    def static_json = JsonOutput.toJson(static_map)
    def train_batch_key = batch_key?.toString()?.trim()
    if (!train_batch_key) {
        throw new IllegalArgumentException("TrainGCNV batch key is required")
    }

    def asList = { value -> value instanceof List ? value : [value] }
    def samples = asList(sample_ids).collect { it.toString() }
    def counts = asList(count_files).collect { it.toRealPath().toString() }
    if (!samples) {
        throw new IllegalArgumentException("TrainGCNV received an empty sample list")
    }
    if (samples.size() != counts.size()) {
        throw new IllegalArgumentException("TrainGCNV input size mismatch: samples=${samples.size()}, counts=${counts.size()}")
    }

    def dynamic = [
        "TrainGCNV.samples"    : samples,
        "TrainGCNV.count_files": counts,
        "TrainGCNV.cohort"     : train_batch_key
    ]
    if (outlier_sample_ids != null && outlier_sample_ids.toString().trim()) {
        dynamic["TrainGCNV.outlier_sample_ids"] = file(outlier_sample_ids.toString()).toRealPath().toString()
    }
    file("train_gcnv_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK TRAINGCNV] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out train_gcnv_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file train_gcnv_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.traingcnv_wdl} \\
        -i train_gcnv_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p "${train_batch_key}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${train_batch_key}/\$(basename "\$source")"
        done < <(find cromwell-executions/TrainGCNV -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected TrainGCNV output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${train_batch_key}-contig-ploidy-model.tar.gz" 1
    copy_outputs "${train_batch_key}-gcnv-model-shard-*.tar.gz" 1
    copy_outputs "${train_batch_key}-contig-ploidy-calls.tar.gz" 1
    copy_outputs "${train_batch_key}-gcnv-calls-shard-*.tar.gz" 1
    copy_outputs "genotyped-segments-*.vcf.gz" 1
    copy_outputs "${train_batch_key}-gcnv-tracking-shard-*.tar.gz" 1
    copy_outputs "genotyped-intervals-*.vcf.gz" 1
    copy_outputs "denoised_copy_ratios-*.tsv" 1
    copy_outputs "*.annotated.tsv" 0
    copy_outputs "*.filtered.interval_list" 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${batch_key}
    touch ${batch_key}/${batch_key}-contig-ploidy-model.tar.gz
    touch ${batch_key}/${batch_key}-gcnv-model-shard-1.tar.gz
    touch ${batch_key}/${batch_key}-contig-ploidy-calls.tar.gz
    touch ${batch_key}/${batch_key}-gcnv-calls-shard-1.tar.gz
    touch ${batch_key}/genotyped-segments-1.vcf.gz
    touch ${batch_key}/${batch_key}-gcnv-tracking-shard-1.tar.gz
    touch ${batch_key}/genotyped-intervals-1.vcf.gz
    touch ${batch_key}/denoised_copy_ratios-1.tsv
    touch ${batch_key}/${batch_key}.annotated.tsv
    touch ${batch_key}/${batch_key}.filtered.interval_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
