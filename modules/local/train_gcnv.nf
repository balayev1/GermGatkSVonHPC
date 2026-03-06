#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATK_TRAINGCNV {
    tag "${cohort}"
    label 'process_medium'

    input:
    tuple val(cohort), val(sample_ids), path(count_files)  // (Mandatory!) channel: [cohort, sample_ids, count_files]

    output:
    tuple val(cohort), path("**/${cohort}-contig-ploidy-model.tar.gz"), emit: cohort_contig_ploidy_model_tar
    tuple val(cohort), path("**/${cohort}-gcnv-model-shard-*.tar.gz"), emit: cohort_gcnv_model_tars
    tuple val(cohort), path("**/${cohort}-contig-ploidy-calls.tar.gz"), emit: cohort_contig_ploidy_calls_tar
    tuple val(cohort), path("**/${cohort}-gcnv-calls-shard-*.tar.gz"), emit: cohort_gcnv_calls_tars
    tuple val(cohort), path("**/genotyped-segments-*.vcf.gz"), emit: cohort_genotyped_segments_vcfs
    tuple val(cohort), path("**/${cohort}-gcnv-tracking-shard*.tar.gz"), emit: cohort_gcnv_tracking_tars
    tuple val(cohort), path("**/genotyped-intervals-*.vcf.gz"), emit: cohort_genotyped_intervals_vcfs
    tuple val(cohort), path("**/denoised_copy_ratios-*.tsv"), emit: cohort_denoised_copy_ratios
    tuple val(cohort), path("**/*.annotated.tsv"), emit: annotated_intervals, optional: true
    tuple val(cohort), path("**/*.filtered.interval_list"), emit: filtered_intervals_cnv, optional: true
    tuple val(cohort), path("**/*.filtered.interval_list"), emit: filtered_intervals_ploidy, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.traingcnv_template).toAbsolutePath()
    def static_map = params.tool_inputs?.traingcnv ?: [:]
    def static_json = JsonOutput.toJson(static_map)
    def cohort_name = cohort?.toString()?.trim()
    if (!cohort_name) {
        throw new IllegalArgumentException("TrainGCNV cohort is required")
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
        "TrainGCNV.cohort"     : cohort_name
    ]
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

    mkdir -p train_gcnv_results
    cp train_gcnv_inputs.json train_gcnv_results/
    find cromwell-executions/TrainGCNV/ -name "call-*" -type d -exec cp -r {} train_gcnv_results/ \\;
    touch train_gcnv_results/placeholder.gcnv_model_tars.tar
    touch train_gcnv_results/placeholder.contig_ploidy_model_tar.tar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p train_gcnv_results/call-stub
    touch train_gcnv_results/call-stub/.stub
    touch train_gcnv_results/placeholder.gcnv_model_tars.tar
    touch train_gcnv_results/placeholder.contig_ploidy_model_tar.tar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
