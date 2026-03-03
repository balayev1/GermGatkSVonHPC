#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATK_TRAINGCNV {
    tag "${cohort}"
    label 'process_medium'

    input:
    tuple val(cohort), val(sample_ids), path(count_files)  // (Mandatory!) channel: [cohort, sample_ids, count_files]


    output:
    tuple val(cohort), path("train_gcnv_results"), emit: train_gcnv_results
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.traingcnv_template).toAbsolutePath()
    def static_map = params.tool_inputs?.traingcnv ?: [:]
    def static_json = JsonOutput.toJson(static_map)
    def cohort_name = cohort?.toString()?.trim()
    if (!cohort_name) {
        throw new IllegalArgumentException("TrainGCNV cohort is required")
    }
    def required_static_keys = [
        "TrainGCNV.reference_fasta",
        "TrainGCNV.reference_index",
        "TrainGCNV.reference_dict",
        "TrainGCNV.ref_copy_number_autosomal_contigs",
        "TrainGCNV.allosomal_contigs",
        "TrainGCNV.contig_ploidy_priors",
        "TrainGCNV.exclude_intervals_for_filter_intervals_ploidy",
        "TrainGCNV.exclude_intervals_for_filter_intervals_cnv",
        "TrainGCNV.num_intervals_per_scatter",
        "TrainGCNV.sv_base_mini_docker",
        "TrainGCNV.linux_docker",
        "TrainGCNV.gatk_docker"
    ]
    required_static_keys.each { key ->
        if (!static_map.containsKey(key) || static_map[key] == null || static_map[key].toString().trim() == "") {
            throw new IllegalArgumentException("TrainGCNV required config missing: ${key}")
        }
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

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.traingcnv_wdl} \\
        -i train_gcnv_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p train_gcnv_results
    cp train_gcnv_inputs.json train_gcnv_results/
    find cromwell-executions/TrainGCNV/ -name "call-*" -type d -exec cp -r {} train_gcnv_results/ \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p train_gcnv_results/call-stub
    touch train_gcnv_results/call-stub/.stub

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
