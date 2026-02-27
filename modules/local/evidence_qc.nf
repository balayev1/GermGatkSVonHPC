#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_EVIDENCEQC {
    tag "${cohort}"
    label 'process_medium'
    
    input:
    tuple val(cohort), val(sample_ids), path(count_files), path(manta_vcf_files), path(wham_vcf_files), path(scramble_vcf_files)

    output:
    tuple val(cohort), val(sample_ids), path("evidence_qc_results"), emit: evidence_qc_results
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.evidqc_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.evidence_qc ?: [:])
    def cohort_name = cohort?.toString()?.trim()
    if (!cohort_name) {
        throw new IllegalArgumentException("EvidenceQC cohort is required and must come from samplesheet metadata (meta.cohort)")
    }
    def asList = { value -> value instanceof List ? value : [value] }
    def samples = asList(sample_ids).collect { it.toString() }
    def counts = asList(count_files).collect { it.toRealPath().toString() }
    def manta = asList(manta_vcf_files).collect { it.toRealPath().toString() }
    def wham = asList(wham_vcf_files).collect { it.toRealPath().toString() }
    def scramble = asList(scramble_vcf_files).collect { it.toRealPath().toString() }

    if (!samples) {
        throw new IllegalArgumentException("EvidenceQC received an empty sample list")
    }
    def n = samples.size()
    if (counts.size() != n || manta.size() != n || wham.size() != n || scramble.size() != n) {
        throw new IllegalArgumentException("EvidenceQC input size mismatch: samples=${n}, counts=${counts.size()}, manta=${manta.size()}, wham=${wham.size()}, scramble=${scramble.size()}")
    }

    def dynamic = [
        "EvidenceQC.batch"        : cohort_name,
        "EvidenceQC.samples"      : samples,
        "EvidenceQC.counts"       : counts,
        "EvidenceQC.manta_vcfs"   : manta,
        "EvidenceQC.wham_vcfs"    : wham,
        "EvidenceQC.scramble_vcfs": scramble
    ]
    file("evidence_qc_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV EVIDENCEQC] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    
    """
    render_json_template.py \\
        --template ${template_path} \\
        --out evidence_qc_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file evidence_qc_dynamic.json

    # Execute Cromwell
    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.evidqc_wdl} \\
        -i evidence_qc_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p evidence_qc_results
    cp evidence_qc_inputs.json evidence_qc_results/
    find cromwell-executions/EvidenceQC/ -name "call-*" -type d -exec cp -r {} evidence_qc_results/ \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p evidence_qc_results/call-stub
    touch evidence_qc_results/call-stub/.stub
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
