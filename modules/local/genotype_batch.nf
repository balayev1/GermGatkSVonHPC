#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GENOTYPEBATCH {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(rf_cutoffs), path(pe_file), path(rd_file), path(sr_file), path(median_coverage), path(ploidy_table), path(vcf)

    output:
    tuple val(batch_name), path("genotype_batch_results"), emit: genotype_batch_results
    tuple val(batch_name), path("genotype_batch_results/exposed/genotyped.vcf*"), emit: genotyped_vcf, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.genotypebatch_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.genotype_batch ?: [:])

    def batch_id = batch_name?.toString()?.trim()
    if (!batch_id) {
        throw new IllegalArgumentException("GenotypeBatch batch is required")
    }
    def cohort_id = cohort?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("GenotypeBatch cohort is required")
    }

    def dynamic = [
        "GenotypeBatch.batch"           : batch_id,
        "GenotypeBatch.rf_cutoffs"      : rf_cutoffs.toRealPath().toString(),
        "GenotypeBatch.pe_file"         : pe_file.toRealPath().toString(),
        "GenotypeBatch.rd_file"         : rd_file.toRealPath().toString(),
        "GenotypeBatch.sr_file"         : sr_file.toRealPath().toString(),
        "GenotypeBatch.median_coverage" : median_coverage.toRealPath().toString(),
        "GenotypeBatch.vcf"             : vcf.toRealPath().toString(),
        "GenotypeBatch.ploidy_table"    : ploidy_table.toRealPath().toString()
    ]
    file("genotype_batch_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV GENOTYPEBATCH] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out genotype_batch_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file genotype_batch_dynamic.json

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.genotypebatch_wdl} \\
        -i genotype_batch_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p genotype_batch_results
    cp genotype_batch_inputs.json genotype_batch_results/
    find cromwell-executions/GenotypeBatch/ -name "call-*" -type d -exec cp -r {} genotype_batch_results/ \\;

    mkdir -p genotype_batch_results/exposed
    genotyped_vcf=\$(find genotype_batch_results -type f \\( -name "*.vcf.gz" -o -name "*.vcf" \\) | grep -E -i "genotyp|batch|cohort" | head -n 1 || true)
    if [[ -n "\${genotyped_vcf}" ]]; then
        suffix=".vcf"
        [[ "\${genotyped_vcf}" == *.vcf.gz ]] && suffix=".vcf.gz"
        cp "\${genotyped_vcf}" "genotype_batch_results/exposed/genotyped\${suffix}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p genotype_batch_results/call-stub
    mkdir -p genotype_batch_results/exposed
    touch genotype_batch_results/call-stub/.stub
    touch genotype_batch_results/exposed/genotyped.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
