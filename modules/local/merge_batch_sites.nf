#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_MERGEBATCHSITES {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(ploidy_tables), path(depth_vcfs), path(pesr_vcfs)

    output:
    tuple val(cohort_name), path("**/*.merge_batch_sites.vcf.gz"), emit: merge_batch_sites_vcf
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.mergebatchsites_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.merge_batch_sites ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("MergeBatchSites cohort is required")
    }

    def asList = { value ->
        if (value == null) {
            return []
        }
        value instanceof List ? value : [value]
    }
    def toRealPathList = { value -> asList(value).collect { it.toRealPath().toString() } }

    def dynamic = [
        "MergeBatchSites.cohort"       : cohort_id,
        "MergeBatchSites.ploidy_tables": toRealPathList(ploidy_tables),
        "MergeBatchSites.depth_vcfs"   : toRealPathList(depth_vcfs),
        "MergeBatchSites.pesr_vcfs"    : toRealPathList(pesr_vcfs)
    ]
    file("merge_batch_sites_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV MERGEBATCHSITES] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out merge_batch_sites_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file merge_batch_sites_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.mergebatchsites_wdl} \\
        -i merge_batch_sites_inputs.json \\
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
    touch ${cohort_name}.merge_batch_sites.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
