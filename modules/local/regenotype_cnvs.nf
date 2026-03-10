#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_REGENOTYPECNVS {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(depth_vcfs), path(merge_batch_sites_vcf), path(batch_depth_vcfs), path(coveragefiles), path(coveragefile_idxs), path(medianfiles), path(genotyping_rd_table), path(ploidy_tables), val(batches), path(regeno_coverage_medians)

    output:
    tuple val(cohort_name), path("**/*.depth.regeno_final.vcf.gz"), emit: regenotyped_depth_vcfs
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.regenotypecnvs_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.regenotype_cnvs ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("RegenotypeCNVs cohort is required")
    }

    def asList = { value ->
        if (value == null) {
            return []
        }
        value instanceof List ? value : [value]
    }
    def toRealPathList = { value -> asList(value).collect { it.toRealPath().toString() } }

    def batches_list = asList(batches).collect { it.toString() }
    if (!batches_list) {
        throw new IllegalArgumentException("RegenotypeCNVs requires a non-empty batch list")
    }

    def dynamic = [
        "RegenotypeCNVs.cohort"                  : cohort_id,
        "RegenotypeCNVs.depth_vcfs"              : toRealPathList(depth_vcfs),
        "RegenotypeCNVs.merge_batch_sites_vcf"   : merge_batch_sites_vcf.toRealPath().toString(),
        "RegenotypeCNVs.batch_depth_vcfs"        : toRealPathList(batch_depth_vcfs),
        "RegenotypeCNVs.coveragefiles"           : toRealPathList(coveragefiles),
        "RegenotypeCNVs.coveragefile_idxs"       : toRealPathList(coveragefile_idxs),
        "RegenotypeCNVs.medianfiles"             : toRealPathList(medianfiles),
        "RegenotypeCNVs.genotyping_rd_table"     : toRealPathList(genotyping_rd_table),
        "RegenotypeCNVs.ploidy_tables"           : toRealPathList(ploidy_tables),
        "RegenotypeCNVs.batches"                 : batches_list,
        "RegenotypeCNVs.regeno_coverage_medians" : toRealPathList(regeno_coverage_medians)
    ]
    file("regenotype_cnvs_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV REGENOTYPECNVS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out regenotype_cnvs_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file regenotype_cnvs_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.regenotypecnvs_wdl} \\
        -i regenotype_cnvs_inputs.json \\
        -p ${params.deps_zip}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p regenotype_cnvs_results/
    touch regenotype_cnvs_results/${cohort_name}.depth.regeno_final.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
