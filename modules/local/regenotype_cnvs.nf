#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_REGENOTYPECNVS {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(depth_vcfs), path(merge_batch_sites_vcf), path(batch_depth_vcfs), path(coveragefiles), path(coveragefile_idxs), path(medianfiles), path(genotyping_rd_table), path(ploidy_tables), val(batches), path(regeno_coverage_medians)

    output:
    tuple val(cohort_name), path("${cohort_name}/*.depth.regeno_final.vcf.gz"), emit: regenotyped_depth_vcfs
    tuple val(cohort_name), path("${cohort_name}/*.depth.regeno_final.vcf.gz.tbi"), emit: regenotyped_depth_vcf_indexes
    tuple val(cohort_name), path("${cohort_name}/regeno_num_lines.txt"), emit: number_regenotyped_file
    tuple val(cohort_name), path("${cohort_name}/regeno_filtered_num_lines.txt"), emit: number_regenotyped_filtered_file
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
    render_json_template.py \
        --template ${template_path} \
        --out regenotype_cnvs_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file regenotype_cnvs_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.regenotypecnvs_wdl} \
        -i regenotype_cnvs_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/RegenotypeCNVs -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected RegenotypeCNVs output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "*.depth.regeno_final.vcf.gz" 1
    copy_outputs "*.depth.regeno_final.vcf.gz.tbi" 1
    copy_outputs "regeno_num_lines.txt" 1
    copy_outputs "regeno_filtered_num_lines.txt" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.depth.regeno_final.vcf.gz
    touch ${cohort_name}/${cohort_name}.depth.regeno_final.vcf.gz.tbi
    touch ${cohort_name}/regeno_num_lines.txt
    touch ${cohort_name}/regeno_filtered_num_lines.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
