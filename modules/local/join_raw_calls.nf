#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_JOINRAWCALLS {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(ped_file), val(clustered_depth_vcfs), val(clustered_depth_vcf_indexes), val(clustered_dragen_vcfs), val(clustered_dragen_vcf_indexes), val(clustered_manta_vcfs), val(clustered_manta_vcf_indexes), val(clustered_melt_vcfs), val(clustered_melt_vcf_indexes), val(clustered_scramble_vcfs), val(clustered_scramble_vcf_indexes), val(clustered_wham_vcfs), val(clustered_wham_vcf_indexes)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.join_raw_calls.vcf.gz"), emit: joined_raw_calls_vcf
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.join_raw_calls.vcf.gz.tbi"), emit: joined_raw_calls_vcf_index
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.ploidy.tsv"), emit: ploidy_table
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.joinrawcalls_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.join_raw_calls ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("JoinRawCalls cohort name is required")
    }

    def asList = { value ->
        if (value == null) {
            return []
        }
        value instanceof List ? value : [value]
    }
    def toOptionalRealPathList = { value ->
        def items = asList(value).findAll { it != null && it.toString().trim() }
        if (!items) {
            return null
        }
        items.collect { file(it.toString()).toRealPath().toString() }
    }

    def dynamic = [
        "JoinRawCalls.prefix"                        : cohort_id,
        "JoinRawCalls.ped_file"                      : ped_file.toRealPath().toString(),
        "JoinRawCalls.clustered_depth_vcfs"          : toOptionalRealPathList(clustered_depth_vcfs),
        "JoinRawCalls.clustered_depth_vcf_indexes"   : toOptionalRealPathList(clustered_depth_vcf_indexes),
        "JoinRawCalls.clustered_dragen_vcfs"         : toOptionalRealPathList(clustered_dragen_vcfs),
        "JoinRawCalls.clustered_dragen_vcf_indexes"  : toOptionalRealPathList(clustered_dragen_vcf_indexes),
        "JoinRawCalls.clustered_manta_vcfs"          : toOptionalRealPathList(clustered_manta_vcfs),
        "JoinRawCalls.clustered_manta_vcf_indexes"   : toOptionalRealPathList(clustered_manta_vcf_indexes),
        "JoinRawCalls.clustered_melt_vcfs"           : toOptionalRealPathList(clustered_melt_vcfs),
        "JoinRawCalls.clustered_melt_vcf_indexes"    : toOptionalRealPathList(clustered_melt_vcf_indexes),
        "JoinRawCalls.clustered_scramble_vcfs"       : toOptionalRealPathList(clustered_scramble_vcfs),
        "JoinRawCalls.clustered_scramble_vcf_indexes": toOptionalRealPathList(clustered_scramble_vcf_indexes),
        "JoinRawCalls.clustered_wham_vcfs"           : toOptionalRealPathList(clustered_wham_vcfs),
        "JoinRawCalls.clustered_wham_vcf_indexes"    : toOptionalRealPathList(clustered_wham_vcf_indexes)
    ]
    file("join_raw_calls_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV JOINRAWCALLS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out join_raw_calls_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file join_raw_calls_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.joinrawcalls_wdl} \\
        -i join_raw_calls_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/JoinRawCalls -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected JoinRawCalls output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.join_raw_calls.vcf.gz" 1
    copy_outputs "${cohort_id}.join_raw_calls.vcf.gz.tbi" 1
    copy_outputs "${cohort_id}.ploidy.tsv" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.join_raw_calls.vcf.gz
    touch ${cohort_name}/${cohort_name}.join_raw_calls.vcf.gz.tbi
    touch ${cohort_name}/${cohort_name}.ploidy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
