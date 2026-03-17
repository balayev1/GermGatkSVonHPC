#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_MERGEBATCHSITES {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(ploidy_tables), path(depth_vcfs), path(pesr_vcfs)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.merge_batch_sites.vcf.gz"), emit: merge_batch_sites_vcf
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.merge_batch_sites.vcf.gz.tbi"), emit: merge_batch_sites_vcf_index
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

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/MergeBatchSites -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected MergeBatchSites output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.merge_batch_sites.vcf.gz" 1
    copy_outputs "${cohort_id}.merge_batch_sites.vcf.gz.tbi" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.merge_batch_sites.vcf.gz
    touch ${cohort_name}/${cohort_name}.merge_batch_sites.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
