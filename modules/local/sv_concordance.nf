#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_SVCONCORDANCE {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(eval_vcf), path(truth_vcf)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.concordance.vcf.gz"), emit: concordance_vcf
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.concordance.vcf.gz.tbi"), emit: concordance_vcf_index
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.svconcordance_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.sv_concordance ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("SVConcordance cohort name is required")
    }

    def dynamic = [
        "SVConcordance.eval_vcf"      : eval_vcf.toRealPath().toString(),
        "SVConcordance.truth_vcf"     : truth_vcf.toRealPath().toString(),
        "SVConcordance.output_prefix" : cohort_id
    ]
    file("sv_concordance_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV SVCONCORDANCE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \
        --template ${template_path} \
        --out sv_concordance_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file sv_concordance_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.svconcordance_wdl} \
        -i sv_concordance_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/SVConcordance -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected SVConcordance output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.concordance.vcf.gz" 1
    copy_outputs "${cohort_id}.concordance.vcf.gz.tbi" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.concordance.vcf.gz
    touch ${cohort_name}/${cohort_name}.concordance.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
