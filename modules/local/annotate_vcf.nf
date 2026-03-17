#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_ANNOTATEVCF {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(vcf), path(ped_file)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.annotated.vcf.gz"), emit: annotated_vcf
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.annotated.vcf.gz.tbi"), emit: annotated_vcf_index
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.annotatevcf_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.annotate_vcf ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("AnnotateVcf cohort name is required")
    }

    def dynamic = [
        "AnnotateVcf.vcf"       : vcf.toRealPath().toString(),
        "AnnotateVcf.prefix"    : cohort_id,
        "AnnotateVcf.ped_file"  : ped_file.toRealPath().toString()
    ]
    file("annotate_vcf_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV ANNOTATEVCF] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \
        --template ${template_path} \
        --out annotate_vcf_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file annotate_vcf_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.annotatevcf_wdl} \
        -i annotate_vcf_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/AnnotateVcf -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected AnnotateVcf output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.annotated.vcf.gz" 1
    copy_outputs "${cohort_id}.annotated.vcf.gz.tbi" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.annotated.vcf.gz
    touch ${cohort_name}/${cohort_name}.annotated.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
