#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_CLEANVCF {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(complex_genotype_vcfs), path(complex_resolve_bothside_pass_list), path(complex_resolve_background_fail_list), path(ped_file)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.cleaned.vcf.gz"), emit: cleaned_vcf
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.cleaned.vcf.gz.tbi"), emit: cleaned_vcf_index
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.metrics.tsv"), emit: metrics_file_makecohortvcf, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.cleanvcf_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.clean_vcf ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("CleanVcf cohort is required")
    }

    def asList = { value ->
        if (value == null) {
            return []
        }
        value instanceof List ? value : [value]
    }
    def toRealPathList = { value -> asList(value).collect { it.toRealPath().toString() } }

    def complexGenotypeVcfs = toRealPathList(complex_genotype_vcfs)
    if (!complexGenotypeVcfs) {
        throw new IllegalArgumentException("CleanVcf requires non-empty complex genotype VCF inputs")
    }

    def dynamic = [
        "CleanVcf.cohort_name"                         : cohort_id,
        "CleanVcf.ped_file"                            : ped_file.toRealPath().toString(),
        "CleanVcf.complex_genotype_vcfs"               : complexGenotypeVcfs,
        "CleanVcf.complex_resolve_bothside_pass_list"  : complex_resolve_bothside_pass_list.toRealPath().toString(),
        "CleanVcf.complex_resolve_background_fail_list": complex_resolve_background_fail_list.toRealPath().toString()
    ]
    file("clean_vcf_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV CLEANVCF] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \
        --template ${template_path} \
        --out clean_vcf_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file clean_vcf_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.cleanvcf_wdl} \
        -i clean_vcf_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/CleanVcf -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected CleanVcf output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.cleaned.vcf.gz" 1
    copy_outputs "${cohort_id}.cleaned.vcf.gz.tbi" 1
    copy_outputs "${cohort_id}.metrics.tsv" 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.cleaned.vcf.gz
    touch ${cohort_name}/${cohort_name}.cleaned.vcf.gz.tbi
    touch ${cohort_name}/${cohort_name}.metrics.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
