#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_FILTERGENOTYPES {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(vcf), path(ploidy_table), path(ped_file), val(optimized_sl_cutoff_table)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.filter_genotypes.sanitized.vcf.gz"), emit: filtered_vcf
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.filter_genotypes.sanitized.vcf.gz.tbi"), emit: filtered_vcf_index
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.filter_genotypes_SV_VCF_QC_output.tar.gz"), emit: main_vcf_qc_tarball, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.filtergenotypes_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.filter_genotypes ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("FilterGenotypes cohort name is required")
    }

    def optimizedSlCutoffTable = optimized_sl_cutoff_table?.toString()?.trim()
    def dynamic = [
        "FilterGenotypes.vcf"          : vcf.toRealPath().toString(),
        "FilterGenotypes.output_prefix": cohort_id,
        "FilterGenotypes.ploidy_table" : ploidy_table.toRealPath().toString(),
        "FilterGenotypes.ped_file"     : ped_file.toRealPath().toString()
    ]
    if (optimizedSlCutoffTable) {
        dynamic["FilterGenotypes.optimized_sl_cutoff_table"] = file(optimizedSlCutoffTable).toRealPath().toString()
    }
    file("filter_genotypes_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV FILTERGENOTYPES] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \
        --template ${template_path} \
        --out filter_genotypes_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file filter_genotypes_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.filtergenotypes_wdl} \
        -i filter_genotypes_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/FilterGenotypes -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected FilterGenotypes output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.filter_genotypes.sanitized.vcf.gz" 1
    copy_outputs "${cohort_id}.filter_genotypes.sanitized.vcf.gz.tbi" 1
    copy_outputs "${cohort_id}.filter_genotypes_SV_VCF_QC_output.tar.gz" 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.filter_genotypes.sanitized.vcf.gz
    touch ${cohort_name}/${cohort_name}.filter_genotypes.sanitized.vcf.gz.tbi
    touch ${cohort_name}/${cohort_name}.filter_genotypes_SV_VCF_QC_output.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
