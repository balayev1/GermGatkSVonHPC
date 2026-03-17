#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_SCOREGENOTYPES {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(vcf), val(truth_json)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.concordance.gq_recalibrated.vcf.gz"), emit: unfiltered_recalibrated_vcf
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.concordance.gq_recalibrated.vcf.gz.tbi"), emit: unfiltered_recalibrated_vcf_index
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.vcf_table.tsv.gz"), emit: vcf_optimization_table, optional: true
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.sl_optimization.sl_cutoff_table.tsv"), emit: sl_cutoff_table, optional: true
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.sl_optimization.tar.gz"), emit: sl_cutoff_qc_tarball, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.scoregenotypes_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.score_genotypes ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("ScoreGenotypes cohort name is required")
    }

    def truthJson = truth_json?.toString()?.trim()
    def dynamic = [
        "ScoreGenotypes.vcf"          : vcf.toRealPath().toString(),
        "ScoreGenotypes.output_prefix": cohort_id
    ]
    if (truthJson) {
        dynamic["ScoreGenotypes.truth_json"] = file(truthJson).toRealPath().toString()
    }
    file("score_genotypes_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def vcfName = vcf.getFileName().toString()
    def recalibratedVcf = vcfName.replaceFirst(/\.vcf\.gz$/, '.gq_recalibrated.vcf.gz')
    if (recalibratedVcf == vcfName) {
        throw new IllegalArgumentException("ScoreGenotypes expects a .vcf.gz input, got: ${vcfName}")
    }
    def recalibratedVcfIndex = "${recalibratedVcf}.tbi"

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV SCOREGENOTYPES] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \
        --template ${template_path} \
        --out score_genotypes_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file score_genotypes_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.scoregenotypes_wdl} \
        -i score_genotypes_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/ScoreGenotypes -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected ScoreGenotypes output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${recalibratedVcf}" 1
    copy_outputs "${recalibratedVcfIndex}" 1
    copy_outputs "${cohort_id}.vcf_table.tsv.gz" 0
    copy_outputs "${cohort_id}.sl_optimization.sl_cutoff_table.tsv" 0
    copy_outputs "${cohort_id}.sl_optimization.tar.gz" 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.concordance.gq_recalibrated.vcf.gz
    touch ${cohort_name}/${cohort_name}.concordance.gq_recalibrated.vcf.gz.tbi
    touch ${cohort_name}/${cohort_name}.vcf_table.tsv.gz
    touch ${cohort_name}/${cohort_name}.sl_optimization.sl_cutoff_table.tsv
    touch ${cohort_name}/${cohort_name}.sl_optimization.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
