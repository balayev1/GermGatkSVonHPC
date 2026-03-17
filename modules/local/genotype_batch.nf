#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GENOTYPEBATCH {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(cohort), path(rf_cutoffs), path(pe_file), path(rd_file), path(sr_file), path(median_coverage), path(ploidy_table), path(vcf)

    output:
    tuple val(batch_name), path("${batch_name}/${batch_name}.genotype_batch.depth.vcf.gz"), emit: genotyped_depth_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.genotype_batch.depth.vcf.gz.tbi"), emit: genotyped_depth_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.genotype_batch.pesr.vcf.gz"), emit: genotyped_pesr_vcf
    tuple val(batch_name), path("${batch_name}/${batch_name}.genotype_batch.pesr.vcf.gz.tbi"), emit: genotyped_pesr_vcf_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.rd_geno_params.tsv"), emit: genotyping_rd_table
    tuple val(batch_name), path("${batch_name}/${batch_name}.pe_geno_params.tsv"), emit: genotyping_pe_table
    tuple val(batch_name), path("${batch_name}/${batch_name}.sr_geno_params.tsv"), emit: genotyping_sr_table
    tuple val(batch_name), path("${batch_name}/${batch_name}.regeno_coverage_medians.tsv.gz"), emit: regeno_coverage_medians
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.genotypebatch_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.genotype_batch ?: [:])

    def batch_id = batch_name?.toString()?.trim()
    if (!batch_id) {
        throw new IllegalArgumentException("GenotypeBatch batch is required")
    }
    def cohort_id = cohort?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("GenotypeBatch cohort is required")
    }

    def dynamic = [
        "GenotypeBatch.batch"           : batch_id,
        "GenotypeBatch.rf_cutoffs"      : rf_cutoffs.toRealPath().toString(),
        "GenotypeBatch.pe_file"         : pe_file.toRealPath().toString(),
        "GenotypeBatch.rd_file"         : rd_file.toRealPath().toString(),
        "GenotypeBatch.sr_file"         : sr_file.toRealPath().toString(),
        "GenotypeBatch.median_coverage" : median_coverage.toRealPath().toString(),
        "GenotypeBatch.vcf"             : vcf.toRealPath().toString(),
        "GenotypeBatch.ploidy_table"    : ploidy_table.toRealPath().toString()
    ]
    file("genotype_batch_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV GENOTYPEBATCH] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out genotype_batch_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file genotype_batch_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.genotypebatch_wdl} \\
        -i genotype_batch_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p "${batch_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${batch_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/GenotypeBatch -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected GenotypeBatch output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${batch_id}.genotype_batch.depth.vcf.gz" 1
    copy_outputs "${batch_id}.genotype_batch.depth.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.genotype_batch.pesr.vcf.gz" 1
    copy_outputs "${batch_id}.genotype_batch.pesr.vcf.gz.tbi" 1
    copy_outputs "${batch_id}.rd_geno_params.tsv" 1
    copy_outputs "${batch_id}.pe_geno_params.tsv" 1
    copy_outputs "${batch_id}.sr_geno_params.tsv" 1
    copy_outputs "${batch_id}.regeno_coverage_medians.tsv.gz" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${batch_name}
    touch ${batch_name}/${batch_name}.genotype_batch.depth.vcf.gz
    touch ${batch_name}/${batch_name}.genotype_batch.depth.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.genotype_batch.pesr.vcf.gz
    touch ${batch_name}/${batch_name}.genotype_batch.pesr.vcf.gz.tbi
    touch ${batch_name}/${batch_name}.rd_geno_params.tsv
    touch ${batch_name}/${batch_name}.pe_geno_params.tsv
    touch ${batch_name}/${batch_name}.sr_geno_params.tsv
    touch ${batch_name}/${batch_name}.regeno_coverage_medians.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
