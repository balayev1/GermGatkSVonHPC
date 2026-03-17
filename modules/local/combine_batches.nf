#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_COMBINEBATCHES {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), val(batches), path(ped_file), path(pesr_vcfs), path(depth_vcfs)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.combine_batches.*.svtk_formatted.vcf.gz"), emit: combined_vcfs
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.combine_batches.*.svtk_formatted.vcf.gz.tbi"), emit: combined_vcf_indexes
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.combine_batches.*.high_sr_background.txt"), emit: cluster_background_fail_lists
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.combine_batches.*.bothsides_sr_support.txt"), emit: cluster_bothside_pass_lists
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.combine_batches.concat_all_contigs.vcf.gz"), emit: combine_batches_merged_vcf, optional: true
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.combine_batches.concat_all_contigs.vcf.gz.tbi"), emit: combine_batches_merged_vcf_index, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.combinebatches_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.combine_batches ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    
    def asList = { value ->
        if (value == null) {
            return []
        }
        value instanceof List ? value : [value]
    }
    def toRealPathList = { value -> asList(value).collect { it.toRealPath().toString() } }

    def dynamic = [
        "CombineBatches.cohort_name" : cohort_id,
        "CombineBatches.batches"     : asList(batches).collect { it.toString() },
        "CombineBatches.ped_file"    : ped_file.toRealPath().toString(),
        "CombineBatches.pesr_vcfs"   : toRealPathList(pesr_vcfs),
        "CombineBatches.depth_vcfs"  : toRealPathList(depth_vcfs)
    ]
    file("combine_batches_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV COMBINEBATCHES] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \
        --template ${template_path} \
        --out combine_batches_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file combine_batches_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.combinebatches_wdl} \
        -i combine_batches_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/CombineBatches -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected CombineBatches output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.combine_batches.*.svtk_formatted.vcf.gz" 1
    copy_outputs "${cohort_id}.combine_batches.*.svtk_formatted.vcf.gz.tbi" 1
    copy_outputs "${cohort_id}.combine_batches.*.high_sr_background.txt" 1
    copy_outputs "${cohort_id}.combine_batches.*.bothsides_sr_support.txt" 1
    copy_outputs "${cohort_id}.combine_batches.concat_all_contigs.vcf.gz" 0
    copy_outputs "${cohort_id}.combine_batches.concat_all_contigs.vcf.gz.tbi" 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.combine_batches.chr1.svtk_formatted.vcf.gz
    touch ${cohort_name}/${cohort_name}.combine_batches.chr1.svtk_formatted.vcf.gz.tbi
    touch ${cohort_name}/${cohort_name}.combine_batches.chr1.high_sr_background.txt
    touch ${cohort_name}/${cohort_name}.combine_batches.chr1.bothsides_sr_support.txt
    touch ${cohort_name}/${cohort_name}.combine_batches.concat_all_contigs.vcf.gz
    touch ${cohort_name}/${cohort_name}.combine_batches.concat_all_contigs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
