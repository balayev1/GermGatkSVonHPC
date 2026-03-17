#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_REFINECOMPLEXVARIANTS {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(vcf), val(batch_name_list), path(batch_sample_lists), path(pe_metrics), path(pe_metrics_indexes), path(depth_del_beds), path(depth_dup_beds)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.cpx_refined.vcf.gz"), emit: cpx_refined_vcf
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.cpx_refined.vcf.gz.tbi"), emit: cpx_refined_vcf_index
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.CPX_evidence.txt"), emit: cpx_evidences
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.refinecomplexvariants_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.refine_complex_variants ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("RefineComplexVariants prefix is required")
    }

    def asList = { value ->
        if (value == null) {
            return []
        }
        value instanceof List ? value : [value]
    }
    def toRealPathList = { value -> asList(value).collect { it.toRealPath().toString() } }

    def batchNames = asList(batch_name_list).collect { it.toString() }
    def batchSampleLists = toRealPathList(batch_sample_lists)
    def peMetrics = toRealPathList(pe_metrics)
    def peMetricsIndexes = toRealPathList(pe_metrics_indexes)
    def depthDelBeds = toRealPathList(depth_del_beds)
    def depthDupBeds = toRealPathList(depth_dup_beds)

    if (!batchNames || !batchSampleLists || !peMetrics || !peMetricsIndexes || !depthDelBeds || !depthDupBeds) {
        throw new IllegalArgumentException("RefineComplexVariants requires non-empty batch names, sample lists, PE metrics, PE metric indexes, and depth support BED inputs")
    }

    def dynamic = [
        "RefineComplexVariants.vcf"              : vcf.toRealPath().toString(),
        "RefineComplexVariants.prefix"           : cohort_id,
        "RefineComplexVariants.batch_name_list"  : batchNames,
        "RefineComplexVariants.batch_sample_lists": batchSampleLists,
        "RefineComplexVariants.PE_metrics"       : peMetrics,
        "RefineComplexVariants.PE_metrics_indexes": peMetricsIndexes,
        "RefineComplexVariants.Depth_DEL_beds"   : depthDelBeds,
        "RefineComplexVariants.Depth_DUP_beds"   : depthDupBeds
    ]
    file("refine_complex_variants_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV REFINECOMPLEXVARIANTS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \
        --template ${template_path} \
        --out refine_complex_variants_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file refine_complex_variants_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.refinecomplexvariants_wdl} \
        -i refine_complex_variants_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/RefineComplexVariants -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected RefineComplexVariants output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.cpx_refined.vcf.gz" 1
    copy_outputs "${cohort_id}.cpx_refined.vcf.gz.tbi" 1
    copy_outputs "${cohort_id}.CPX_evidence.txt" 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.cpx_refined.vcf.gz
    touch ${cohort_name}/${cohort_name}.cpx_refined.vcf.gz.tbi
    touch ${cohort_name}/${cohort_name}.CPX_evidence.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
