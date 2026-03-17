#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_RESOLVECOMPLEXVARIANTS {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), path(cluster_vcfs), path(cluster_bothside_pass_lists), path(cluster_background_fail_lists), path(disc_files), path(rf_cutoff_files)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.reshard_vcf.*.resharded.vcf.gz"), emit: complex_resolve_vcfs
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.reshard_vcf.*.resharded.vcf.gz.tbi"), emit: complex_resolve_vcf_indexes
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.all.sr_bothside_pass.updated3.txt"), emit: complex_resolve_bothside_pass_list
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.all.sr_background_fail.updated3.txt"), emit: complex_resolve_background_fail_list
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.*.breakpoint_overlap.dropped_records.vcf.gz"), emit: breakpoint_overlap_dropped_record_vcfs
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.*.breakpoint_overlap.dropped_records.vcf.gz.tbi"), emit: breakpoint_overlap_dropped_record_vcf_indexes
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.complex_resolve.vcf.gz"), emit: complex_resolve_merged_vcf, optional: true
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.complex_resolve.vcf.gz.tbi"), emit: complex_resolve_merged_vcf_index, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.resolvecomplexvariants_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.resolve_complex_variants ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("ResolveComplexVariants cohort is required")
    }

    def asList = { value ->
        if (value == null) {
            return []
        }
        value instanceof List ? value : [value]
    }
    def toRealPathList = { value -> asList(value).collect { it.toRealPath().toString() } }

    def clusterVcfs = toRealPathList(cluster_vcfs)
    def bothsidePassLists = toRealPathList(cluster_bothside_pass_lists)
    def backgroundFailLists = toRealPathList(cluster_background_fail_lists)
    def discFiles = toRealPathList(disc_files)
    def rfCutoffFiles = toRealPathList(rf_cutoff_files)

    if (!clusterVcfs || !bothsidePassLists || !backgroundFailLists || !discFiles || !rfCutoffFiles) {
        throw new IllegalArgumentException("ResolveComplexVariants requires non-empty cluster VCF, SR list, disc file, and cutoff file inputs")
    }

    def dynamic = [
        "ResolveComplexVariants.cohort_name"                  : cohort_id,
        "ResolveComplexVariants.cluster_vcfs"                : clusterVcfs,
        "ResolveComplexVariants.cluster_bothside_pass_lists" : bothsidePassLists,
        "ResolveComplexVariants.cluster_background_fail_lists": backgroundFailLists,
        "ResolveComplexVariants.disc_files"                  : discFiles,
        "ResolveComplexVariants.rf_cutoff_files"             : rfCutoffFiles
    ]
    file("resolve_complex_variants_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV RESOLVECOMPLEXVARIANTS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \
        --template ${template_path} \
        --out resolve_complex_variants_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file resolve_complex_variants_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.resolvecomplexvariants_wdl} \
        -i resolve_complex_variants_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/ResolveComplexVariants -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected ResolveComplexVariants output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.reshard_vcf.*.resharded.vcf.gz" 1
    copy_outputs "${cohort_id}.reshard_vcf.*.resharded.vcf.gz.tbi" 1
    copy_outputs "${cohort_id}.all.sr_bothside_pass.updated3.txt" 1
    copy_outputs "${cohort_id}.all.sr_background_fail.updated3.txt" 1
    copy_outputs "${cohort_id}.*.breakpoint_overlap.dropped_records.vcf.gz" 1
    copy_outputs "${cohort_id}.*.breakpoint_overlap.dropped_records.vcf.gz.tbi" 1
    copy_outputs "${cohort_id}.complex_resolve.vcf.gz" 0
    copy_outputs "${cohort_id}.complex_resolve.vcf.gz.tbi" 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.reshard_vcf.shard_1.resharded.vcf.gz
    touch ${cohort_name}/${cohort_name}.reshard_vcf.shard_1.resharded.vcf.gz.tbi
    touch ${cohort_name}/${cohort_name}.all.sr_bothside_pass.updated3.txt
    touch ${cohort_name}/${cohort_name}.all.sr_background_fail.updated3.txt
    touch ${cohort_name}/${cohort_name}.chr1.breakpoint_overlap.dropped_records.vcf.gz
    touch ${cohort_name}/${cohort_name}.chr1.breakpoint_overlap.dropped_records.vcf.gz.tbi
    touch ${cohort_name}/${cohort_name}.complex_resolve.vcf.gz
    touch ${cohort_name}/${cohort_name}.complex_resolve.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
