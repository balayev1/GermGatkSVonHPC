#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GENOTYPECOMPLEXVARIANTS {
    tag "${cohort_name}"
    label 'process_medium'

    input:
    tuple val(cohort_name), val(batches), path(depth_vcfs), path(complex_resolve_vcfs), path(complex_resolve_vcf_indexes), path(ped_file), path(bincov_files), path(genotyping_rd_tables), path(median_coverage_files)

    output:
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.*.regenotyped.vcf.gz"), emit: complex_genotype_vcfs
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.*.regenotyped.vcf.gz.tbi"), emit: complex_genotype_vcf_indexes
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.complex_genotype.vcf.gz"), emit: complex_genotype_merged_vcf, optional: true
    tuple val(cohort_name), path("${cohort_name}/${cohort_name}.complex_genotype.vcf.gz.tbi"), emit: complex_genotype_merged_vcf_index, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.genotypecomplexvariants_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.genotype_complex_variants ?: [:])

    def cohort_id = cohort_name?.toString()?.trim()
    if (!cohort_id) {
        throw new IllegalArgumentException("GenotypeComplexVariants cohort is required")
    }

    def asList = { value ->
        if (value == null) {
            return []
        }
        value instanceof List ? value : [value]
    }
    def toRealPathList = { value -> asList(value).collect { it.toRealPath().toString() } }

    def batchIds = asList(batches).collect { it.toString() }
    def depthVcfs = toRealPathList(depth_vcfs)
    def complexResolveVcfs = toRealPathList(complex_resolve_vcfs)
    def complexResolveVcfIndexes = toRealPathList(complex_resolve_vcf_indexes)
    def bincovFiles = toRealPathList(bincov_files)
    def genotypingRdTables = toRealPathList(genotyping_rd_tables)
    def medianCoverageFiles = toRealPathList(median_coverage_files)

    if (!batchIds || !depthVcfs || !complexResolveVcfs || !complexResolveVcfIndexes || !bincovFiles || !genotypingRdTables || !medianCoverageFiles) {
        throw new IllegalArgumentException("GenotypeComplexVariants requires non-empty batch, depth VCF, complex-resolve VCF, bincov, RD table, and median coverage inputs")
    }

    def dynamic = [
        "GenotypeComplexVariants.cohort_name"                 : cohort_id,
        "GenotypeComplexVariants.batches"                     : batchIds,
        "GenotypeComplexVariants.depth_vcfs"                  : depthVcfs,
        "GenotypeComplexVariants.complex_resolve_vcfs"        : complexResolveVcfs,
        "GenotypeComplexVariants.complex_resolve_vcf_indexes" : complexResolveVcfIndexes,
        "GenotypeComplexVariants.ped_file"                    : ped_file.toRealPath().toString(),
        "GenotypeComplexVariants.bincov_files"                : bincovFiles,
        "GenotypeComplexVariants.genotyping_rd_tables"        : genotypingRdTables,
        "GenotypeComplexVariants.median_coverage_files"       : medianCoverageFiles
    ]
    file("genotype_complex_variants_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV GENOTYPECOMPLEXVARIANTS] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \
        --template ${template_path} \
        --out genotype_complex_variants_inputs.json \
        --static-json '${static_json}' \
        --merge-json-file genotype_complex_variants_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \
        run ${params.genotypecomplexvariants_wdl} \
        -i genotype_complex_variants_inputs.json \
        -p ${params.deps_zip}

    mkdir -p "${cohort_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${cohort_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/GenotypeComplexVariants -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected GenotypeComplexVariants output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${cohort_id}.*.regenotyped.vcf.gz" 1
    copy_outputs "${cohort_id}.*.regenotyped.vcf.gz.tbi" 1
    copy_outputs "${cohort_id}.complex_genotype.vcf.gz" 0
    copy_outputs "${cohort_id}.complex_genotype.vcf.gz.tbi" 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*"//; s/".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort_name}
    touch ${cohort_name}/${cohort_name}.chr1.regenotyped.vcf.gz
    touch ${cohort_name}/${cohort_name}.chr1.regenotyped.vcf.gz.tbi
    touch ${cohort_name}/${cohort_name}.complex_genotype.vcf.gz
    touch ${cohort_name}/${cohort_name}.complex_genotype.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
