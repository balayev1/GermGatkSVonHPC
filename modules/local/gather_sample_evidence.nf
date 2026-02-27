#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GATHERSAMPLEEVIDENCE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}/call-*"), emit: gse_outdir
    tuple val(meta), path("${meta.id}/**/call-CountsMetrics/**/${meta.id}.counts.tsv.gz"), path("${meta.id}/**/execution/**/${meta.id}.manta.std.vcf.gz"), path("${meta.id}/**/execution/**/${meta.id}.wham.std.vcf.gz"), path("${meta.id}/**/execution/**/${meta.id}.scramble.vcf.gz"), emit: gse_outfiles
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_json = "${prefix}_inputs.json"
    def static_json = JsonOutput.toJson(params.tool_inputs?.gse ?: [:])
    def bam_path = bam.toRealPath().toString()
    def bai_path = bai.toRealPath().toString()

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV GATHERSAMPLEEVIDENCE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    mkdir -p ${prefix}

    render_json_template.py \
        --template ${params.gse_template} \
        --out ${out_json} \
        --static-json '${static_json}' \
        --set "GatherSampleEvidence.sample_id=${prefix}" \
        --set "GatherSampleEvidence.bam_or_cram_file=${bam_path}" \
        --set "GatherSampleEvidence.bam_or_cram_index=${bai_path}"

    # Run GSE using Cromwell
    java -Xmx${avail_mem}M \\
        -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.gse_wdl} \\
        -i ${out_json} \\
        -p ${params.deps_zip}

    # Move results into the output directory
    mv cromwell-executions/GatherSampleEvidence/*/call-* ${prefix}/
    find ${prefix}/ -type d -name 'tmpVcfs' -exec rm -rf {} +

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/call-stub
    touch ${prefix}/call-stub/.stub

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
