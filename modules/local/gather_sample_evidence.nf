#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GATHERSAMPLEEVIDENCE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path(sd_locs_vcf)

    output:
    tuple val(meta), path("${meta.id}/call-CollectCounts/**/${meta.id}.counts.tsv.gz"), emit: counts
    tuple val(meta), path("${meta.id}/call-CollectSVEvidence/**/${meta.id}.pe.txt.gz"), emit: pe_file
    tuple val(meta), path("${meta.id}/call-CollectSVEvidence/**/${meta.id}.sr.txt.gz"), emit: sr_file
    tuple val(meta), path("${meta.id}/call-CollectSVEvidence/**/${meta.id}.sd.txt.gz"), emit: sd_file
    tuple val(meta), path("${meta.id}/call-Manta/**/${meta.id}.manta.vcf.gz"), emit: manta_vcf
    tuple val(meta), path("${meta.id}/call-Whamg/**/${meta.id}.wham.vcf.gz"), emit: wham_vcf
    tuple val(meta), path("${meta.id}/call-Scramble/**/${meta.id}.scramble.vcf.gz"), emit: scramble_vcf
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_json = "${prefix}_inputs.json"
    def static_json = JsonOutput.toJson(params.tool_inputs?.gse ?: [:])
    def bam_path = bam.toRealPath().toString()
    def bai_path = bai.toRealPath().toString()
    def sd_locs_vcf_path = sd_locs_vcf.toRealPath().toString()

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
        --set "GatherSampleEvidence.bam_or_cram_index=${bai_path}" \
        --set "GatherSampleEvidence.sd_locs_vcf=${sd_locs_vcf_path}"

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

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
