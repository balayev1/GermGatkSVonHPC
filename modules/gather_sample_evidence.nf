#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GATHER_SAMPLE_EVIDENCE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}"), emit: full_evidence_dir

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}

    # Create sample JSON
    sed -e "s|PLACEHOLDER_ID|${prefix}|g" \
        -e "s|PLACEHOLDER_BAM|${bam.toRealPath()}|g" \
        -e "s|PLACEHOLDER_BAI|${bai.toRealPath()}|g" \
        ${params.gse_template} > ${prefix}_inputs.json

    # Run GSE using Cromwell
    java -Xmx12G -Dconfig.file=${params.cromwell_conf} -jar ${params.cromwell_jar} \
        run ${params.gse_wdl} \
        -i ${prefix}_inputs.json \
        -p ${params.deps_zip}

    # Move results into the output directory
    mv cromwell-executions/GatherSampleEvidence/*/call-* ${prefix}/
    find ${prefix}/ -type d -name 'tmpVcfs' -exec rm -rf {} +
    """
}