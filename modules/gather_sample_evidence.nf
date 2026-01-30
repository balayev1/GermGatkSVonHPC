#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process GATHER_SAMPLE_EVIDENCE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}"), emit: full_evidence_dir
    // tuple val(meta), path("*.counts.tsv.gz"), emit: read_counts
    // tuple val(meta), path("*.pe.txt.gz"), emit: discordant_pairs
    // tuple val(meta), path("*.sr.txt.gz"), emit: split_reads
    // tuple val(meta), path("*.sd.txt.gz"), emit: site_depth
    // tuple val(meta), path("*.manta.std.vcf.gz"), emit: manta_vcf
    // tuple val(meta), path("*.wham.std.vcf.gz"), emit: wham_vcf
    // tuple val(meta), path("*.scramble.vcf.gz"), emit: scramble_vcf

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """    
    mkdir -p ${prefix}

    # Create sample JSON
    sed -e "s|PLACEHOLDER_ID|${prefix}|g" \
        -e "s|PLACEHOLDER_BAM|${bam.toAbsolutePath()}|g" \
        -e "s|PLACEHOLDER_BAI|${bai.toAbsolutePath()}|g" \
        ${params.gse_template} > ${prefix}_inputs.json
    
    # Run GSE using Cromwell
    java -Xmx12G -Dconfig.file=${params.cromwell_conf} -jar ${params.cromwell_jar} \
        run ${params.gse_wdl} \
        -i ${prefix}_inputs.json \
        -p ${params.deps_zip}

    # Move results to the specific output folder
    mv cromwell-executions/GatherSampleEvidence/*/call-* ${prefix}/ && find ${prefix}/ -type d -name 'tmpVcfs' -exec rm -rf {} +
    """
}