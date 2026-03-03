process GATK_UPDATEVCFSEQUENCEDICTIONARY {
    tag "${meta.id}"
    label 'process_low'

    container "${workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e3d753d93f57969fe76b8628a8dfcd23ef44bccd08c4ced7089c1f94bf47c89f/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel_htslib_samtools:d3becb6465454c35'}"

    input:
    tuple val(meta), path(vcf)
    path(dict)

    output:
    tuple val(meta), path("${prefix}.vcf"), emit: vcf
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def base = vcf.baseName.replaceAll(/\.vcf$/, "")
    prefix   = task.ext.prefix ?: "${base}_updated"
    """
    gatk UpdateVCFSequenceDictionary \\
        --variant ${vcf} \\
        --source-dictionary ${dict} \\
        --replace true \\
        --output ${prefix}.vcf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/.* //')
    END_VERSIONS
    """
}
