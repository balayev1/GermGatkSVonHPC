
process COLLECT_INSERT_SIZE {
    tag "${meta.id}"

    container 'docker://broadinstitute/picard:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.insert_size_metrics.txt"), emit: metrics
    tuple val(meta), path("${meta.id}.insert_size_histogram.pdf"), emit: histogram

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: [:]
    def prefix = task.ext.prefix ?: "${meta.id}"

    def min_pct = args.min_pct ?: "0.05"
    def ref_arg = args.ref_genome ? "R=${args.ref_genome}" : ""

    """
    java -Xmx6G -jar /usr/picard/picard.jar CollectInsertSizeMetrics \\
        I=${bam} \\
        O=${prefix}.insert_size_metrics.txt \\
        H=${prefix}.insert_size_histogram.pdf \\
        M=${min_pct} \\
        ${ref_arg}
    """
}