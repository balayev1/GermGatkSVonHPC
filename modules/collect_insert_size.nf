process COLLECT_INSERT_SIZE {
    tag "${meta.id}"

    container: 'docker://broadinstitute/picard:latest'
    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.insert_size_metrics.txt"), emit: metrics
    tuple val(meta), path("${meta.id}.insert_size_histogram.pdf"), emit: histogram

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def min_pct_args = args.min_pct ? "${args.min_pct}" : "0.05"
    def ref_genome_args = args.ref_genome ? "${args.ref_genome}" : "null"

    """    
    picard CollectInsertSizeMetrics \\
        I=${bam} \\
        O=${meta.id}.insert_size_metrics.txt \\
        H=${meta.id}.insert_size_histogram.pdf \\
        M=${min_pct_args} \\
        R=${ref_genome_args}
    """
}