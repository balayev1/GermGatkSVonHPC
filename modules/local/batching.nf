#!/usr/bin/env nextflow

process BATCHING {
    tag "${cohort}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity'
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e3d753d93f57969fe76b8628a8dfcd23ef44bccd08c4ced7089c1f94bf47c89f/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel_htslib_samtools:d3becb6465454c35'}"


    input:
    tuple val(cohort), path(passing_samples_metadata), path(ped_file)

    output:
    tuple val(cohort), path("batching_out/batch_assignments.tsv"), emit: batch_assignments
    tuple val(cohort), path("batching_out/batching_metadata.tsv"), emit: batching_metadata
    tuple val(cohort), path("batching_out/metric_plots/*.png"), emit: metric_plots
    tuple val(cohort), path("batching_out/cluster_plots/*.png"), emit: cluster_plots
    tuple val(cohort), path("batching_out/pannelled_cluster_plots/*.png"), emit: panelled_cluster_plots
    path "versions.yml", emit: versions

    script:
    if (!params.batching_include_metrics) {
        throw new IllegalArgumentException("Set --batching_include_metrics when --run_batching true (comma-separated list).")
    }
    def include_bins = params.batching_include_bins ?: ''
    def batch_prefix = params.batching_batch_prefix ?: "${cohort}_batch_"
    def batch_suffix = params.batching_batch_suffix ?: ''
    """
    mkdir -p batching_out/input

    merged_passing_samples_metadata="batching_out/input/passing_samples_metadata.tsv"
    first=1
    for f in ${passing_samples_metadata}; do
        if [[ "\$first" -eq 1 ]]; then
            cat "\$f" > "\$merged_passing_samples_metadata"
            first=0
        else
            tail -n +2 "\$f" >> "\$merged_passing_samples_metadata"
        fi
    done

    python ${projectDir}/bin/batching.py \\
        --pass-metadata "\$merged_passing_samples_metadata" \\
        --include-metrics "${params.batching_include_metrics}" \\
        --include-bins "${include_bins}" \\
        --target-batch-size ${params.batching_target_batch_size} \\
        --min-batch-size ${params.batching_min_batch_size} \\
        --max-batch-size ${params.batching_max_batch_size} \\
        --batch-prefix "${batch_prefix}" \\
        --batch-suffix "${batch_suffix}" \\
        --ped-file ${ped_file} \\
        --outdir batching_out \\
        --plot-prefix "${cohort}"

    if [[ -d batching_out/batching ]]; then
        shopt -s dotglob nullglob
        mv batching_out/batching/* batching_out/
        rmdir batching_out/batching || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p batching_out/metric_plots
    mkdir -p batching_out/cluster_plots
    mkdir -p batching_out/pannelled_cluster_plots
    touch batching_out/batch_assignments.tsv
    touch batching_out/batching_metadata.tsv
    touch batching_out/metric_plots/${cohort}_distribution_1_batches.png
    touch batching_out/cluster_plots/${cohort}_distribution_1_batches.png
    touch batching_out/pannelled_cluster_plots/${cohort}_distribution_1_batches.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
