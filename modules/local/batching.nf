#!/usr/bin/env nextflow

process BATCHING {
    tag "${cohort}"
    label 'process_single'

    input:
    tuple val(cohort), path(evidence_qc_results), path(ped_file)

    output:
    tuple val(cohort), path("batching_out"), emit: batching_out
    path "versions.yml", emit: versions

    script:
    if (!params.batching_include_metrics) {
        throw new IllegalArgumentException("Set --batching_include_metrics when --run_batching true (comma-separated list).")
    }
    def include_bins = params.batching_include_bins ?: ''
    def batch_prefix = params.batching_batch_prefix ?: "${cohort}_batch_"
    def batch_suffix = params.batching_batch_suffix ?: ''
    """
    pass_metadata=\$(find ${evidence_qc_results} -type f -name "passing_samples_metadata.tsv" | head -n 1)
    if [[ -z "\${pass_metadata}" ]]; then
        echo "ERROR: passing_samples_metadata.tsv not found under ${evidence_qc_results}" >&2
        exit 1
    fi

    python ${projectDir}/bin/Batching.py \\
        --pass-metadata "\${pass_metadata}" \\
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p batching_out/batching
    touch batching_out/batching/batch_assignments.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}

