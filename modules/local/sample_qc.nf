#!/usr/bin/env nextflow

process SAMPLE_QC {
    tag "${cohort}"
    label 'process_single'

    container 'docker://rocker/tidyverse:latest'

    input:
    tuple val(cohort), val(sample_ids), path(insert_size_metrics), path(ped_file), path(evidence_qc_results)

    output:
    path "*.pdf", emit: sample_qc_plots
    path "*.tsv", emit: sample_qc_reports
    path "versions.yml", emit: versions

    script:
    def ids = (sample_ids instanceof List ? sample_ids : [sample_ids]).collect { it.toString() }
    if (!ids) {
        throw new IllegalArgumentException("SAMPLE_QC received an empty sample list for cohort '${cohort}'")
    }
    def num_samples = ids.size()
    file("sample_ids.list").text = ids.join('\n') + '\n'

    """
    awk 'NR==FNR { keep[\$1]=1; next } (\$2 in keep) { sex=(\$5=="1" ? "XY" : (\$5=="2" ? "XX" : "XX")); print \$2 "\\t" sex }' \\
        sample_ids.list ${ped_file} > sample_qc_metadata.tsv

    mkdir -p insert_size_files
    for f in ${insert_size_metrics}; do
        sample_id=\$(basename "\$f")
        sample_id=\${sample_id%%.*}
        cp "\$f" "insert_size_files/\${sample_id}.insert_size_metrics.txt"
    done

    Rscript sample_qc.R \\
        sample_qc_metadata.tsv \\
        ${evidence_qc_results} \\
        insert_size_files \\
        ${num_samples}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript --version 2>&1 | sed 's/^Rscript (R) version //')
    END_VERSIONS
    """

    stub:
    """
    touch WGD_Score_Distribution.pdf
    touch Excluded_Samples_Report.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: "stub"
    END_VERSIONS
    """
}
