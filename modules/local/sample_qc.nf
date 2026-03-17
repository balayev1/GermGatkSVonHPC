#!/usr/bin/env nextflow

process SAMPLE_QC {
    tag "${cohort}"
    label 'process_single'

    container 'docker://rocker/tidyverse:4.5.3'

    input:
    tuple val(cohort), val(sample_ids), path(insert_size_metrics), path(ped_file), path(evidence_qc_table), path(sample_sex_assignments)

    output:
    tuple val(cohort), path("*.pdf"), emit: sample_qc_plots
    tuple val(cohort), path("*.tsv"), emit: sample_qc_reports
    tuple val(cohort), val(sample_ids), path("passing_samples_metadata.tsv"), emit: passing_samples_metadata
    tuple val(cohort), path("*.ped"), emit: updated_ped
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

    mkdir -p sample_qc_inputs

    evidence_qc_merged="sample_qc_inputs/evidence_qc_table.tsv"
    first=1
    for f in ${evidence_qc_table}; do
        if [[ "\$first" -eq 1 ]]; then
            cat "\$f" > "\$evidence_qc_merged"
            first=0
        else
            tail -n +2 "\$f" >> "\$evidence_qc_merged"
        fi
    done

    sample_sex_sources="sample_qc_inputs/sample_sex_sources.list"
    : > "\$sample_sex_sources"
    for f in ${sample_sex_assignments}; do
        printf "%s\\n" "\$f" >> "\$sample_sex_sources"
    done

    mkdir -p insert_size_files
    for f in ${insert_size_metrics}; do
        sample_id=\$(basename "\$f")
        sample_id=\${sample_id%%.*}
        cp "\$f" "insert_size_files/\${sample_id}.insert_size_metrics.txt"
    done

    Rscript sample_qc.R \\
        sample_qc_metadata.tsv \\
        "\$evidence_qc_merged" \\
        "\$sample_sex_sources" \\
        insert_size_files \\
        ${num_samples} \\
        "${ped_file}" \\
        "${cohort}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript --version 2>&1 | sed 's/^Rscript (R) version //')
    END_VERSIONS
    """

    stub:
    """
    touch WGD_Score_Distribution.pdf
    touch Excluded_Samples_Report.tsv
    touch Excluded_Sample_ID_only.tsv
    touch passing_samples_metadata.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: "stub"
    END_VERSIONS
    """
}
