#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_EVIDENCEQC {
    tag "${cohort}"
    label 'process_medium'

    input:
    tuple val(cohort), val(sample_ids), path(count_files), path(manta_vcf_files), path(wham_vcf_files), path(scramble_vcf_files)

    output:
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.evidence_qc_table.tsv"), emit: qc_table
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.RD.txt.gz"), emit: bincov_matrix
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.RD.txt.gz.tbi"), emit: bincov_matrix_index
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}_medianCov.transposed.bed"), emit: bincov_median
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}_ploidy_matrix.bed.gz"), emit: ploidy_matrix
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}_ploidy_est.tar.gz"), emit: ploidy_plots
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}_WGD_scoring_matrix_output.bed.gz"), emit: WGD_matrix
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}_WGD_scores.txt.gz"), emit: WGD_scores
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}_WGD_score_distributions.pdf"), emit: WGD_dist
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.manta.variant_counts.tsv"), emit: manta_variant_counts 
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.wham.variant_counts.tsv"), emit: wham_variant_counts 
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.scramble.variant_counts.tsv"), emit: scramble_variant_counts    
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.manta.QC.outlier.low"), emit: manta_qc_low    
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.manta.QC.outlier.high"), emit: manta_qc_high    
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.scramble.QC.outlier.low"), emit: scramble_qc_low    
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.scramble.QC.outlier.high"), emit: scramble_qc_high    
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.wham.QC.outlier.low"), emit: wham_qc_low    
    tuple val(cohort), val(sample_ids), path("${cohort}/${cohort}.wham.QC.outlier.high"), emit: wham_qc_high    
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.evidqc_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.evidence_qc ?: [:])
    def cohort_name = cohort?.toString()?.trim()
    if (!cohort_name) {
        throw new IllegalArgumentException("EvidenceQC cohort is required and must come from samplesheet metadata (meta.cohort)")
    }
    def asList = { value -> value instanceof List ? value : [value] }
    def samples = asList(sample_ids).collect { it.toString() }
    def counts = asList(count_files).collect { it.toRealPath().toString() }
    def manta = asList(manta_vcf_files).collect { it.toRealPath().toString() }
    def wham = asList(wham_vcf_files).collect { it.toRealPath().toString() }
    def scramble = asList(scramble_vcf_files).collect { it.toRealPath().toString() }

    if (!samples) {
        throw new IllegalArgumentException("EvidenceQC received an empty sample list")
    }
    def n = samples.size()
    if (counts.size() != n || manta.size() != n || wham.size() != n || scramble.size() != n) {
        throw new IllegalArgumentException("EvidenceQC input size mismatch: samples=${n}, counts=${counts.size()}, manta=${manta.size()}, wham=${wham.size()}, scramble=${scramble.size()}")
    }

    def dynamic = [
        "EvidenceQC.batch"        : cohort_name,
        "EvidenceQC.samples"      : samples,
        "EvidenceQC.counts"       : counts,
        "EvidenceQC.manta_vcfs"   : manta,
        "EvidenceQC.wham_vcfs"    : wham,
        "EvidenceQC.scramble_vcfs": scramble
    ]
    file("evidence_qc_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV EVIDENCEQC] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    
    """
    render_json_template.py \\
        --template ${template_path} \\
        --out evidence_qc_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file evidence_qc_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL
    
    # Execute Cromwell
    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.evidqc_wdl} \\
        -i evidence_qc_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p "${cohort_name}"

    copy_output_file() {
        local filename="\$1"
        local source
        source=\$(find cromwell-executions/EvidenceQC -type f -name "\${filename}" | head -n 1 || true)
        if [[ -z "\$source" ]]; then
            echo "ERROR: Expected EvidenceQC output not found: \${filename}" >&2
            exit 1
        fi
        cp -L "\$source" "${cohort_name}/\${filename}"
    }

    copy_output_file "${cohort_name}.evidence_qc_table.tsv"
    copy_output_file "${cohort_name}.RD.txt.gz"
    copy_output_file "${cohort_name}.RD.txt.gz.tbi"
    copy_output_file "${cohort_name}_medianCov.transposed.bed"
    copy_output_file "${cohort_name}_ploidy_matrix.bed.gz"
    copy_output_file "${cohort_name}_ploidy_est.tar.gz"
    copy_output_file "${cohort_name}_WGD_scoring_matrix_output.bed.gz"
    copy_output_file "${cohort_name}_WGD_scores.txt.gz"
    copy_output_file "${cohort_name}_WGD_score_distributions.pdf"
    copy_output_file "${cohort_name}.manta.variant_counts.tsv"
    copy_output_file "${cohort_name}.wham.variant_counts.tsv"
    copy_output_file "${cohort_name}.scramble.variant_counts.tsv"
    copy_output_file "${cohort_name}.manta.QC.outlier.low"
    copy_output_file "${cohort_name}.manta.QC.outlier.high"
    copy_output_file "${cohort_name}.scramble.QC.outlier.low"
    copy_output_file "${cohort_name}.scramble.QC.outlier.high"
    copy_output_file "${cohort_name}.wham.QC.outlier.low"
    copy_output_file "${cohort_name}.wham.QC.outlier.high"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${cohort}
    touch ${cohort}/${cohort}.evidence_qc_table.tsv
    touch ${cohort}/${cohort}.RD.txt.gz
    touch ${cohort}/${cohort}.RD.txt.gz.tbi
    touch ${cohort}/${cohort}_medianCov.transposed.bed
    touch ${cohort}/${cohort}_ploidy_matrix.bed.gz
    touch ${cohort}/${cohort}_ploidy_est.tar.gz
    touch ${cohort}/${cohort}_WGD_scoring_matrix_output.bed.gz
    touch ${cohort}/${cohort}_WGD_scores.txt.gz
    touch ${cohort}/${cohort}_WGD_score_distributions.pdf
    touch ${cohort}/${cohort}.manta.variant_counts.tsv
    touch ${cohort}/${cohort}.wham.variant_counts.tsv
    touch ${cohort}/${cohort}.scramble.variant_counts.tsv
    touch ${cohort}/${cohort}.manta.QC.outlier.low
    touch ${cohort}/${cohort}.manta.QC.outlier.high
    touch ${cohort}/${cohort}.scramble.QC.outlier.low
    touch ${cohort}/${cohort}.scramble.QC.outlier.high
    touch ${cohort}/${cohort}.wham.QC.outlier.low
    touch ${cohort}/${cohort}.wham.QC.outlier.high
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
