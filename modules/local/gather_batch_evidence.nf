#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GATHERBATCHEVIDENCE {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(sample_ids), path(ped_file), path(count_files), path(pe_files), path(sr_files), path(sd_files), path(manta_vcf_files), path(wham_vcf_files), path(scramble_vcf_files), path(gcnv_model_tars), path(contig_ploidy_model_tar)

    output:
    tuple val(batch_name), path("**/${batch_name}.baf.txt.gz"), emit: merged_BAF
    tuple val(batch_name), path("**/${batch_name}.sr.txt.gz"), emit: merged_SR
    tuple val(batch_name), path("**/${batch_name}.pe.txt.gz"), emit: merged_PE
    tuple val(batch_name), path("**/${batch_name}.RD.txt.gz"), emit: merged_bincov
    tuple val(batch_name), path("**/*.del.bed"), emit: merged_dels
    tuple val(batch_name), path("**/*.dup.bed"), emit: merged_dups
    tuple val(batch_name), path("**/*_medianCov.transposed.bed"), emit: median_cov
    tuple val(batch_name), path("**/*.manta_std.tar.gz"), emit: std_manta_vcf_tar
    tuple val(batch_name), path("**/*.wham_std.tar.gz"), emit: std_wham_vcf_tar
    tuple val(batch_name), path("**/*.scramble_std.tar.gz"), emit: std_scramble_vcf_tar
    tuple val(batch_name), path("**/${batch_name}_ploidy_matrix.bed.gz"), emit: batch_ploidy_matrix, optional: true
    tuple val(batch_name), path("**/${batch_name}_ploidy_plots.tar.gz"), emit: batch_ploidy_plots, optional: true
    tuple val(batch_name), path("**/${batch_name}.BAF.QC_matrix.txt"), emit: baf_qc_stats, optional: true
    tuple val(batch_name), path("**/${batch_name}.RD.QC_matrix.txt"), emit: rd_qc_stats, optional: true
    tuple val(batch_name), path("**/${batch_name}.PE.QC_matrix.txt"), emit: pe_qc_stats, optional: true
    tuple val(batch_name), path("**/${batch_name}.SR.QC_matrix.txt"), emit: sr_qc_stats, optional: true
    tuple val(batch_name), path("**/${batch_name}.00_matrix_FC_QC.png"), emit: plot_matrix_qc, optional: true
    tuple val(batch_name), path("**/tloc_*.vcf.gz"), emit: manta_tloc, optional: true
    path "versions.yml", emit: versions

    script:
    def template_path = file(params.gbe_template).toAbsolutePath()
    def static_json = JsonOutput.toJson(params.tool_inputs?.gather_batch_evidence ?: [:])

    def batch_id = batch_name?.toString()?.trim()
    if (!batch_id) {
        throw new IllegalArgumentException("GatherBatchEvidence batch name is required")
    }
    if (!(batch_id ==~ /^[A-Za-z0-9_]+$/)) {
        throw new IllegalArgumentException("GatherBatchEvidence batch '${batch_id}' is invalid. Batch IDs must contain only letters, numbers, and underscores.")
    }

    def asList = { value -> value instanceof List ? value : [value] }
    def toRealPathList = { value -> asList(value).collect { it.toRealPath().toString() } }

    def samples = asList(sample_ids).collect { it.toString() }
    def counts = toRealPathList(count_files)
    def pe = toRealPathList(pe_files)
    def sr = toRealPathList(sr_files)
    def sd = toRealPathList(sd_files)
    def manta = toRealPathList(manta_vcf_files)
    def wham = toRealPathList(wham_vcf_files)
    def scramble = toRealPathList(scramble_vcf_files)
    def gcnvModels = toRealPathList(gcnv_model_tars)
    def contigPloidyModel = contig_ploidy_model_tar.toRealPath().toString()
    def pedPath = ped_file.toRealPath().toString()

    if (!samples) {
        throw new IllegalArgumentException("GatherBatchEvidence received an empty sample list")
    }
    def n = samples.size()
    if (counts.size() != n || pe.size() != n || sr.size() != n || sd.size() != n || manta.size() != n || wham.size() != n || scramble.size() != n) {
        throw new IllegalArgumentException("GatherBatchEvidence input size mismatch: samples=${n}, counts=${counts.size()}, pe=${pe.size()}, sr=${sr.size()}, sd=${sd.size()}, manta=${manta.size()}, wham=${wham.size()}, scramble=${scramble.size()}")
    }
    if (!gcnvModels) {
        throw new IllegalArgumentException("GatherBatchEvidence received an empty gcnv_model_tars list")
    }

    def dynamic = [
        "GatherBatchEvidence.batch"                  : batch_id,
        "GatherBatchEvidence.ped_file"               : pedPath,
        "GatherBatchEvidence.samples"                : samples,
        "GatherBatchEvidence.counts"                 : counts,
        "GatherBatchEvidence.PE_files"               : pe,
        "GatherBatchEvidence.SR_files"               : sr,
        "GatherBatchEvidence.SD_files"               : sd,
        "GatherBatchEvidence.manta_vcfs"             : manta,
        "GatherBatchEvidence.wham_vcfs"              : wham,
        "GatherBatchEvidence.scramble_vcfs"          : scramble,
        "GatherBatchEvidence.gcnv_model_tars"        : gcnvModels,
        "GatherBatchEvidence.contig_ploidy_model_tar": contigPloidyModel
    ]
    file("gather_batch_evidence_dynamic.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(dynamic))

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV GATHERBATCHEVIDENCE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    render_json_template.py \\
        --template ${template_path} \\
        --out gather_batch_evidence_inputs.json \\
        --static-json '${static_json}' \\
        --merge-json-file gather_batch_evidence_dynamic.json

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    # Execute Cromwell
    java -Xmx${avail_mem}M -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.gbe_wdl} \\
        -i gather_batch_evidence_inputs.json \\
        -p ${params.deps_zip}

    mkdir -p gather_batch_evidence_results
    cp gather_batch_evidence_inputs.json gather_batch_evidence_results/
    find cromwell-executions/GatherBatchEvidence/ -name "call-*" -type d -exec cp -r {} gather_batch_evidence_results/ \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p gather_batch_evidence_results/call-stub/execution
    touch gather_batch_evidence_results/call-stub/.stub
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}.baf.txt.gz
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}.sr.txt.gz
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}.pe.txt.gz
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}.RD.txt.gz
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}.del.bed
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}.dup.bed
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}_medianCov.transposed.bed
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}.manta_std.tar.gz
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}.wham_std.tar.gz
    touch gather_batch_evidence_results/call-stub/execution/${batch_name}.scramble_std.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
