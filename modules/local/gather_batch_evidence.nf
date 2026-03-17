#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GATHERBATCHEVIDENCE {
    tag "${batch_name}"
    label 'process_medium'

    input:
    tuple val(batch_name), val(sample_ids), path(ped_file), path(count_files), path(pe_files), path(sr_files), path(sd_files), path(manta_vcf_files), path(wham_vcf_files), path(scramble_vcf_files), path(gcnv_model_tars), path(contig_ploidy_model_tar)

    output:
    tuple val(batch_name), path("${batch_name}/${batch_name}.baf.txt.gz"), emit: merged_BAF
    tuple val(batch_name), path("${batch_name}/${batch_name}.baf.txt.gz.tbi"), emit: merged_BAF_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.sr.txt.gz"), emit: merged_SR
    tuple val(batch_name), path("${batch_name}/${batch_name}.sr.txt.gz.tbi"), emit: merged_SR_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.pe.txt.gz"), emit: merged_PE
    tuple val(batch_name), path("${batch_name}/${batch_name}.pe.txt.gz.tbi"), emit: merged_PE_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.RD.txt.gz"), emit: merged_bincov
    tuple val(batch_name), path("${batch_name}/${batch_name}.RD.txt.gz.tbi"), emit: merged_bincov_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.DEL.bed.gz"), emit: merged_dels
    tuple val(batch_name), path("${batch_name}/${batch_name}.DUP.bed.gz"), emit: merged_dups
    tuple val(batch_name), path("${batch_name}/${batch_name}.DEL.header.bed.gz"), emit: cnmops_del
    tuple val(batch_name), path("${batch_name}/${batch_name}.DEL.header.bed.gz.tbi"), emit: cnmops_del_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.DUP.header.bed.gz"), emit: cnmops_dup
    tuple val(batch_name), path("${batch_name}/${batch_name}.DUP.header.bed.gz.tbi"), emit: cnmops_dup_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.DEL.large.bed.gz"), emit: cnmops_large_del
    tuple val(batch_name), path("${batch_name}/${batch_name}.DEL.large.bed.gz.tbi"), emit: cnmops_large_del_index
    tuple val(batch_name), path("${batch_name}/${batch_name}.DUP.large.bed.gz"), emit: cnmops_large_dup
    tuple val(batch_name), path("${batch_name}/${batch_name}.DUP.large.bed.gz.tbi"), emit: cnmops_large_dup_index
    tuple val(batch_name), path("${batch_name}/${batch_name}_medianCov.transposed.bed"), emit: median_cov
    tuple val(batch_name), path("${batch_name}/${batch_name}.manta.std.tar.gz"), emit: std_manta_vcf_tar
    tuple val(batch_name), path("${batch_name}/${batch_name}.wham.std.tar.gz"), emit: std_wham_vcf_tar
    tuple val(batch_name), path("${batch_name}/${batch_name}.scramble.std.tar.gz"), emit: std_scramble_vcf_tar
    tuple val(batch_name), path("${batch_name}/${batch_name}_ploidy_matrix.bed.gz"), emit: batch_ploidy_matrix, optional: true
    tuple val(batch_name), path("${batch_name}/${batch_name}_ploidy_plots.tar.gz"), emit: batch_ploidy_plots, optional: true
    tuple val(batch_name), path("${batch_name}/${batch_name}.BAF.QC_matrix.txt"), emit: BAF_stats, optional: true
    tuple val(batch_name), path("${batch_name}/${batch_name}.RD.QC_matrix.txt"), emit: RD_stats, optional: true
    tuple val(batch_name), path("${batch_name}/${batch_name}.PE.QC_matrix.txt"), emit: PE_stats, optional: true
    tuple val(batch_name), path("${batch_name}/${batch_name}.SR.QC_matrix.txt"), emit: SR_stats, optional: true
    tuple val(batch_name), path("${batch_name}/${batch_name}.00_matrix_FC_QC.png"), emit: Matrix_QC_plot, optional: true
    tuple val(batch_name), path("${batch_name}/GatherBatchEvidence.${batch_name}.metrics.tsv"), emit: metrics_file_batchevidence, optional: true
    tuple val(batch_name), path("${batch_name}/tloc_*.manta.complex.vcf.gz"), emit: manta_tloc, optional: true
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

    mkdir -p "${batch_id}"

    copy_outputs() {
        local pattern="\$1"
        local required="\${2:-1}"
        local found=0
        while IFS= read -r -d '' source; do
            found=1
            cp -L "\$source" "${batch_id}/\$(basename "\$source")"
        done < <(find cromwell-executions/GatherBatchEvidence -type f -name "\${pattern}" -print0)

        if [[ "\$required" -eq 1 && "\$found" -eq 0 ]]; then
            echo "ERROR: Expected GatherBatchEvidence output(s) not found for pattern: \${pattern}" >&2
            exit 1
        fi
    }

    copy_outputs "${batch_id}.baf.txt.gz" 1
    copy_outputs "${batch_id}.baf.txt.gz.tbi" 1
    copy_outputs "${batch_id}.sr.txt.gz" 1
    copy_outputs "${batch_id}.sr.txt.gz.tbi" 1
    copy_outputs "${batch_id}.pe.txt.gz" 1
    copy_outputs "${batch_id}.pe.txt.gz.tbi" 1
    copy_outputs "${batch_id}.RD.txt.gz" 1
    copy_outputs "${batch_id}.RD.txt.gz.tbi" 1
    copy_outputs "${batch_id}.DEL.bed.gz" 1
    copy_outputs "${batch_id}.DUP.bed.gz" 1
    copy_outputs "${batch_id}.DEL.header.bed.gz" 1
    copy_outputs "${batch_id}.DEL.header.bed.gz.tbi" 1
    copy_outputs "${batch_id}.DUP.header.bed.gz" 1
    copy_outputs "${batch_id}.DUP.header.bed.gz.tbi" 1
    copy_outputs "${batch_id}.DEL.large.bed.gz" 1
    copy_outputs "${batch_id}.DEL.large.bed.gz.tbi" 1
    copy_outputs "${batch_id}.DUP.large.bed.gz" 1
    copy_outputs "${batch_id}.DUP.large.bed.gz.tbi" 1
    copy_outputs "${batch_id}_medianCov.transposed.bed" 1
    copy_outputs "${batch_id}.manta.std.tar.gz" 1
    copy_outputs "${batch_id}.wham.std.tar.gz" 1
    copy_outputs "${batch_id}.scramble.std.tar.gz" 1
    copy_outputs "${batch_id}_ploidy_matrix.bed.gz" 0
    copy_outputs "${batch_id}_ploidy_plots.tar.gz" 0
    copy_outputs "${batch_id}.BAF.QC_matrix.txt" 0
    copy_outputs "${batch_id}.RD.QC_matrix.txt" 0
    copy_outputs "${batch_id}.PE.QC_matrix.txt" 0
    copy_outputs "${batch_id}.SR.QC_matrix.txt" 0
    copy_outputs "${batch_id}.00_matrix_FC_QC.png" 0
    copy_outputs "GatherBatchEvidence.${batch_id}.metrics.tsv" 0
    copy_outputs "tloc_*.manta.complex.vcf.gz" 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${batch_name}
    touch ${batch_name}/${batch_name}.baf.txt.gz
    touch ${batch_name}/${batch_name}.baf.txt.gz.tbi
    touch ${batch_name}/${batch_name}.sr.txt.gz
    touch ${batch_name}/${batch_name}.sr.txt.gz.tbi
    touch ${batch_name}/${batch_name}.pe.txt.gz
    touch ${batch_name}/${batch_name}.pe.txt.gz.tbi
    touch ${batch_name}/${batch_name}.RD.txt.gz
    touch ${batch_name}/${batch_name}.RD.txt.gz.tbi
    touch ${batch_name}/${batch_name}.DEL.bed.gz
    touch ${batch_name}/${batch_name}.DUP.bed.gz
    touch ${batch_name}/${batch_name}.DEL.header.bed.gz
    touch ${batch_name}/${batch_name}.DEL.header.bed.gz.tbi
    touch ${batch_name}/${batch_name}.DUP.header.bed.gz
    touch ${batch_name}/${batch_name}.DUP.header.bed.gz.tbi
    touch ${batch_name}/${batch_name}.DEL.large.bed.gz
    touch ${batch_name}/${batch_name}.DEL.large.bed.gz.tbi
    touch ${batch_name}/${batch_name}.DUP.large.bed.gz
    touch ${batch_name}/${batch_name}.DUP.large.bed.gz.tbi
    touch ${batch_name}/${batch_name}_medianCov.transposed.bed
    touch ${batch_name}/${batch_name}.manta.std.tar.gz
    touch ${batch_name}/${batch_name}.wham.std.tar.gz
    touch ${batch_name}/${batch_name}.scramble.std.tar.gz
    touch ${batch_name}/${batch_name}_ploidy_matrix.bed.gz
    touch ${batch_name}/${batch_name}_ploidy_plots.tar.gz
    touch ${batch_name}/${batch_name}.BAF.QC_matrix.txt
    touch ${batch_name}/${batch_name}.RD.QC_matrix.txt
    touch ${batch_name}/${batch_name}.PE.QC_matrix.txt
    touch ${batch_name}/${batch_name}.SR.QC_matrix.txt
    touch ${batch_name}/${batch_name}.00_matrix_FC_QC.png
    touch ${batch_name}/GatherBatchEvidence.${batch_name}.metrics.tsv
    touch ${batch_name}/tloc_001.manta.complex.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
