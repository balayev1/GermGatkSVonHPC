#!/usr/bin/env nextflow

import groovy.json.JsonOutput

process GATKSV_GATHERSAMPLEEVIDENCE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path(sd_locs_vcf)

    output:
    tuple val(meta), path("${meta.id}/${meta.id}.counts.tsv.gz"), emit: coverage_counts
    tuple val(meta), path("${meta.id}/${meta.id}.pe.txt.gz"), emit: pesr_disc
    tuple val(meta), path("${meta.id}/${meta.id}.pe.txt.gz.tbi"), emit: pesr_disc_index
    tuple val(meta), path("${meta.id}/${meta.id}.sr.txt.gz"), emit: pesr_split
    tuple val(meta), path("${meta.id}/${meta.id}.sr.txt.gz.tbi"), emit: pesr_split_index
    tuple val(meta), path("${meta.id}/${meta.id}.sd.txt.gz"), emit: pesr_sd
    tuple val(meta), path("${meta.id}/${meta.id}.sd.txt.gz.tbi"), emit: pesr_sd_index
    tuple val(meta), path("${meta.id}/${meta.id}.manta.vcf.gz"), emit: manta_vcf
    tuple val(meta), path("${meta.id}/${meta.id}.manta.vcf.gz.tbi"), emit: manta_index
    tuple val(meta), path("${meta.id}/${meta.id}.wham.vcf.gz"), emit: wham_vcf
    tuple val(meta), path("${meta.id}/${meta.id}.wham.vcf.gz.tbi"), emit: wham_index
    tuple val(meta), path("${meta.id}/${meta.id}.scramble.vcf.gz"), emit: scramble_vcf
    tuple val(meta), path("${meta.id}/${meta.id}.scramble.vcf.gz.tbi"), emit: scramble_index
    tuple val(meta), path("${meta.id}/qc_metrics/${meta.id}.pe-file.tsv"), emit: pe_metrics, optional: true
    tuple val(meta), path("${meta.id}/qc_metrics/${meta.id}.sr-file.tsv"), emit: sr_metrics, optional: true
    tuple val(meta), path("${meta.id}/qc_metrics/${meta.id}.raw-counts.tsv"), emit: counts_metrics, optional: true
    tuple val(meta), path("${meta.id}/qc_metrics/manta_${meta.id}.vcf.tsv"), emit: manta_metrics, optional: true
    tuple val(meta), path("${meta.id}/qc_metrics/wham_${meta.id}.vcf.tsv"), emit: wham_metrics, optional: true
    tuple val(meta), path("${meta.id}/qc_metrics/scramble_${meta.id}.vcf.tsv"), emit: scramble_metrics, optional: true
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_json = "${prefix}_inputs.json"
    def static_json = JsonOutput.toJson(params.tool_inputs?.gse ?: [:])
    def bam_path = bam.toRealPath().toString()
    def bai_path = bai.toRealPath().toString()
    def sd_locs_vcf_path = sd_locs_vcf.toRealPath().toString()

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATKSV GATHERSAMPLEEVIDENCE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = Math.min(8192, (task.memory.mega * 0.8).intValue())
    }

    """
    mkdir -p ${prefix}

    render_json_template.py \
        --template ${params.gse_template} \
        --out ${out_json} \
        --static-json '${static_json}' \
        --set "GatherSampleEvidence.sample_id=${prefix}" \
        --set "GatherSampleEvidence.bam_or_cram_file=${bam_path}" \
        --set "GatherSampleEvidence.bam_or_cram_index=${bai_path}" \
        --set "GatherSampleEvidence.sd_locs_vcf=${sd_locs_vcf_path}"

    unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL

    # Run GSE using Cromwell
    java -Xmx${avail_mem}M \\
        -Dconfig.file=${params.cromwell_config} -jar ${params.cromwell_jar} \\
        run ${params.gse_wdl} \\
        -i ${out_json} \\
        -p ${params.deps_zip}

    # Copy output files required for downstream analysis to sample folder
    counts_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-CollectCounts/execution/${prefix}.counts.tsv.gz" || true)
    cp -L "\$counts_source" "${prefix}/${prefix}.counts.tsv.gz"

    pe_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-CollectSVEvidence/CollectSVEvidence/*/call-RunCollectSVEvidence/execution/${prefix}.pe.txt.gz" || true)
    pe_index_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-CollectSVEvidence/CollectSVEvidence/*/call-RunCollectSVEvidence/execution/${prefix}.pe.txt.gz.tbi" || true)
    cp -L "\$pe_source" "${prefix}/${prefix}.pe.txt.gz"
    cp -L "\$pe_index_source" "${prefix}/${prefix}.pe.txt.gz.tbi"

    sr_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-CollectSVEvidence/CollectSVEvidence/*/call-RunCollectSVEvidence/execution/${prefix}.sr.txt.gz" || true)
    sr_index_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-CollectSVEvidence/CollectSVEvidence/*/call-RunCollectSVEvidence/execution/${prefix}.sr.txt.gz.tbi" || true)
    cp -L "\$sr_source" "${prefix}/${prefix}.sr.txt.gz"
    cp -L "\$sr_index_source" "${prefix}/${prefix}.sr.txt.gz.tbi"

    sd_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-CollectSVEvidence/CollectSVEvidence/*/call-RunCollectSVEvidence/execution/${prefix}.sd.txt.gz" || true)
    sd_index_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-CollectSVEvidence/CollectSVEvidence/*/call-RunCollectSVEvidence/execution/${prefix}.sd.txt.gz.tbi" || true)
    cp -L "\$sd_source" "${prefix}/${prefix}.sd.txt.gz"
    cp -L "\$sd_index_source" "${prefix}/${prefix}.sd.txt.gz.tbi"

    manta_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-Manta/Manta/*/call-RunManta/execution/${prefix}.manta.vcf.gz" || true)
    manta_index_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-Manta/Manta/*/call-RunManta/execution/${prefix}.manta.vcf.gz.tbi" || true)
    cp -L "\$manta_source" "${prefix}/${prefix}.manta.vcf.gz"
    cp -L "\$manta_index_source" "${prefix}/${prefix}.manta.vcf.gz.tbi"

    wham_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-Whamg/Whamg/*/call-RunWhamgOnBam/execution/${prefix}.wham.vcf.gz" || true)
    wham_index_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-Whamg/Whamg/*/call-RunWhamgOnBam/execution/${prefix}.wham.vcf.gz.tbi" || true)
    cp -L "\$wham_source" "${prefix}/${prefix}.wham.vcf.gz"
    cp -L "\$wham_index_source" "${prefix}/${prefix}.wham.vcf.gz.tbi"

    scramble_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-Scramble/Scramble/*/call-MakeScrambleVcf/execution/${prefix}.scramble.vcf.gz" || true)
    scramble_index_source=\$(find cromwell-executions/GatherSampleEvidence -type f -path "*/call-Scramble/Scramble/*/call-MakeScrambleVcf/execution/${prefix}.scramble.vcf.gz.tbi" || true)
    cp -L "\$scramble_source" "${prefix}/${prefix}.scramble.vcf.gz"
    cp -L "\$scramble_index_source" "${prefix}/${prefix}.scramble.vcf.gz.tbi"

    metrics_found=0
    while IFS= read -r -d '' metric_file; do
        if [[ \$metrics_found -eq 0 ]]; then
            mkdir -p "${prefix}/qc_metrics"
            metrics_found=1
        fi
        cp -L "\$metric_file" "${prefix}/qc_metrics/"
    done < <(
        find cromwell-executions/GatherSampleEvidence \
            -type f \
            -path "*/call-GatherSampleEvidenceMetrics/GatherSampleEvidenceMetrics/*/call-*/execution/*" \
            -name "*.tsv" \
            -print0 2>/dev/null
    )


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/^.*version[[:space:]]*\"//; s/\".*\$//')
    END_VERSIONS
    """

    stub:
    def sample_id = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${sample_id}
    touch ${sample_id}/${sample_id}.counts.tsv.gz
    touch ${sample_id}/${sample_id}.pe.txt.gz
    touch ${sample_id}/${sample_id}.pe.txt.gz.tbi
    touch ${sample_id}/${sample_id}.sr.txt.gz
    touch ${sample_id}/${sample_id}.sr.txt.gz.tbi
    touch ${sample_id}/${sample_id}.sd.txt.gz
    touch ${sample_id}/${sample_id}.sd.txt.gz.tbi
    touch ${sample_id}/${sample_id}.manta.vcf.gz
    touch ${sample_id}/${sample_id}.manta.vcf.gz.tbi
    touch ${sample_id}/${sample_id}.wham.vcf.gz
    touch ${sample_id}/${sample_id}.wham.vcf.gz.tbi
    touch ${sample_id}/${sample_id}.scramble.vcf.gz
    touch ${sample_id}/${sample_id}.scramble.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: "stub"
    END_VERSIONS
    """
}
