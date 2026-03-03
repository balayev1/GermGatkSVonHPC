#!/usr/bin/env nextflow

process PREPARE_SD_LOCS_VCF {
    tag "sd_locs_vcf"
    label 'process_single'

    input:
    path(sd_locs_vcf)

    output:
    path("sd_locs.prepared.vcf"), emit: vcf
    path("versions.yml"), emit: versions

    script:
    def sd_locs_vcf_path = sd_locs_vcf.toRealPath().toString()

    """
    if [[ "${sd_locs_vcf_path}" == *.gz ]]; then
        gzip -cd "${sd_locs_vcf_path}" > sd_locs.prepared.vcf
    else
        cp "${sd_locs_vcf_path}" sd_locs.prepared.vcf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version 2>/dev/null | head -n 1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch sd_locs.prepared.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: "stub"
    END_VERSIONS
    """
}
