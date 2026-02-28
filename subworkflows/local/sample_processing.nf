//
// SAMPLE_PROCESSING
//

params.options = [:]
include { PICARD_COLLECTINSERTSIZEMETRICS   } from '../../modules/nf-core/picard/collectinsertsizemetrics/main.nf'
include { GATKSV_GATHERSAMPLEEVIDENCE       } from '../../modules/local/gather_sample_evidence.nf'
include { GATKSV_EVIDENCEQC as EVIDENCE_QC  } from '../../modules/local/evidence_qc.nf'
include { SAMPLE_QC                         } from '../../modules/local/sample_qc.nf'
include { BATCHING                          } from '../../modules/local/batching.nf'

workflow SAMPLE_PROCESSING {
    take:
    ch_sample       // channel: [ meta, file(row.bam_or_cram), file(row.bai_or_crai) ]

    main:
    versions     = Channel.empty()

    //
    // Collect bam/cram insert size info using PICARD_COLLECTINSERTSIZEMETRICS
    //
    PICARD_COLLECTINSERTSIZEMETRICS (
        ch_sample
    )
    versions = versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions_picard)

    //
    // Collect per-sample SV evidence using GATKSV_GATHERSAMPLEEVIDENCE
    //
    GATKSV_GATHERSAMPLEEVIDENCE (
        ch_sample
    )
    versions = versions.mix(GATKSV_GATHERSAMPLEEVIDENCE.out.versions)

    evidqc_input = GATKSV_GATHERSAMPLEEVIDENCE.out.gse_outfiles
        .map { meta, counts_file, manta_file, wham_file, scramble_file ->
            tuple(meta.cohort, meta.id, counts_file, manta_file, wham_file, scramble_file)
        }
        .groupTuple()

    //
    // Ploidy estimation, dosage scoring and sex assignment using GATKSV_EVIDENCEQC      
    //
    EVIDENCE_QC (
        evidqc_input
    )
    versions = versions.mix(EVIDENCE_QC.out.versions)

    insert_metrics_by_cohort = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics
        .map { meta, metrics_file -> tuple(meta.cohort, meta.id, metrics_file) }
        .groupTuple()

    ped_by_cohort = ch_sample
        .map { meta, bam, bai -> tuple(meta.cohort, meta.ped) }
        .groupTuple()
        .map { cohort, ped_files ->
            def unique_peds = ped_files.collect { it.toString() }.unique()
            if (unique_peds.size() != 1) {
                exit 1, "ERROR: Cohort '${cohort}' has multiple PED files: ${unique_peds.join(', ')}"
            }
            tuple(cohort, file(unique_peds[0]))
        }

    sample_qc_input = EVIDENCE_QC.out.evidence_qc_results
        .join(insert_metrics_by_cohort)
        .map { cohort, sample_ids, evidence_qc_results, picard_sample_ids, insert_size_metrics ->
            tuple(cohort, sample_ids, insert_size_metrics, evidence_qc_results)
        }
        .join(ped_by_cohort)
        .map { cohort, sample_ids, insert_size_metrics, evidence_qc_results, ped_file ->
            tuple(cohort, sample_ids, insert_size_metrics, ped_file, evidence_qc_results)
        }

    SAMPLE_QC (
        sample_qc_input
    )
    versions = versions.mix(SAMPLE_QC.out.versions)

    if (params.run_batching) {
        batching_input = EVIDENCE_QC.out.evidence_qc_results
            .map { cohort, sample_ids, evidence_qc_results -> tuple(cohort, evidence_qc_results) }
            .join(ped_by_cohort)
            .map { cohort, evidence_qc_results, ped_file -> tuple(cohort, evidence_qc_results, ped_file) }

        BATCHING (
            batching_input
        )
        versions = versions.mix(BATCHING.out.versions)
    }

    emit:
    versions
}
