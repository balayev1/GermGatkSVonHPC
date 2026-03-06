//
// SAMPLE_PROCESSING
//

params.options = [:]
include { PICARD_COLLECTINSERTSIZEMETRICS   } from '../../modules/nf-core/picard/collectinsertsizemetrics/main.nf'
include { BCFTOOLS_ANNOTATE                 } from '../../modules/nf-core/bcftools/annotate/main.nf'
include { GATK_UPDATEVCFSEQUENCEDICTIONARY  } from '../../modules/local/update_vcf_dict.nf'
include { GATKSV_GATHERSAMPLEEVIDENCE       } from '../../modules/local/gather_sample_evidence.nf'
include { GATKSV_EVIDENCEQC as EVIDENCE_QC  } from '../../modules/local/evidence_qc.nf'
include { SAMPLE_QC                         } from '../../modules/local/sample_qc.nf'
include { BATCHING                          } from '../../modules/local/batching.nf'

workflow SAMPLE_PROCESSING {
    take:
    ch_sample       // channel: [ meta, file(row.bam_or_cram), file(row.bai_or_crai) ]
    ch_sd_locs_vcf
    sd_locs_vcf_index
    ch_reference_dict

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
    // Add 'chr' to chromosome names as in reference genome fasta 
    //
    def map_string = (1..22).collect { "${it} chr${it}" }.join('\n') + '\nX chrX\nY chrY\n'
    ch_rename_map = Channel.of(map_string).collectFile(name: 'chromosome_map.txt', newLine: true)

    ch_bcftools_input = ch_sd_locs_vcf
        .combine(sd_locs_vcf_index)
        .combine(ch_rename_map)
        .map { vcf, index, rename_map ->
            tuple(
                [ id: 'dbsnp_common' ], // meta
                vcf,                    // input
                index,                  // index
                [],                     // annotations
                [],                     // annotations_index
                [],                     // columns
                [],                     // header_lines
                rename_map              // rename_chrs
            )
        }

    BCFTOOLS_ANNOTATE (
        ch_bcftools_input
    )
    versions = versions.mix(BCFTOOLS_ANNOTATE.out.versions_bcftools)

    GATK_UPDATEVCFSEQUENCEDICTIONARY (
        BCFTOOLS_ANNOTATE.out.vcf,
        ch_reference_dict
    )
    versions = versions.mix(GATK_UPDATEVCFSEQUENCEDICTIONARY.out.versions)

    ch_prepared_vcf = GATK_UPDATEVCFSEQUENCEDICTIONARY.out.vcf
        .map { meta, vcf -> vcf }
        .first()

    ch_sample_with_vcf = ch_sample.combine(ch_prepared_vcf)

    //
    // Collect per-sample SV evidence using GATKSV_GATHERSAMPLEEVIDENCE
    //
    GATKSV_GATHERSAMPLEEVIDENCE (
        ch_sample_with_vcf.map { meta, bam, bai, vcf -> tuple(meta, bam, bai) },
        ch_sample_with_vcf.map { meta, bam, bai, vcf -> vcf }
    )
    versions = versions.mix(GATKSV_GATHERSAMPLEEVIDENCE.out.versions)

    def sampleKeyForMeta = { meta -> "${meta.cohort}::${meta.id}" }
    def keyedEvidence = { ch -> ch.map { meta, evidence_file -> tuple(sampleKeyForMeta(meta), evidence_file) } }
    def asList = { value -> value instanceof List ? value : [value] }
    def requireSingleSampleFile = { sampleKey, label, value ->
        def files = asList(value).findAll { it != null }
        if (files.size() != 1) {
            throw new IllegalStateException("${label} resolution mismatch for '${sampleKey}': expected 1 file, found ${files.size()}")
        }
        files[0]
    }
    def keyedSingleSampleEvidence = { ch, label ->
        ch.map { meta, evidence_file ->
            def sampleKey = sampleKeyForMeta(meta)
            tuple(sampleKey, requireSingleSampleFile(sampleKey, label, evidence_file))
        }
    }
    def keyedEvidenceQcFile = { ch, label ->
        ch.map { cohort, sample_ids, file_value ->
            tuple(cohort, sample_ids, requireSingleSampleFile(cohort.toString(), label, file_value))
        }
    }

    gse_outfiles = GATKSV_GATHERSAMPLEEVIDENCE.out.counts
        .map { meta, counts_file -> tuple(sampleKeyForMeta(meta), meta, counts_file) }
        .join(keyedEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.manta_vcf))
        .join(keyedEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.wham_vcf))
        .join(keyedEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.scramble_vcf))
        .map { _sample_key, meta, counts_file, manta_file, wham_file, scramble_file -> tuple(meta, counts_file, manta_file, wham_file, scramble_file) }

    evidqc_input = gse_outfiles
        .map { meta, counts_file, manta_file, wham_file, scramble_file -> tuple(meta.cohort, meta.id, counts_file, manta_file, wham_file, scramble_file) }
        .groupTuple()

    counts_by_sample = gse_outfiles
        .map { meta, counts_file, _manta_file, _wham_file, _scramble_file -> tuple(meta.cohort.toString(), meta.batch.toString(), meta.id.toString(), counts_file) }

    evidence_by_sample_keyed = gse_outfiles
        .map { meta, counts_file, manta_file, wham_file, scramble_file ->
            tuple(sampleKeyForMeta(meta), meta.cohort.toString(), meta.batch.toString(), meta.id.toString(), counts_file, manta_file, wham_file, scramble_file)
        }

    evidence_with_required_files_keyed = evidence_by_sample_keyed
        .join(keyedSingleSampleEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.pe_file, "PE file"))
        .join(keyedSingleSampleEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.sr_file, "SR file"))
        .join(keyedSingleSampleEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.sd_file, "SD file"))
        .map { sample_key, cohort, batch, sample_id, counts_file, manta_file, wham_file, scramble_file, pe_file, sr_file, sd_file ->
            tuple(sample_key, cohort, batch, sample_id, counts_file, pe_file, sr_file, sd_file, manta_file, wham_file, scramble_file)
        }

    //
    // Ploidy estimation, dosage scoring and sex assignment using GATKSV_EVIDENCEQC      
    //
    EVIDENCE_QC (
        evidqc_input
    )
    versions = versions.mix(EVIDENCE_QC.out.versions)
    evidence_qc_table = keyedEvidenceQcFile(EVIDENCE_QC.out.evidence_qc_table, "EvidenceQC table")
    sample_sex_assignments = keyedEvidenceQcFile(EVIDENCE_QC.out.sample_sex_assignments, "sample_sex_assignments")
    passing_samples_metadata = keyedEvidenceQcFile(EVIDENCE_QC.out.passing_samples_metadata, "passing_samples_metadata")

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

    sample_qc_input = evidence_qc_table
        .join(sample_sex_assignments)
        .map { cohort, sample_ids, evidence_qc_table_file, _sample_ids2, sample_sex_assignments_file ->
            tuple(cohort, sample_ids, evidence_qc_table_file, sample_sex_assignments_file)
        }
        .join(insert_metrics_by_cohort)
        .join(ped_by_cohort)
        .map { cohort, sample_ids, evidence_qc_table_file, sample_sex_assignments_file, _picard_sample_ids, insert_size_metrics, ped_file ->
            tuple(cohort, sample_ids, insert_size_metrics, ped_file, evidence_qc_table_file, sample_sex_assignments_file)
        }

    SAMPLE_QC (
        sample_qc_input
    )
    versions = versions.mix(SAMPLE_QC.out.versions)

    counts_by_sample_keyed = counts_by_sample
        .map { cohort, batch, sample_id, counts_file -> tuple("${cohort}::${sample_id}", cohort, batch, sample_id, counts_file) }

    train_gcnv_input = Channel.empty()
    model_batch_assignments_keyed = Channel.empty()
    if (params.run_batching) {
        batching_input = passing_samples_metadata
            .map { cohort, sample_ids, passing_samples_metadata_file -> tuple(cohort, passing_samples_metadata_file) }
            .join(ped_by_cohort)
            .map { cohort, passing_samples_metadata_file, ped_file -> tuple(cohort, passing_samples_metadata_file, ped_file) }

        BATCHING (
            batching_input
        )
        versions = versions.mix(BATCHING.out.versions)

        batch_assignments_keyed = BATCHING.out.batch_assignments
            .flatMap { cohort, batch_assignments_file ->
                def cohort_name = cohort.toString()
                def rows = batch_assignments_file.readLines()
                if (rows.size() <= 1) {
                    return []
                }
                rows.drop(1).findAll { it?.trim() }.collect { row ->
                    def fields = row.split('\t')
                    if (fields.size() < 2) {
                        throw new IllegalStateException("Malformed row in ${batch_assignments_file}: ${row}")
                    }
                    def batch = fields[0].trim()
                    def sample = fields[1].trim()
                    tuple("${cohort_name}::${sample}", cohort_name, batch, sample)
                }
            }

        model_batch_assignments_keyed = batch_assignments_keyed

        train_gcnv_input = batch_assignments_keyed
            .join(counts_by_sample_keyed)
            .map { sample_key, cohort, batch, sample_id, cohort2, batch2, sample_id2, counts_file ->
                if (cohort != cohort2 || sample_id != sample_id2) {
                    throw new IllegalStateException("Batch assignment join mismatch for '${sample_key}'")
                }
                tuple("${cohort}__${batch}", sample_id2, counts_file)
            }
            .groupTuple()
    } else {
        // No automatic batching requested: train one model per input cohort+batch.
        model_batch_assignments_keyed = counts_by_sample_keyed
            .map { sample_key, cohort, batch, sample_id, counts_file ->
                tuple(sample_key, cohort, batch, sample_id)
            }

        train_gcnv_input = counts_by_sample_keyed
            .map { sample_key, cohort, batch, sample_id, counts_file ->
                tuple("${cohort}__${batch}", sample_id, counts_file)
            }
            .groupTuple()
    }

    batch_evidence_input = model_batch_assignments_keyed
        .join(evidence_with_required_files_keyed)
        .map { sample_key, cohort, batch, sample_id, cohort2, batch2, sample_id2, counts_file, pe_file, sr_file, sd_file, manta_file, wham_file, scramble_file ->
            if (cohort != cohort2 || sample_id != sample_id2) {
                throw new IllegalStateException("Batch evidence join mismatch for '${sample_key}'")
            }
            tuple("${cohort}__${batch}", cohort, sample_id2, counts_file, pe_file, sr_file, sd_file, manta_file, wham_file, scramble_file)
        }
        .groupTuple()
        .map { batch_key, cohorts, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files ->
            def cohort_values = cohorts.collect { it.toString() }.unique()
            if (cohort_values.size() != 1) {
                throw new IllegalStateException("Batch '${batch_key}' resolved to multiple cohorts: ${cohort_values.join(', ')}")
            }
            tuple(cohort_values[0], batch_key, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files)
        }
        .join(ped_by_cohort)
        .map { cohort, batch_key, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files, ped_file ->
            tuple(batch_key, sample_ids, ped_file, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files)
        }

    emit:
    versions
    train_gcnv_input
    batch_evidence_input
}
