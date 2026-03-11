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
    def evidenceQcBatchKeyForMeta = { meta -> "${meta.cohort}__${meta.batch}" }
    def cohortFromBatchKey = { batchKey ->
        def key = batchKey.toString()
        def sep = '__'
        if (!key.contains(sep)) {
            throw new IllegalStateException("EvidenceQC batch key '${key}' does not contain expected '${sep}' separator")
        }
        key.split(sep, 2)[0]
    }
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
    gse_outfiles = GATKSV_GATHERSAMPLEEVIDENCE.out.coverage_counts
        .map { meta, counts_file -> tuple(sampleKeyForMeta(meta), meta, counts_file) }
        .join(keyedEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.manta_vcf))
        .join(keyedEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.wham_vcf))
        .join(keyedEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.scramble_vcf))
        .map { _sample_key, meta, counts_file, manta_file, wham_file, scramble_file -> tuple(meta, counts_file, manta_file, wham_file, scramble_file) }

    evidqc_input = gse_outfiles
        .map { meta, counts_file, manta_file, wham_file, scramble_file -> tuple(evidenceQcBatchKeyForMeta(meta), meta.id, counts_file, manta_file, wham_file, scramble_file) }
        .groupTuple()

    counts_by_sample = gse_outfiles
        .map { meta, counts_file, _manta_file, _wham_file, _scramble_file -> tuple(meta.cohort.toString(), meta.batch.toString(), meta.id.toString(), counts_file) }

    evidence_by_sample_keyed = gse_outfiles
        .map { meta, counts_file, manta_file, wham_file, scramble_file ->
            tuple(sampleKeyForMeta(meta), meta.cohort.toString(), meta.batch.toString(), meta.id.toString(), counts_file, manta_file, wham_file, scramble_file)
        }

    evidence_with_required_files_keyed = evidence_by_sample_keyed
        .join(keyedSingleSampleEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.pesr_disc, "PE file"))
        .join(keyedSingleSampleEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.pesr_split, "SR file"))
        .join(keyedSingleSampleEvidence(GATKSV_GATHERSAMPLEEVIDENCE.out.pesr_sd, "SD file"))
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

    def aggregateEvidenceQcByBatch = { ch, label ->
        ch.map { batch_key, sample_ids, file_value ->
                tuple(
                    batch_key.toString(),
                    asList(sample_ids).collect { it.toString() },
                    requireSingleSampleFile(batch_key.toString(), label, file_value)
                )
            }
            .groupTuple()
            .map { batch_key, sample_id_groups, files ->
                def merged_ids = sample_id_groups.collectMany { asList(it) }.collect { it.toString() }.unique()
                tuple(batch_key.toString(), merged_ids, asList(files))
            }
    }

    evidence_qc_table_by_batch = aggregateEvidenceQcByBatch(EVIDENCE_QC.out.qc_table, "EvidenceQC table")
    sample_sex_assignments_by_batch = aggregateEvidenceQcByBatch(EVIDENCE_QC.out.ploidy_plots, "sample_sex_assignments")

    insert_metrics_by_batch = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics
        .map { meta, metrics_file -> tuple(evidenceQcBatchKeyForMeta(meta), meta.id, metrics_file) }
        .groupTuple()

    ped_by_batch = ch_sample
        .map { meta, bam, bai -> tuple(evidenceQcBatchKeyForMeta(meta), meta.cohort.toString(), meta.ped) }
        .groupTuple()
        .map { batch_key, cohorts, ped_files ->
            def cohort_values = cohorts.collect { it.toString() }.unique()
            if (cohort_values.size() != 1) {
                throw new IllegalStateException("Batch '${batch_key}' resolved to multiple cohorts: ${cohort_values.join(', ')}")
            }
            def unique_peds = ped_files.collect { it.toString() }.unique()
            if (unique_peds.size() != 1) {
                exit 1, "ERROR: Batch '${batch_key}' has multiple PED files: ${unique_peds.join(', ')}"
            }
            tuple(batch_key.toString(), cohort_values[0], file(unique_peds[0]))
        }

    sample_qc_input = evidence_qc_table_by_batch
        .join(sample_sex_assignments_by_batch)
        .map { batch_key, sample_ids_1, evidence_qc_table_files, _sample_ids_2, sample_sex_assignments_files ->
            def merged_sample_ids = (asList(sample_ids_1) + asList(_sample_ids_2)).collect { it.toString() }.unique()
            tuple(batch_key.toString(), merged_sample_ids, evidence_qc_table_files, sample_sex_assignments_files)
        }
        .join(insert_metrics_by_batch)
        .join(ped_by_batch)
        .map { batch_key, sample_ids, evidence_qc_table_files, sample_sex_assignments_files, _insert_metric_sample_ids, insert_size_metrics, _cohort_name, ped_file ->
            tuple(batch_key, sample_ids, insert_size_metrics, ped_file, evidence_qc_table_files, sample_sex_assignments_files)
        }

    SAMPLE_QC (
        sample_qc_input
    )
    versions = versions.mix(SAMPLE_QC.out.versions)

    ch_updated_ped_by_batch = SAMPLE_QC.out.updated_ped
        .map { batch_key, ped_file -> tuple(batch_key.toString(), ped_file) }

    sample_qc_outlier_ids_by_batch = SAMPLE_QC.out.sample_qc_reports
        .filter { batch_key, report_file -> report_file.getName() == "Excluded_Sample_ID_only.tsv" }
        .map { batch_key, report_file -> tuple(batch_key.toString(), report_file) }

    passing_samples_metadata_by_batch = SAMPLE_QC.out.passing_samples_metadata
        .map { batch_key, sample_ids, passing_samples_metadata_file ->
            tuple(
                batch_key.toString(),
                asList(sample_ids).collect { it.toString() },
                [passing_samples_metadata_file]
            )
        }

    batch_key_to_cohort = ped_by_batch
        .map { batch_key, cohort, ped_file -> tuple(batch_key.toString(), cohort.toString()) }
        .unique()

    // Keep a cohort-keyed PED channel for downstream cohort-level optional steps.
    ch_updated_ped_by_cohort = ch_updated_ped_by_batch
        .join(batch_key_to_cohort)
        .map { batch_key, ped_file, cohort -> tuple(cohort.toString(), ped_file) }
        .groupTuple()
        .map { cohort, ped_files ->
            def unique_peds = asList(ped_files).collect { file(it.toString()) }.unique { it.toString() }
            tuple(cohort.toString(), unique_peds[0])
        }

    counts_by_sample_keyed = counts_by_sample
        .map { cohort, batch, sample_id, counts_file ->
            def source_batch_key = "${cohort}__${batch}"
            tuple("${cohort}::${sample_id}", source_batch_key, cohort, batch, sample_id, counts_file)
        }

    train_gcnv_input_base = Channel.empty()
    model_batch_assignments_keyed = Channel.empty()
    if (params.run_batching) {
        batching_input = passing_samples_metadata_by_batch
            .map { batch_key, sample_ids, passing_samples_metadata_files -> tuple(batch_key, passing_samples_metadata_files) }
            .join(ped_by_batch)
            .map { batch_key, passing_samples_metadata_files, _cohort_name, ped_file -> tuple(batch_key, passing_samples_metadata_files, ped_file) }

        BATCHING (
            batching_input
        )
        versions = versions.mix(BATCHING.out.versions)

        batch_assignments_keyed = BATCHING.out.batch_assignments
            .flatMap { source_batch_key, batch_assignments_file ->
                def cohort_name = cohortFromBatchKey(source_batch_key)
                def rows = batch_assignments_file.readLines()
                if (rows.size() <= 1) {
                    return []
                }
                rows.drop(1).findAll { it?.trim() }.collect { row ->
                    def fields = row.split('\t')
                    if (fields.size() < 2) {
                        throw new IllegalStateException("Malformed row in ${batch_assignments_file}: ${row}")
                    }
                    def model_batch_label = fields[0].trim()
                    def sample = fields[1].trim()
                    tuple("${cohort_name}::${sample}", source_batch_key.toString(), model_batch_label, sample, cohort_name)
                }
            }

        model_batch_assignments_keyed = batch_assignments_keyed
            .join(counts_by_sample_keyed)
            .map { sample_key, source_batch_key, model_batch_label, sample_id, cohort_from_batching, source_batch_key2, cohort, batch, sample_id2, counts_file ->
                if (source_batch_key != source_batch_key2 || cohort_from_batching != cohort || sample_id != sample_id2) {
                    throw new IllegalStateException("Batch assignment join mismatch for '${sample_key}'")
                }
                def model_batch_key = "${source_batch_key}__${model_batch_label}"
                tuple(sample_key, model_batch_key, source_batch_key, cohort, sample_id2, counts_file)
            }

        train_gcnv_input_base = model_batch_assignments_keyed
            .map { sample_key, model_batch_key, source_batch_key, cohort, sample_id, counts_file ->
                tuple(model_batch_key, sample_id, counts_file)
            }
            .groupTuple()
    } else {
        // No automatic batching requested: one model per source cohort+batch key.
        model_batch_assignments_keyed = counts_by_sample_keyed
            .map { sample_key, source_batch_key, cohort, batch, sample_id, counts_file ->
                tuple(sample_key, source_batch_key, source_batch_key, cohort, sample_id, counts_file)
            }

        train_gcnv_input_base = model_batch_assignments_keyed
            .map { sample_key, model_batch_key, source_batch_key, cohort, sample_id, counts_file ->
                tuple(model_batch_key, sample_id, counts_file)
            }
            .groupTuple()
    }

    model_batch_to_source_batch = model_batch_assignments_keyed
        .map { sample_key, model_batch_key, source_batch_key, cohort, sample_id, counts_file ->
            tuple(model_batch_key, source_batch_key)
        }
        .unique()

    outlier_sample_ids_by_batch = model_batch_to_source_batch
        .map { model_batch_key, source_batch_key -> tuple(source_batch_key, model_batch_key) }
        .join(sample_qc_outlier_ids_by_batch)
        .map { source_batch_key, model_batch_key, outlier_sample_ids_file -> tuple(model_batch_key, outlier_sample_ids_file) }

    train_gcnv_input = train_gcnv_input_base
        .join(outlier_sample_ids_by_batch)
        .map { batch_key, sample_ids, count_files, outlier_sample_ids_file ->
            tuple(batch_key, sample_ids, count_files, outlier_sample_ids_file)
        }

    batch_evidence_input = model_batch_assignments_keyed
        .join(evidence_with_required_files_keyed)
        .map { sample_key, model_batch_key, source_batch_key, cohort, sample_id, counts_file, cohort2, batch2, sample_id2, counts_file2, pe_file, sr_file, sd_file, manta_file, wham_file, scramble_file ->
            if (cohort != cohort2 || sample_id != sample_id2) {
                throw new IllegalStateException("Batch evidence join mismatch for '${sample_key}'")
            }
            def expected_source_batch_key = "${cohort2}__${batch2}"
            if (source_batch_key != expected_source_batch_key) {
                throw new IllegalStateException("Source batch key mismatch for '${sample_key}': expected '${expected_source_batch_key}', got '${source_batch_key}'")
            }
            tuple(model_batch_key, source_batch_key, cohort, sample_id2, counts_file2, pe_file, sr_file, sd_file, manta_file, wham_file, scramble_file)
        }
        .groupTuple()
        .map { model_batch_key, source_batch_keys, cohorts, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files ->
            def source_values = source_batch_keys.collect { it.toString() }.unique()
            if (source_values.size() != 1) {
                throw new IllegalStateException("Model batch '${model_batch_key}' resolved to multiple source batches: ${source_values.join(', ')}")
            }
            def cohort_values = cohorts.collect { it.toString() }.unique()
            if (cohort_values.size() != 1) {
                throw new IllegalStateException("Model batch '${model_batch_key}' resolved to multiple cohorts: ${cohort_values.join(', ')}")
            }
            tuple(model_batch_key, source_values[0], cohort_values[0], sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files)
        }
        .map { model_batch_key, source_batch_key, cohort, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files ->
            tuple(source_batch_key, model_batch_key, cohort, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files)
        }
        .join(ch_updated_ped_by_batch)
        .map { source_batch_key, model_batch_key, cohort, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files, ped_file ->
            tuple(model_batch_key, sample_ids, ped_file, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files)
        }

    emit:
    versions
    train_gcnv_input
    batch_evidence_input
    outlier_sample_ids_by_batch
    updated_ped = ch_updated_ped_by_cohort
    evidence_qc_manta_variant_counts = EVIDENCE_QC.out.manta_variant_counts
    evidence_qc_wham_variant_counts = EVIDENCE_QC.out.wham_variant_counts
    evidence_qc_scramble_variant_counts = EVIDENCE_QC.out.scramble_variant_counts
}
