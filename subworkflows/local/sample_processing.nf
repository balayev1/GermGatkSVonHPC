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
include { GATK_TRAINGCNV                    } from '../../modules/local/train_gcnv.nf'

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
        ch_reference_dict.first()
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

    gse_outfiles = GATKSV_GATHERSAMPLEEVIDENCE.out.gse_outfiles
    gse_outfiles.into { gse_for_evidqc; gse_for_traingcnv }

    evidqc_input = gse_for_evidqc
        .map { meta, counts_file, manta_file, wham_file, scramble_file ->
            tuple(meta.cohort, meta.id, counts_file, manta_file, wham_file, scramble_file)
        }
        .groupTuple()

    counts_by_sample = gse_for_traingcnv
        .map { meta, counts_file, manta_file, wham_file, scramble_file ->
            tuple(meta.cohort.toString(), meta.id.toString(), counts_file)
        }

    //
    // Ploidy estimation, dosage scoring and sex assignment using GATKSV_EVIDENCEQC      
    //
    EVIDENCE_QC (
        evidqc_input
    )
    versions = versions.mix(EVIDENCE_QC.out.versions)
    evidence_qc_results = EVIDENCE_QC.out.evidence_qc_results
    evidence_qc_results.into { evidence_qc_for_sample_qc; evidence_qc_for_batching }

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

    sample_qc_input = evidence_qc_for_sample_qc
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

    counts_by_sample_keyed = counts_by_sample
        .map { cohort, sample_id, counts_file -> tuple("${cohort}::${sample_id}", cohort, sample_id, counts_file) }

    def train_gcnv_input
    if (params.run_batching) {
        batching_input = evidence_qc_for_batching
            .map { cohort, sample_ids, evidence_qc_results -> tuple(cohort, evidence_qc_results) }
            .join(ped_by_cohort)
            .map { cohort, evidence_qc_results, ped_file -> tuple(cohort, evidence_qc_results, ped_file) }

        BATCHING (
            batching_input
        )
        versions = versions.mix(BATCHING.out.versions)

        batch_assignments_keyed = BATCHING.out.batching_out
            .flatMap { cohort, batching_out ->
                def cohort_name = cohort.toString()
                def batch_assignments_file = file("${batching_out}/batching/batch_assignments.tsv")
                if (!batch_assignments_file.exists()) {
                    throw new IllegalStateException("Batch assignments file not found: ${batch_assignments_file}")
                }
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

        train_gcnv_input = batch_assignments_keyed
            .join(counts_by_sample_keyed)
            .map { sample_key, cohort, batch, sample_id, cohort2, sample_id2, counts_file ->
                tuple("${cohort}__${batch}", sample_id, counts_file)
            }
            .groupTuple()
    } else {
        // No batching requested: train one model across all samples.
        train_gcnv_input = counts_by_sample_keyed
            .map { sample_key, cohort, sample_id, counts_file ->
                tuple("single_batch", sample_id, counts_file)
            }
            .groupTuple()
    }

    GATK_TRAINGCNV (
        train_gcnv_input
    )
    versions = versions.mix(GATK_TRAINGCNV.out.versions)

    emit:
    versions
    train_gcnv_results = GATK_TRAINGCNV.out.train_gcnv_results
}
