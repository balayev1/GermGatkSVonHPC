//
// JOINT_COHORT_CALLING
//

include { GATKSV_MERGEBATCHSITES } from '../../modules/local/merge_batch_sites.nf'
include { GATKSV_GENOTYPEBATCH } from '../../modules/local/genotype_batch.nf'

workflow JOINT_COHORT_CALLING {
    take:
    batch_cohort
    filtered_depth_vcf
    filtered_pesr_vcf
    ploidy_table
    cutoffs
    merged_PE
    merged_bincov
    merged_SR
    median_cov

    main:
    versions = Channel.empty()

    def merge_batch_sites_results = Channel.empty()
    def cohort_pesr_vcf = Channel.empty()
    def cohort_depth_vcf = Channel.empty()
    def merge_batch_sites_vcf = Channel.empty()
    def genotype_batch_results = Channel.empty()
    def genotyped_vcf = Channel.empty()

    def asList = { value ->
        if (value == null) {
            return []
        }
        value instanceof List ? value : [value]
    }
    def selectSinglePath = { value, label ->
        def candidates = asList(value).collect { file(it.toString()) }
        if (!candidates) {
            throw new IllegalStateException("Expected at least one file for ${label}")
        }
        if (candidates.size() == 1) {
            return candidates[0]
        }
        def preferred = candidates.findAll { it.getName().toLowerCase().contains("merged") || it.getName().toLowerCase().contains("cohort") }
        return preferred ? preferred[0] : candidates[0]
    }

    if (params.run_merge_batch_sites) {
        merge_batch_sites_input = batch_cohort
            .join(filtered_depth_vcf)
            .join(filtered_pesr_vcf)
            .join(ploidy_table)
            .map { batch_key, cohort, depth_vcf, pesr_vcf, ploidy_table_file ->
                tuple(cohort, depth_vcf, pesr_vcf, ploidy_table_file)
            }
            .groupTuple()
            .map { cohort, depth_vcfs, pesr_vcfs, ploidy_tables ->
                tuple(
                    cohort,
                    asList(ploidy_tables).collect { file(it.toString()) },
                    asList(depth_vcfs).collect { file(it.toString()) },
                    asList(pesr_vcfs).collect { file(it.toString()) }
                )
            }

        GATKSV_MERGEBATCHSITES (
            merge_batch_sites_input
        )
        versions = versions.mix(GATKSV_MERGEBATCHSITES.out.versions)
        merge_batch_sites_results = GATKSV_MERGEBATCHSITES.out.merge_batch_sites_results
        cohort_pesr_vcf = GATKSV_MERGEBATCHSITES.out.cohort_pesr_vcf
        cohort_depth_vcf = GATKSV_MERGEBATCHSITES.out.cohort_depth_vcf
        merge_batch_sites_vcf = GATKSV_MERGEBATCHSITES.out.merge_batch_sites_vcf

        if (params.run_genotype_batch) {
            genotype_batch_per_batch = batch_cohort
                .join(cutoffs)
                .join(merged_PE)
                .join(merged_bincov)
                .join(merged_SR)
                .join(median_cov)
                .join(ploidy_table)
                .map { batch_key, cohort, cutoffs_file, merged_pe, merged_bincov_file, merged_sr_file, median_cov_file, ploidy_table_file ->
                    tuple(
                        cohort,
                        batch_key,
                        selectSinglePath(cutoffs_file, "cutoffs"),
                        selectSinglePath(merged_pe, "merged_PE"),
                        selectSinglePath(merged_bincov_file, "merged_bincov"),
                        selectSinglePath(merged_sr_file, "merged_SR"),
                        selectSinglePath(median_cov_file, "median_cov"),
                        selectSinglePath(ploidy_table_file, "ploidy_table")
                    )
                }
                .join(GATKSV_MERGEBATCHSITES.out.merge_batch_sites_vcf)
                .map { cohort, batch_key, cutoffs_file, merged_pe, merged_bincov_file, merged_sr_file, median_cov_file, ploidy_table_file, merged_sites_vcf ->
                    tuple(
                        batch_key,
                        cohort,
                        cutoffs_file,
                        merged_pe,
                        merged_bincov_file,
                        merged_sr_file,
                        median_cov_file,
                        ploidy_table_file,
                        selectSinglePath(merged_sites_vcf, "merge_batch_sites_vcf")
                    )
                }

            GATKSV_GENOTYPEBATCH (
                genotype_batch_per_batch
            )
            versions = versions.mix(GATKSV_GENOTYPEBATCH.out.versions)
            genotype_batch_results = GATKSV_GENOTYPEBATCH.out.genotype_batch_results
            genotyped_vcf = GATKSV_GENOTYPEBATCH.out.genotyped_vcf
        }
    }

    emit:
    versions
    merge_batch_sites_results
    cohort_pesr_vcf
    cohort_depth_vcf
    merge_batch_sites_vcf
    genotype_batch_results
    genotyped_vcf
}
