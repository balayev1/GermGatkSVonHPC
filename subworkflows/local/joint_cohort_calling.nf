//
// JOINT_COHORT_CALLING
//

include { GATKSV_MERGEBATCHSITES } from '../../modules/local/merge_batch_sites.nf'
include { GATKSV_GENOTYPEBATCH } from '../../modules/local/genotype_batch.nf'
include { GATKSV_REGENOTYPECNVS } from '../../modules/local/regenotype_cnvs.nf'
include { GATKSV_COMBINEBATCHES } from '../../modules/local/combine_batches.nf'

workflow JOINT_COHORT_CALLING {
    take:
    batch_cohort
    filtered_depth_vcf
    filtered_pesr_vcf
    ploidy_table
    cutoffs
    merged_PE
    merged_bincov
    merged_bincov_index
    merged_SR
    median_cov
    ch_ped_file

    main:
    versions = Channel.empty()

    merge_batch_sites_vcf = Channel.empty()
    genotyped_depth_vcf = Channel.empty()
    genotyped_pesr_vcf = Channel.empty()
    genotyping_rd_table = Channel.empty()
    genotyping_pe_table = Channel.empty()
    genotyping_sr_table = Channel.empty()
    regeno_coverage_medians = Channel.empty()
    regenotyped_depth_vcfs = Channel.empty()
    combine_batches_vcf = Channel.empty()

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

    if (params.run_regenotype_cnvs && !params.run_merge_batch_sites) {
        throw new IllegalStateException("run_regenotype_cnvs=true requires run_merge_batch_sites=true")
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
                .join(merge_batch_sites_vcf)
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
            genotyped_depth_vcf = GATKSV_GENOTYPEBATCH.out.genotyped_depth_vcf
            genotyped_pesr_vcf = GATKSV_GENOTYPEBATCH.out.genotyped_pesr_vcf
            genotyping_rd_table = GATKSV_GENOTYPEBATCH.out.genotyping_rd_table
            genotyping_pe_table = GATKSV_GENOTYPEBATCH.out.genotyping_pe_table
            genotyping_sr_table = GATKSV_GENOTYPEBATCH.out.genotyping_sr_table
            regeno_coverage_medians = GATKSV_GENOTYPEBATCH.out.regeno_coverage_medians
        }

        if (params.run_regenotype_cnvs) {
            if (!params.run_genotype_batch) {
                throw new IllegalStateException("run_regenotype_cnvs=true requires run_genotype_batch=true")
            }

            regenotype_cnvs_input = batch_cohort
                .join(filtered_depth_vcf)
                .join(genotyped_depth_vcf)
                .join(merged_bincov)
                .join(merged_bincov_index)
                .join(median_cov)
                .join(genotyping_rd_table)
                .join(ploidy_table)
                .join(regeno_coverage_medians)
                .map { batch_key, cohort, batch_depth_vcf, depth_vcf, coveragefile, coveragefile_idx, medianfile, rd_table, ploidy_table_file, regeno_cov_median ->
                    tuple(
                        cohort,
                        batch_key.toString(),
                        selectSinglePath(batch_depth_vcf, "batch_depth_vcf"),
                        selectSinglePath(depth_vcf, "depth_vcf"),
                        selectSinglePath(coveragefile, "coveragefile"),
                        selectSinglePath(coveragefile_idx, "coveragefile_idx"),
                        selectSinglePath(medianfile, "medianfile"),
                        selectSinglePath(rd_table, "genotyping_rd_table"),
                        selectSinglePath(ploidy_table_file, "ploidy_table"),
                        selectSinglePath(regeno_cov_median, "regeno_coverage_medians")
                    )
                }
                .groupTuple()
                .join(merge_batch_sites_vcf)
                .map { cohort, batches, batch_depth_vcfs, depth_vcfs, coveragefiles, coveragefile_idxs, medianfiles, rd_tables, ploidy_tables, regeno_cov_medians, merge_sites_vcf ->
                    tuple(
                        cohort,
                        asList(depth_vcfs).collect { file(it.toString()) },
                        selectSinglePath(merge_sites_vcf, "merge_batch_sites_vcf"),
                        asList(batch_depth_vcfs).collect { file(it.toString()) },
                        asList(coveragefiles).collect { file(it.toString()) },
                        asList(coveragefile_idxs).collect { file(it.toString()) },
                        asList(medianfiles).collect { file(it.toString()) },
                        asList(rd_tables).collect { file(it.toString()) },
                        asList(ploidy_tables).collect { file(it.toString()) },
                        asList(batches).collect { it.toString() },
                        asList(regeno_cov_medians).collect { file(it.toString()) }
                    )
                }

            GATKSV_REGENOTYPECNVS (
                regenotype_cnvs_input
            )
            versions = versions.mix(GATKSV_REGENOTYPECNVS.out.versions)
            regenotyped_depth_vcfs = GATKSV_REGENOTYPECNVS.out.regenotyped_depth_vcfs
        }

        if (params.run_combine_batches) {
            if (!params.run_regenotype_cnvs) {
                throw new IllegalStateException("run_combine_batches=true requires run_regenotype_cnvs=true")
            }

            combine_batches_input = batch_cohort
                .join(genotyped_pesr_vcf)
                .join(regenotyped_depth_vcfs)
                .map { batch_key, cohort, pesr_vcf, depth_vcf ->
                    tuple(
                        cohort,
                        batch_key.toString(),
                        selectSinglePath(pesr_vcf, "pesr_vcf"),
                        selectSinglePath(depth_vcf, "depth_vcf")
                    )
                }
                .groupTuple()
                .join(ch_ped_file)
                .map { cohort, batches, pesr_vcfs, depth_vcfs, ped_file ->
                    tuple(
                        cohort,
                        asList(batches).collect { it.toString() },
                        ped_file,
                        asList(pesr_vcfs).collect { file(it.toString()) },
                        asList(depth_vcfs).collect { file(it.toString()) }
                    )
                }

            GATKSV_COMBINEBATCHES (
                combine_batches_input
            )
            versions = versions.mix(GATKSV_COMBINEBATCHES.out.versions)
            combine_batches_vcf = GATKSV_COMBINEBATCHES.out.combined_vcf
        }
    }

    emit:
    versions
    merge_batch_sites_vcf
    genotyped_depth_vcf
    genotyped_pesr_vcf
    genotyping_rd_table
    genotyping_pe_table
    genotyping_sr_table
    regeno_coverage_medians
    regenotyped_depth_vcfs
    combine_batches_vcf
}
