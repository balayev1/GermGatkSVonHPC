//
// JOINT_COHORT_CALLING
//

include { GATKSV_MERGEBATCHSITES } from '../../modules/local/merge_batch_sites.nf'
include { GATKSV_GENOTYPEBATCH } from '../../modules/local/genotype_batch.nf'
include { GATKSV_REGENOTYPECNVS } from '../../modules/local/regenotype_cnvs.nf'
include { GATKSV_COMBINEBATCHES } from '../../modules/local/combine_batches.nf'
include { GATKSV_RESOLVECOMPLEXVARIANTS } from '../../modules/local/resolve_complex_variants.nf'
include { GATKSV_GENOTYPECOMPLEXVARIANTS } from '../../modules/local/genotype_complex_variants.nf'
include { GATKSV_CLEANVCF } from '../../modules/local/clean_vcf.nf'

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
    merge_batch_sites_vcf_index = Channel.empty()
    genotyped_depth_vcf = Channel.empty()
    genotyped_depth_vcf_index = Channel.empty()
    genotyped_pesr_vcf = Channel.empty()
    genotyped_pesr_vcf_index = Channel.empty()
    genotyping_rd_table = Channel.empty()
    genotyping_pe_table = Channel.empty()
    genotyping_sr_table = Channel.empty()
    regeno_coverage_medians = Channel.empty()
    regenotyped_depth_vcfs = Channel.empty()
    regenotyped_depth_vcf_indexes = Channel.empty()
    number_regenotyped_file = Channel.empty()
    number_regenotyped_filtered_file = Channel.empty()
    combined_vcfs = Channel.empty()
    combined_vcf_indexes = Channel.empty()
    cluster_background_fail_lists = Channel.empty()
    cluster_bothside_pass_lists = Channel.empty()
    combine_batches_merged_vcf = Channel.empty()
    combine_batches_merged_vcf_index = Channel.empty()
    combine_batches_vcf = Channel.empty()
    combine_batches_vcf_index = Channel.empty()
    complex_resolve_vcfs = Channel.empty()
    complex_resolve_vcf_indexes = Channel.empty()
    complex_resolve_bothside_pass_list = Channel.empty()
    complex_resolve_background_fail_list = Channel.empty()
    breakpoint_overlap_dropped_record_vcfs = Channel.empty()
    breakpoint_overlap_dropped_record_vcf_indexes = Channel.empty()
    complex_resolve_merged_vcf = Channel.empty()
    complex_resolve_merged_vcf_index = Channel.empty()
    complex_genotype_vcfs = Channel.empty()
    complex_genotype_vcf_indexes = Channel.empty()
    complex_genotype_merged_vcf = Channel.empty()
    complex_genotype_merged_vcf_index = Channel.empty()
    cleaned_vcf = Channel.empty()
    cleaned_vcf_index = Channel.empty()
    metrics_file_makecohortvcf = Channel.empty()

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
    def asPathList = { value -> asList(value).collect { file(it.toString()) } }
    def genome_resources = params.genomes?.get(params.genome as String) ?: [:]
    def contigOrder = genome_resources.primary_contigs_fai ? file(genome_resources.primary_contigs_fai.toString(), checkIfExists: true).readLines().collect { it.tokenize('\t')[0] } : []
    def matchesKey = { pathValue, key ->
        def name = pathValue.getName()
        def token = key.toString()
        name == token ||
            name.startsWith("${token}.") ||
            name.startsWith("${token}_") ||
            name.contains(".${token}.") ||
            name.contains(".${token}_") ||
            name.contains("_${token}.") ||
            name.contains("_${token}_") ||
            name.endsWith(".${token}")
    }
    def orderFilesByKeys = { keys, value, label ->
        def files = asPathList(value)
        def expectedKeys = asList(keys).collect { it.toString() }
        if (files.size() != expectedKeys.size()) {
            throw new IllegalStateException("Expected ${expectedKeys.size()} ${label} file(s) but found ${files.size()}")
        }
        expectedKeys.collect { key ->
            def matches = files.findAll { matchesKey(it, key) }
            if (matches.size() != 1) {
                throw new IllegalStateException("Could not uniquely match ${label} for ${key}")
            }
            matches[0]
        }
    }
    def orderContigFiles = { value, label ->
        if (!contigOrder) {
            return asPathList(value)
        }
        orderFilesByKeys(contigOrder, value, label)
    }

    if (params.run_regenotype_cnvs && !params.run_merge_batch_sites) {
        throw new IllegalStateException("run_regenotype_cnvs=true requires run_merge_batch_sites=true")
    }
    if (params.run_resolve_complex_variants && !params.run_combine_batches) {
        throw new IllegalStateException("run_resolve_complex_variants=true requires run_combine_batches=true")
    }
    if (params.run_genotype_complex_variants && !params.run_resolve_complex_variants) {
        throw new IllegalStateException("run_genotype_complex_variants=true requires run_resolve_complex_variants=true")
    }
    if (params.run_clean_vcf && !params.run_genotype_complex_variants) {
        throw new IllegalStateException("run_clean_vcf=true requires run_genotype_complex_variants=true")
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
        merge_batch_sites_vcf_index = GATKSV_MERGEBATCHSITES.out.merge_batch_sites_vcf_index

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
            genotyped_depth_vcf_index = GATKSV_GENOTYPEBATCH.out.genotyped_depth_vcf_index
            genotyped_pesr_vcf = GATKSV_GENOTYPEBATCH.out.genotyped_pesr_vcf
            genotyped_pesr_vcf_index = GATKSV_GENOTYPEBATCH.out.genotyped_pesr_vcf_index
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
            regenotyped_depth_vcf_indexes = GATKSV_REGENOTYPECNVS.out.regenotyped_depth_vcf_indexes
            number_regenotyped_file = GATKSV_REGENOTYPECNVS.out.number_regenotyped_file
            number_regenotyped_filtered_file = GATKSV_REGENOTYPECNVS.out.number_regenotyped_filtered_file
        }

        if (params.run_combine_batches) {
            if (!params.run_regenotype_cnvs) {
                throw new IllegalStateException("run_combine_batches=true requires run_regenotype_cnvs=true")
            }

            combine_batches_input = batch_cohort
                .join(genotyped_pesr_vcf)
                .map { batch_key, cohort, pesr_vcf ->
                    tuple(
                        cohort,
                        batch_key.toString(),
                        selectSinglePath(pesr_vcf, "pesr_vcf")
                    )
                }
                .groupTuple()
                .join(regenotyped_depth_vcfs)
                .join(ch_ped_file)
                .map { cohort, batches, pesr_vcfs, depth_vcfs, ped_file ->
                    def batchIds = asList(batches).collect { it.toString() }
                    tuple(
                        cohort,
                        batchIds,
                        ped_file,
                        orderFilesByKeys(batchIds, pesr_vcfs, "genotyped_pesr_vcf"),
                        orderFilesByKeys(batchIds, depth_vcfs, "regenotyped_depth_vcf")
                    )
                }

            GATKSV_COMBINEBATCHES (
                combine_batches_input
            )
            versions = versions.mix(GATKSV_COMBINEBATCHES.out.versions)
            combined_vcfs = GATKSV_COMBINEBATCHES.out.combined_vcfs
            combined_vcf_indexes = GATKSV_COMBINEBATCHES.out.combined_vcf_indexes
            cluster_background_fail_lists = GATKSV_COMBINEBATCHES.out.cluster_background_fail_lists
            cluster_bothside_pass_lists = GATKSV_COMBINEBATCHES.out.cluster_bothside_pass_lists
            combine_batches_merged_vcf = GATKSV_COMBINEBATCHES.out.combine_batches_merged_vcf
            combine_batches_merged_vcf_index = GATKSV_COMBINEBATCHES.out.combine_batches_merged_vcf_index
            combine_batches_vcf = combine_batches_merged_vcf
            combine_batches_vcf_index = combine_batches_merged_vcf_index

            if (params.run_resolve_complex_variants) {
                resolve_complex_variants_input = batch_cohort
                    .join(merged_PE)
                    .join(cutoffs)
                    .map { batch_key, cohort, merged_pe, cutoff_file ->
                        tuple(
                            cohort,
                            batch_key.toString(),
                            selectSinglePath(merged_pe, "merged_PE"),
                            selectSinglePath(cutoff_file, "cutoffs")
                        )
                    }
                    .groupTuple()
                    .join(combined_vcfs)
                    .join(cluster_bothside_pass_lists)
                    .join(cluster_background_fail_lists)
                    .map { cohort, batches, disc_files, rf_cutoff_files, cluster_vcf_files, bothside_pass_lists, background_fail_lists ->
                        def batchIds = asList(batches).collect { it.toString() }
                        tuple(
                            cohort,
                            orderContigFiles(cluster_vcf_files, "combined_vcf"),
                            orderContigFiles(bothside_pass_lists, "cluster_bothside_pass_list"),
                            orderContigFiles(background_fail_lists, "cluster_background_fail_list"),
                            orderFilesByKeys(batchIds, disc_files, "merged_PE"),
                            orderFilesByKeys(batchIds, rf_cutoff_files, "cutoffs")
                        )
                    }

                GATKSV_RESOLVECOMPLEXVARIANTS (
                    resolve_complex_variants_input
                )
                versions = versions.mix(GATKSV_RESOLVECOMPLEXVARIANTS.out.versions)
                complex_resolve_vcfs = GATKSV_RESOLVECOMPLEXVARIANTS.out.complex_resolve_vcfs
                complex_resolve_vcf_indexes = GATKSV_RESOLVECOMPLEXVARIANTS.out.complex_resolve_vcf_indexes
                complex_resolve_bothside_pass_list = GATKSV_RESOLVECOMPLEXVARIANTS.out.complex_resolve_bothside_pass_list
                complex_resolve_background_fail_list = GATKSV_RESOLVECOMPLEXVARIANTS.out.complex_resolve_background_fail_list
                breakpoint_overlap_dropped_record_vcfs = GATKSV_RESOLVECOMPLEXVARIANTS.out.breakpoint_overlap_dropped_record_vcfs
                breakpoint_overlap_dropped_record_vcf_indexes = GATKSV_RESOLVECOMPLEXVARIANTS.out.breakpoint_overlap_dropped_record_vcf_indexes
                complex_resolve_merged_vcf = GATKSV_RESOLVECOMPLEXVARIANTS.out.complex_resolve_merged_vcf
                complex_resolve_merged_vcf_index = GATKSV_RESOLVECOMPLEXVARIANTS.out.complex_resolve_merged_vcf_index
            }

            if (params.run_genotype_complex_variants) {
                genotype_complex_variants_input = batch_cohort
                    .join(merged_bincov)
                    .join(genotyping_rd_table)
                    .join(median_cov)
                    .map { batch_key, cohort, bincov_file, rd_table, median_cov_file ->
                        tuple(
                            cohort,
                            batch_key.toString(),
                            selectSinglePath(bincov_file, "merged_bincov"),
                            selectSinglePath(rd_table, "genotyping_rd_table"),
                            selectSinglePath(median_cov_file, "median_cov")
                        )
                    }
                    .groupTuple()
                    .join(regenotyped_depth_vcfs)
                    .join(complex_resolve_vcfs)
                    .join(complex_resolve_vcf_indexes)
                    .join(ch_ped_file)
                    .map { cohort, batches, bincov_files, genotyping_rd_tables, median_coverage_files, depth_vcfs, complex_resolve_vcfs_files, complex_resolve_vcf_index_files, ped_file ->
                        def batchIds = asList(batches).collect { it.toString() }
                        tuple(
                            cohort,
                            batchIds,
                            orderFilesByKeys(batchIds, depth_vcfs, "regenotyped_depth_vcf"),
                            orderContigFiles(complex_resolve_vcfs_files, "complex_resolve_vcf"),
                            orderContigFiles(complex_resolve_vcf_index_files, "complex_resolve_vcf_index"),
                            ped_file,
                            orderFilesByKeys(batchIds, bincov_files, "merged_bincov"),
                            orderFilesByKeys(batchIds, genotyping_rd_tables, "genotyping_rd_table"),
                            orderFilesByKeys(batchIds, median_coverage_files, "median_cov")
                        )
                    }

                GATKSV_GENOTYPECOMPLEXVARIANTS (
                    genotype_complex_variants_input
                )
                versions = versions.mix(GATKSV_GENOTYPECOMPLEXVARIANTS.out.versions)
                complex_genotype_vcfs = GATKSV_GENOTYPECOMPLEXVARIANTS.out.complex_genotype_vcfs
                complex_genotype_vcf_indexes = GATKSV_GENOTYPECOMPLEXVARIANTS.out.complex_genotype_vcf_indexes
                complex_genotype_merged_vcf = GATKSV_GENOTYPECOMPLEXVARIANTS.out.complex_genotype_merged_vcf
                complex_genotype_merged_vcf_index = GATKSV_GENOTYPECOMPLEXVARIANTS.out.complex_genotype_merged_vcf_index
            }

            if (params.run_clean_vcf) {
                clean_vcf_input = complex_genotype_vcfs
                    .join(complex_resolve_bothside_pass_list)
                    .join(complex_resolve_background_fail_list)
                    .join(ch_ped_file)
                    .map { cohort, complex_genotype_vcf_files, bothside_pass_file, background_fail_file, ped_file ->
                        tuple(
                            cohort,
                            orderContigFiles(complex_genotype_vcf_files, "complex_genotype_vcf"),
                            selectSinglePath(bothside_pass_file, "complex_resolve_bothside_pass_list"),
                            selectSinglePath(background_fail_file, "complex_resolve_background_fail_list"),
                            ped_file
                        )
                    }

                GATKSV_CLEANVCF (
                    clean_vcf_input
                )
                versions = versions.mix(GATKSV_CLEANVCF.out.versions)
                cleaned_vcf = GATKSV_CLEANVCF.out.cleaned_vcf
                cleaned_vcf_index = GATKSV_CLEANVCF.out.cleaned_vcf_index
                metrics_file_makecohortvcf = GATKSV_CLEANVCF.out.metrics_file_makecohortvcf
            }
        }
    }

    emit:
    versions
    merge_batch_sites_vcf
    merge_batch_sites_vcf_index
    genotyped_depth_vcf
    genotyped_depth_vcf_index
    genotyped_pesr_vcf
    genotyped_pesr_vcf_index
    genotyping_rd_table
    genotyping_pe_table
    genotyping_sr_table
    regeno_coverage_medians
    regenotyped_depth_vcfs
    regenotyped_depth_vcf_indexes
    number_regenotyped_file
    number_regenotyped_filtered_file
    combined_vcfs
    combined_vcf_indexes
    cluster_background_fail_lists
    cluster_bothside_pass_lists
    combine_batches_merged_vcf
    combine_batches_merged_vcf_index
    combine_batches_vcf
    combine_batches_vcf_index
    complex_resolve_vcfs
    complex_resolve_vcf_indexes
    complex_resolve_bothside_pass_list
    complex_resolve_background_fail_list
    breakpoint_overlap_dropped_record_vcfs
    breakpoint_overlap_dropped_record_vcf_indexes
    complex_resolve_merged_vcf
    complex_resolve_merged_vcf_index
    complex_genotype_vcfs
    complex_genotype_vcf_indexes
    complex_genotype_merged_vcf
    complex_genotype_merged_vcf_index
    cleaned_vcf
    cleaned_vcf_index
    metrics_file_makecohortvcf
}
