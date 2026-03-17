//
// VCF_REFINEMENT
//

include { GATKSV_REFINECOMPLEXVARIANTS } from '../../modules/local/refine_complex_variants.nf'
include { GATKSV_JOINRAWCALLS } from '../../modules/local/join_raw_calls.nf'
include { GATKSV_SVCONCORDANCE } from '../../modules/local/sv_concordance.nf'
include { GATKSV_SCOREGENOTYPES } from '../../modules/local/score_genotypes.nf'
include { GATKSV_FILTERGENOTYPES } from '../../modules/local/filter_genotypes.nf'
include { GATKSV_ANNOTATEVCF } from '../../modules/local/annotate_vcf.nf'

workflow VCF_REFINEMENT {
    take:
    batch_cohort
    cleaned_vcf
    merged_PE
    merged_PE_index
    merged_dels
    merged_dups
    filtered_batch_samples_file
    clustered_depth_vcf
    clustered_depth_vcf_index
    clustered_manta_vcf
    clustered_manta_vcf_index
    clustered_wham_vcf
    clustered_wham_vcf_index
    clustered_scramble_vcf
    clustered_scramble_vcf_index
    ch_ped_file

    main:
    versions = Channel.empty()
    cpx_refined_vcf = Channel.empty()
    cpx_refined_vcf_index = Channel.empty()
    cpx_evidences = Channel.empty()
    joined_raw_calls_vcf = Channel.empty()
    joined_raw_calls_vcf_index = Channel.empty()
    join_raw_calls_ploidy_table = Channel.empty()
    concordance_vcf = Channel.empty()
    concordance_vcf_index = Channel.empty()
    unfiltered_recalibrated_vcf = Channel.empty()
    unfiltered_recalibrated_vcf_index = Channel.empty()
    vcf_optimization_table = Channel.empty()
    sl_cutoff_table = Channel.empty()
    sl_cutoff_qc_tarball = Channel.empty()
    filtered_vcf = Channel.empty()
    filtered_vcf_index = Channel.empty()
    main_vcf_qc_tarball = Channel.empty()
    annotated_vcf = Channel.empty()
    annotated_vcf_index = Channel.empty()

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
    def cohort_batches = batch_cohort
        .map { batch_key, cohort -> tuple(cohort.toString(), batch_key.toString()) }
        .groupTuple()
        .map { cohort, batch_keys ->
            tuple(cohort.toString(), asList(batch_keys).collect { it.toString() })
        }
    def groupBatchFilesByCohort = { ch, label ->
        batch_cohort
            .join(ch)
            .map { batch_key, cohort, file_value -> tuple(cohort.toString(), batch_key.toString(), selectSinglePath(file_value, label)) }
            .groupTuple()
            .join(cohort_batches)
            .map { cohort, _file_batch_keys, files, batch_ids ->
                tuple(cohort.toString(), orderFilesByKeys(batch_ids, files, label))
            }
    }

    if (params.run_refine_complex_variants && !params.run_clean_vcf) {
        throw new IllegalStateException("run_refine_complex_variants=true requires run_clean_vcf=true")
    }
    if (params.run_join_raw_calls && !params.run_cluster_batch) {
        throw new IllegalStateException("run_join_raw_calls=true requires run_cluster_batch=true")
    }
    if (params.run_sv_concordance && !params.run_refine_complex_variants) {
        throw new IllegalStateException("run_sv_concordance=true requires run_refine_complex_variants=true")
    }
    if (params.run_sv_concordance && !params.run_join_raw_calls) {
        throw new IllegalStateException("run_sv_concordance=true requires run_join_raw_calls=true")
    }
    if (params.run_score_genotypes && !params.run_sv_concordance) {
        throw new IllegalStateException("run_score_genotypes=true requires run_sv_concordance=true")
    }
    if (params.run_filter_genotypes && !params.run_score_genotypes) {
        throw new IllegalStateException("run_filter_genotypes=true requires run_score_genotypes=true")
    }
    if (params.run_filter_genotypes && !params.run_join_raw_calls) {
        throw new IllegalStateException("run_filter_genotypes=true requires run_join_raw_calls=true")
    }
    if (params.run_annotate_vcf && !params.run_filter_genotypes) {
        throw new IllegalStateException("run_annotate_vcf=true requires run_filter_genotypes=true")
    }

    if (params.run_refine_complex_variants) {
        refine_complex_variants_input = batch_cohort
            .join(filtered_batch_samples_file)
            .join(merged_PE)
            .join(merged_PE_index)
            .join(merged_dels)
            .join(merged_dups)
            .map { batch_key, cohort, batch_samples_file, pe_metric_file, pe_metric_index_file, depth_del_bed, depth_dup_bed ->
                tuple(
                    cohort,
                    batch_key.toString(),
                    selectSinglePath(batch_samples_file, "filtered_batch_samples_file"),
                    selectSinglePath(pe_metric_file, "merged_PE"),
                    selectSinglePath(pe_metric_index_file, "merged_PE_index"),
                    selectSinglePath(depth_del_bed, "merged_dels"),
                    selectSinglePath(depth_dup_bed, "merged_dups")
                )
            }
            .groupTuple()
            .join(cleaned_vcf)
            .map { cohort, batches, batch_sample_lists, pe_metric_files, pe_metric_index_files, depth_del_beds, depth_dup_beds, cleaned_vcf_file ->
                def batchIds = asList(batches).collect { it.toString() }
                tuple(
                    cohort,
                    selectSinglePath(cleaned_vcf_file, "cleaned_vcf"),
                    batchIds,
                    orderFilesByKeys(batchIds, batch_sample_lists, "filtered_batch_samples_file"),
                    orderFilesByKeys(batchIds, pe_metric_files, "merged_PE"),
                    orderFilesByKeys(batchIds, pe_metric_index_files, "merged_PE_index"),
                    orderFilesByKeys(batchIds, depth_del_beds, "merged_dels"),
                    orderFilesByKeys(batchIds, depth_dup_beds, "merged_dups")
                )
            }

        GATKSV_REFINECOMPLEXVARIANTS (
            refine_complex_variants_input
        )
        versions = versions.mix(GATKSV_REFINECOMPLEXVARIANTS.out.versions)
        cpx_refined_vcf = GATKSV_REFINECOMPLEXVARIANTS.out.cpx_refined_vcf
        cpx_refined_vcf_index = GATKSV_REFINECOMPLEXVARIANTS.out.cpx_refined_vcf_index
        cpx_evidences = GATKSV_REFINECOMPLEXVARIANTS.out.cpx_evidences
    }

    if (params.run_join_raw_calls) {
        join_raw_calls_input = groupBatchFilesByCohort(clustered_depth_vcf, "clustered_depth_vcf")
            .join(groupBatchFilesByCohort(clustered_depth_vcf_index, "clustered_depth_vcf_index"))
            .join(groupBatchFilesByCohort(clustered_manta_vcf, "clustered_manta_vcf"))
            .join(groupBatchFilesByCohort(clustered_manta_vcf_index, "clustered_manta_vcf_index"))
            .join(groupBatchFilesByCohort(clustered_wham_vcf, "clustered_wham_vcf"))
            .join(groupBatchFilesByCohort(clustered_wham_vcf_index, "clustered_wham_vcf_index"))
            .join(groupBatchFilesByCohort(clustered_scramble_vcf, "clustered_scramble_vcf"))
            .join(groupBatchFilesByCohort(clustered_scramble_vcf_index, "clustered_scramble_vcf_index"))
            .join(ch_ped_file)
            .map { cohort, depth_vcfs, depth_vcf_indexes, manta_vcfs, manta_vcf_indexes, wham_vcfs, wham_vcf_indexes, scramble_vcfs, scramble_vcf_indexes, ped_file ->
                tuple(
                    cohort,
                    ped_file,
                    depth_vcfs,
                    depth_vcf_indexes,
                    null,
                    null,
                    manta_vcfs,
                    manta_vcf_indexes,
                    null,
                    null,
                    scramble_vcfs,
                    scramble_vcf_indexes,
                    wham_vcfs,
                    wham_vcf_indexes
                )
            }

        GATKSV_JOINRAWCALLS (
            join_raw_calls_input
        )
        versions = versions.mix(GATKSV_JOINRAWCALLS.out.versions)
        joined_raw_calls_vcf = GATKSV_JOINRAWCALLS.out.joined_raw_calls_vcf
        joined_raw_calls_vcf_index = GATKSV_JOINRAWCALLS.out.joined_raw_calls_vcf_index
        join_raw_calls_ploidy_table = GATKSV_JOINRAWCALLS.out.ploidy_table
    }

    if (params.run_sv_concordance) {
        sv_concordance_input = cpx_refined_vcf
            .join(joined_raw_calls_vcf)
            .map { cohort, eval_vcf_file, truth_vcf_file ->
                tuple(
                    cohort,
                    selectSinglePath(eval_vcf_file, "cpx_refined_vcf"),
                    selectSinglePath(truth_vcf_file, "joined_raw_calls_vcf")
                )
            }

        GATKSV_SVCONCORDANCE (
            sv_concordance_input
        )
        versions = versions.mix(GATKSV_SVCONCORDANCE.out.versions)
        concordance_vcf = GATKSV_SVCONCORDANCE.out.concordance_vcf
        concordance_vcf_index = GATKSV_SVCONCORDANCE.out.concordance_vcf_index
    }


    if (params.run_score_genotypes) {
        score_genotypes_input = concordance_vcf
            .map { cohort, concordance_vcf_file ->
                tuple(
                    cohort,
                    selectSinglePath(concordance_vcf_file, "concordance_vcf"),
                    params.score_genotypes_truth_json
                )
            }

        GATKSV_SCOREGENOTYPES (
            score_genotypes_input
        )
        versions = versions.mix(GATKSV_SCOREGENOTYPES.out.versions)
        unfiltered_recalibrated_vcf = GATKSV_SCOREGENOTYPES.out.unfiltered_recalibrated_vcf
        unfiltered_recalibrated_vcf_index = GATKSV_SCOREGENOTYPES.out.unfiltered_recalibrated_vcf_index
        vcf_optimization_table = GATKSV_SCOREGENOTYPES.out.vcf_optimization_table
        sl_cutoff_table = GATKSV_SCOREGENOTYPES.out.sl_cutoff_table
        sl_cutoff_qc_tarball = GATKSV_SCOREGENOTYPES.out.sl_cutoff_qc_tarball
    }


    def optimized_sl_cutoff_by_cohort = sl_cutoff_table
        .map { cohort, cutoff_file -> tuple(cohort.toString(), selectSinglePath(cutoff_file, "sl_cutoff_table")) }
        .collect()
        .map { rows ->
            rows.collectEntries { row ->
                [(row[0].toString()): row[1]]
            }
        }

    if (params.run_filter_genotypes) {
        filter_genotypes_input = unfiltered_recalibrated_vcf
            .join(join_raw_calls_ploidy_table)
            .join(ch_ped_file)
            .combine(optimized_sl_cutoff_by_cohort)
            .map { cohort, unfiltered_vcf_file, ploidy_table_file, ped_file, optimizedCutoffMap ->
                tuple(
                    cohort,
                    selectSinglePath(unfiltered_vcf_file, "unfiltered_recalibrated_vcf"),
                    selectSinglePath(ploidy_table_file, "join_raw_calls_ploidy_table"),
                    selectSinglePath(ped_file, "ped_file"),
                    optimizedCutoffMap[cohort.toString()]
                )
            }

        GATKSV_FILTERGENOTYPES (
            filter_genotypes_input
        )
        versions = versions.mix(GATKSV_FILTERGENOTYPES.out.versions)
        filtered_vcf = GATKSV_FILTERGENOTYPES.out.filtered_vcf
        filtered_vcf_index = GATKSV_FILTERGENOTYPES.out.filtered_vcf_index
        main_vcf_qc_tarball = GATKSV_FILTERGENOTYPES.out.main_vcf_qc_tarball
    }


    if (params.run_annotate_vcf) {
        annotate_vcf_input = filtered_vcf
            .join(ch_ped_file)
            .map { cohort, filtered_vcf_file, ped_file ->
                tuple(
                    cohort,
                    selectSinglePath(filtered_vcf_file, "filtered_vcf"),
                    selectSinglePath(ped_file, "ped_file")
                )
            }

        GATKSV_ANNOTATEVCF (
            annotate_vcf_input
        )
        versions = versions.mix(GATKSV_ANNOTATEVCF.out.versions)
        annotated_vcf = GATKSV_ANNOTATEVCF.out.annotated_vcf
        annotated_vcf_index = GATKSV_ANNOTATEVCF.out.annotated_vcf_index
    }

    emit:
    versions
    cpx_refined_vcf
    cpx_refined_vcf_index
    cpx_evidences
    joined_raw_calls_vcf
    joined_raw_calls_vcf_index
    join_raw_calls_ploidy_table
    concordance_vcf
    concordance_vcf_index
    unfiltered_recalibrated_vcf
    unfiltered_recalibrated_vcf_index
    vcf_optimization_table
    sl_cutoff_table
    sl_cutoff_qc_tarball
    filtered_vcf
    filtered_vcf_index
    main_vcf_qc_tarball
    annotated_vcf
    annotated_vcf_index
}
