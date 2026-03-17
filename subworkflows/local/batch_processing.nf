//
// BATCH_PROCESSING
//

include { GATK_TRAINGCNV } from '../../modules/local/train_gcnv.nf'
include { GATKSV_GATHERBATCHEVIDENCE } from '../../modules/local/gather_batch_evidence.nf'
include { GATKSV_CLUSTERBATCH } from '../../modules/local/cluster_batch.nf'
include { GATKSV_GENERATEBATCHMETRICS } from '../../modules/local/generate_batch_metrics.nf'
include { GATKSV_FILTERBATCHSITES } from '../../modules/local/filter_batch_sites.nf'
include { GATKSV_FILTERBATCHSAMPLES } from '../../modules/local/filter_batch_samples.nf'

workflow BATCH_PROCESSING {
    take:
    train_gcnv_input  // channel: [ cohort_or_batch, sample_ids, count_files, outlier_sample_ids ]
    batch_evidence_input  // channel: [ batch, cohort, sample_ids, ped_file, counts, PE_files, SR_files, SD_files, *_vcfs ]
    outlier_sample_ids_by_batch  // channel: [ batch, outlier_sample_ids ]

    main:
    versions = Channel.empty()
    gcnv_model_tars = Channel.empty()
    contig_ploidy_model_tar = Channel.empty()

    if (params.run_train_gcnv) {
        GATK_TRAINGCNV (
            train_gcnv_input
        )
        versions = versions.mix(GATK_TRAINGCNV.out.versions)
        gcnv_model_tars = GATK_TRAINGCNV.out.cohort_gcnv_model_tars
        contig_ploidy_model_tar = GATK_TRAINGCNV.out.cohort_contig_ploidy_model_tar
    }

    generate_batch_metrics_metrics = Channel.empty()
    metrics_file_batchmetrics = Channel.empty()
    generate_batch_metrics_ploidy_table = Channel.empty()
    batch_cohort = Channel.empty()
    merged_BAF = Channel.empty()
    merged_BAF_index = Channel.empty()
    merged_PE = Channel.empty()
    merged_PE_index = Channel.empty()
    merged_bincov = Channel.empty()
    merged_bincov_index = Channel.empty()
    merged_SR = Channel.empty()
    merged_SR_index = Channel.empty()
    merged_dels = Channel.empty()
    merged_dups = Channel.empty()
    cnmops_del = Channel.empty()
    cnmops_del_index = Channel.empty()
    cnmops_dup = Channel.empty()
    cnmops_dup_index = Channel.empty()
    cnmops_large_del = Channel.empty()
    cnmops_large_del_index = Channel.empty()
    cnmops_large_dup = Channel.empty()
    cnmops_large_dup_index = Channel.empty()
    median_cov = Channel.empty()
    std_manta_vcf_tar = Channel.empty()
    std_wham_vcf_tar = Channel.empty()
    std_scramble_vcf_tar = Channel.empty()
    batch_ploidy_matrix = Channel.empty()
    batch_ploidy_plots = Channel.empty()
    BAF_stats = Channel.empty()
    RD_stats = Channel.empty()
    PE_stats = Channel.empty()
    SR_stats = Channel.empty()
    Matrix_QC_plot = Channel.empty()
    metrics_file_batchevidence = Channel.empty()
    manta_tloc = Channel.empty()
    filter_batch_sites_cutoffs = Channel.empty()
    filter_batch_sites_scores = Channel.empty()
    filter_batch_sites_RF_intermediate_files = Channel.empty()
    filter_batch_sites_sv_counts = Channel.empty()
    filter_batch_sites_sv_count_plots = Channel.empty()
    sites_filtered_depth_vcf = Channel.empty()
    sites_filtered_manta_vcf = Channel.empty()
    sites_filtered_wham_vcf = Channel.empty()
    sites_filtered_scramble_vcf = Channel.empty()
    filter_batch_sites_outlier_samples_preview = Channel.empty()
    filter_batch_sites_outlier_samples_with_reason = Channel.empty()
    filter_batch_sites_num_outlier_samples = Channel.empty()
    outlier_filtered_depth_vcf = Channel.empty()
    outlier_filtered_depth_vcf_index = Channel.empty()
    outlier_filtered_manta_vcf = Channel.empty()
    outlier_filtered_manta_vcf_index = Channel.empty()
    outlier_filtered_scramble_vcf = Channel.empty()
    outlier_filtered_scramble_vcf_index = Channel.empty()
    outlier_filtered_wham_vcf = Channel.empty()
    outlier_filtered_wham_vcf_index = Channel.empty()
    outlier_filtered_pesr_vcf = Channel.empty()
    outlier_filtered_pesr_vcf_index = Channel.empty()
    clustered_depth_vcf = Channel.empty()
    clustered_depth_vcf_index = Channel.empty()
    clustered_manta_vcf = Channel.empty()
    clustered_manta_vcf_index = Channel.empty()
    clustered_wham_vcf = Channel.empty()
    clustered_wham_vcf_index = Channel.empty()
    clustered_scramble_vcf = Channel.empty()
    clustered_scramble_vcf_index = Channel.empty()
    clustered_sv_counts = Channel.empty()
    clustered_sv_count_plots = Channel.empty()
    clustered_outlier_samples_preview = Channel.empty()
    clustered_outlier_samples_with_reason = Channel.empty()
    cluster_batch_num_outlier_samples = Channel.empty()
    metrics_file_clusterbatch = Channel.empty()
    cluster_batch_num_outlier_samples_int = Channel.empty()
    filter_batch_sites_num_outlier_samples_int = Channel.empty()
    filtered_depth_vcf = Channel.empty()
    filtered_pesr_vcf = Channel.empty()
    outlier_samples_excluded_file = Channel.empty()
    filtered_batch_samples_file = Channel.empty()
    if (params.run_gather_batch_evidence) {
        if (!params.run_train_gcnv) {
            throw new IllegalStateException("run_gather_batch_evidence=true requires run_train_gcnv=true because GatherBatchEvidence consumes TrainGCNV model outputs")
        }
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
            def nonIndexCandidates = candidates.findAll { f ->
                def n = f.getName().toLowerCase()
                !(n.endsWith('.tbi') || n.endsWith('.csi') || n.endsWith('.idx'))
            }
            def pool = nonIndexCandidates ? nonIndexCandidates : candidates
            if (pool.size() == 1) {
                return pool[0]
            }
            def merged = pool.findAll { it.getName().toLowerCase().contains("merged") }
            if (merged) {
                return merged[0]
            }
            def vcf = pool.findAll { f ->
                def n = f.getName().toLowerCase()
                n.endsWith('.vcf') || n.endsWith('.vcf.gz')
            }
            return vcf ? vcf[0] : pool[0]
        }

        batch_cohort = batch_evidence_input
            .map { batch_key, cohort, sample_ids, ped_file, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files ->
                tuple(batch_key, cohort.toString())
            }
            .unique()

        batch_evidence_with_models = batch_evidence_input
            .join(gcnv_model_tars)
            .map { batch_key, cohort, sample_ids, ped_file, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files, gcnv_model_tars_files ->
                tuple(
                    batch_key,
                    sample_ids,
                    ped_file,
                    counts_files,
                    pe_files,
                    sr_files,
                    sd_files,
                    manta_files,
                    wham_files,
                    scramble_files,
                    gcnv_model_tars_files
                )
            }
            .join(contig_ploidy_model_tar)
            .map { batch_key, sample_ids, ped_file, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files, gcnv_model_tars_files, contig_ploidy_model_tar_file ->
                tuple(
                    batch_key,
                    asList(sample_ids).collect { it.toString() },
                    ped_file,
                    asList(counts_files),
                    asList(pe_files),
                    asList(sr_files),
                    asList(sd_files),
                    asList(manta_files),
                    asList(wham_files),
                    asList(scramble_files),
                    asList(gcnv_model_tars_files),
                    contig_ploidy_model_tar_file
                )
            }

        GATKSV_GATHERBATCHEVIDENCE (
            batch_evidence_with_models
        )
        versions = versions.mix(GATKSV_GATHERBATCHEVIDENCE.out.versions)
        merged_BAF = GATKSV_GATHERBATCHEVIDENCE.out.merged_BAF
        merged_BAF_index = GATKSV_GATHERBATCHEVIDENCE.out.merged_BAF_index
        merged_PE = GATKSV_GATHERBATCHEVIDENCE.out.merged_PE
        merged_PE_index = GATKSV_GATHERBATCHEVIDENCE.out.merged_PE_index
        merged_bincov = GATKSV_GATHERBATCHEVIDENCE.out.merged_bincov
        merged_bincov_index = GATKSV_GATHERBATCHEVIDENCE.out.merged_bincov_index
        merged_SR = GATKSV_GATHERBATCHEVIDENCE.out.merged_SR
        merged_SR_index = GATKSV_GATHERBATCHEVIDENCE.out.merged_SR_index
        merged_dels = GATKSV_GATHERBATCHEVIDENCE.out.merged_dels
        merged_dups = GATKSV_GATHERBATCHEVIDENCE.out.merged_dups
        cnmops_del = GATKSV_GATHERBATCHEVIDENCE.out.cnmops_del
        cnmops_del_index = GATKSV_GATHERBATCHEVIDENCE.out.cnmops_del_index
        cnmops_dup = GATKSV_GATHERBATCHEVIDENCE.out.cnmops_dup
        cnmops_dup_index = GATKSV_GATHERBATCHEVIDENCE.out.cnmops_dup_index
        cnmops_large_del = GATKSV_GATHERBATCHEVIDENCE.out.cnmops_large_del
        cnmops_large_del_index = GATKSV_GATHERBATCHEVIDENCE.out.cnmops_large_del_index
        cnmops_large_dup = GATKSV_GATHERBATCHEVIDENCE.out.cnmops_large_dup
        cnmops_large_dup_index = GATKSV_GATHERBATCHEVIDENCE.out.cnmops_large_dup_index
        median_cov = GATKSV_GATHERBATCHEVIDENCE.out.median_cov
        std_manta_vcf_tar = GATKSV_GATHERBATCHEVIDENCE.out.std_manta_vcf_tar
        std_wham_vcf_tar = GATKSV_GATHERBATCHEVIDENCE.out.std_wham_vcf_tar
        std_scramble_vcf_tar = GATKSV_GATHERBATCHEVIDENCE.out.std_scramble_vcf_tar
        batch_ploidy_matrix = GATKSV_GATHERBATCHEVIDENCE.out.batch_ploidy_matrix
        batch_ploidy_plots = GATKSV_GATHERBATCHEVIDENCE.out.batch_ploidy_plots
        BAF_stats = GATKSV_GATHERBATCHEVIDENCE.out.BAF_stats
        RD_stats = GATKSV_GATHERBATCHEVIDENCE.out.RD_stats
        PE_stats = GATKSV_GATHERBATCHEVIDENCE.out.PE_stats
        SR_stats = GATKSV_GATHERBATCHEVIDENCE.out.SR_stats
        Matrix_QC_plot = GATKSV_GATHERBATCHEVIDENCE.out.Matrix_QC_plot
        metrics_file_batchevidence = GATKSV_GATHERBATCHEVIDENCE.out.metrics_file_batchevidence
        manta_tloc = GATKSV_GATHERBATCHEVIDENCE.out.manta_tloc

        if (params.run_cluster_batch) {
            def nIqrCutoffPlotting = params.cluster_n_iqr_cutoff_plotting
            ped_by_batch = batch_evidence_input
                .map { batch_key, cohort, sample_ids, ped_file, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files ->
                    tuple(batch_key, cohort.toString(), ped_file)
                }

            cluster_batch_input = ped_by_batch
                .join(merged_dels)
                .join(merged_dups)
                .join(std_wham_vcf_tar)
                .join(std_manta_vcf_tar)
                .join(std_scramble_vcf_tar)
                .map { batch_key, cohort, ped_file, merged_dels, merged_dups, std_wham_vcf_tar, std_manta_vcf_tar, std_scramble_vcf_tar ->
                    tuple(
                        batch_key,
                        cohort,
                        ped_file,
                        selectSinglePath(merged_dels, "merged_dels"),
                        selectSinglePath(merged_dups, "merged_dups"),
                        selectSinglePath(std_wham_vcf_tar, "std_wham_vcf_tar"),
                        selectSinglePath(std_manta_vcf_tar, "std_manta_vcf_tar"),
                        selectSinglePath(std_scramble_vcf_tar, "std_scramble_vcf_tar"),
                        nIqrCutoffPlotting
                    )
                }

            GATKSV_CLUSTERBATCH (
                cluster_batch_input
            )
            versions = versions.mix(GATKSV_CLUSTERBATCH.out.versions)
            clustered_depth_vcf = GATKSV_CLUSTERBATCH.out.clustered_depth_vcf
            clustered_depth_vcf_index = GATKSV_CLUSTERBATCH.out.clustered_depth_vcf_index
            clustered_manta_vcf = GATKSV_CLUSTERBATCH.out.clustered_manta_vcf
            clustered_manta_vcf_index = GATKSV_CLUSTERBATCH.out.clustered_manta_vcf_index
            clustered_wham_vcf = GATKSV_CLUSTERBATCH.out.clustered_wham_vcf
            clustered_wham_vcf_index = GATKSV_CLUSTERBATCH.out.clustered_wham_vcf_index
            clustered_scramble_vcf = GATKSV_CLUSTERBATCH.out.clustered_scramble_vcf
            clustered_scramble_vcf_index = GATKSV_CLUSTERBATCH.out.clustered_scramble_vcf_index
            clustered_sv_counts = GATKSV_CLUSTERBATCH.out.clustered_sv_counts
            clustered_sv_count_plots = GATKSV_CLUSTERBATCH.out.clustered_sv_count_plots
            clustered_outlier_samples_preview = GATKSV_CLUSTERBATCH.out.clustered_outlier_samples_preview
            clustered_outlier_samples_with_reason = GATKSV_CLUSTERBATCH.out.clustered_outlier_samples_with_reason
            cluster_batch_num_outlier_samples = GATKSV_CLUSTERBATCH.out.num_outlier_samples
            metrics_file_clusterbatch = GATKSV_CLUSTERBATCH.out.metrics_file_clusterbatch
            cluster_batch_num_outlier_samples_int = cluster_batch_num_outlier_samples
                .map { batch_key, outliers_file ->
                    tuple(batch_key, outliers_file.text.trim().toInteger())
                }

            if (params.run_generate_batch_metrics) {
                generate_batch_metrics_input = ped_by_batch
                    .join(merged_PE)
                    .join(merged_BAF)
                    .join(merged_bincov)
                    .join(merged_SR)
                    .join(median_cov)
                    .join(clustered_depth_vcf)
                    .join(clustered_manta_vcf)
                    .join(clustered_wham_vcf)
                    .join(clustered_scramble_vcf)
                    .join(outlier_sample_ids_by_batch)
                    .map { batch_key, cohort, ped_file, merged_pe, merged_baf, merged_bincov_file, merged_sr, median_cov_file, clustered_depth_vcf, clustered_manta_vcf, clustered_wham_vcf, clustered_scramble_vcf, outlier_sample_ids_file ->
                        tuple(
                            batch_key,
                            ped_file,
                            selectSinglePath(merged_pe, "merged_PE"),
                            selectSinglePath(merged_baf, "merged_BAF"),
                            selectSinglePath(merged_bincov_file, "merged_bincov"),
                            selectSinglePath(merged_sr, "merged_SR"),
                            selectSinglePath(median_cov_file, "median_cov"),
                            selectSinglePath(clustered_depth_vcf, "clustered_depth_vcf"),
                            selectSinglePath(clustered_manta_vcf, "clustered_manta_vcf"),
                            selectSinglePath(clustered_wham_vcf, "clustered_wham_vcf"),
                            selectSinglePath(clustered_scramble_vcf, "clustered_scramble_vcf"),
                            selectSinglePath(outlier_sample_ids_file, "outlier_sample_ids")
                        )
                    }

                GATKSV_GENERATEBATCHMETRICS (
                    generate_batch_metrics_input
                )
                versions = versions.mix(GATKSV_GENERATEBATCHMETRICS.out.versions)
                generate_batch_metrics_metrics = GATKSV_GENERATEBATCHMETRICS.out.metrics
                metrics_file_batchmetrics = GATKSV_GENERATEBATCHMETRICS.out.metrics_file_batchmetrics
                generate_batch_metrics_ploidy_table = GATKSV_GENERATEBATCHMETRICS.out.ploidy_table

                if (params.run_filter_batch_sites) {
                    filter_batch_sites_input = batch_cohort
                        .join(clustered_depth_vcf)
                        .join(clustered_manta_vcf)
                        .join(clustered_wham_vcf)
                        .join(clustered_scramble_vcf)
                        .join(generate_batch_metrics_metrics)
                        .map { batch_key, cohort, clustered_depth_vcf, clustered_manta_vcf, clustered_wham_vcf, clustered_scramble_vcf, evidence_metrics ->
                            tuple(
                                batch_key,
                                cohort,
                                selectSinglePath(clustered_depth_vcf, "clustered_depth_vcf"),
                                selectSinglePath(clustered_manta_vcf, "clustered_manta_vcf"),
                                selectSinglePath(clustered_wham_vcf, "clustered_wham_vcf"),
                                selectSinglePath(clustered_scramble_vcf, "clustered_scramble_vcf"),
                                selectSinglePath(evidence_metrics, "evidence_metrics")
                            )
                        }

                    GATKSV_FILTERBATCHSITES (
                        filter_batch_sites_input
                    )
                    versions = versions.mix(GATKSV_FILTERBATCHSITES.out.versions)
                    sites_filtered_depth_vcf = GATKSV_FILTERBATCHSITES.out.sites_filtered_depth_vcf
                    sites_filtered_manta_vcf = GATKSV_FILTERBATCHSITES.out.sites_filtered_manta_vcf
                    sites_filtered_wham_vcf = GATKSV_FILTERBATCHSITES.out.sites_filtered_wham_vcf
                    sites_filtered_scramble_vcf = GATKSV_FILTERBATCHSITES.out.sites_filtered_scramble_vcf
                    filter_batch_sites_cutoffs = GATKSV_FILTERBATCHSITES.out.cutoffs
                    filter_batch_sites_scores = GATKSV_FILTERBATCHSITES.out.scores
                    filter_batch_sites_RF_intermediate_files = GATKSV_FILTERBATCHSITES.out.RF_intermediate_files
                    filter_batch_sites_sv_counts = GATKSV_FILTERBATCHSITES.out.sv_counts
                    filter_batch_sites_sv_count_plots = GATKSV_FILTERBATCHSITES.out.sv_count_plots
                    filter_batch_sites_outlier_samples_preview = GATKSV_FILTERBATCHSITES.out.sites_filtered_outlier_samples_preview
                    filter_batch_sites_outlier_samples_with_reason = GATKSV_FILTERBATCHSITES.out.sites_filtered_outlier_samples_with_reason
                    filter_batch_sites_num_outlier_samples = GATKSV_FILTERBATCHSITES.out.num_outlier_samples
                    filter_batch_sites_num_outlier_samples_int = filter_batch_sites_num_outlier_samples
                        .map { batch_key, outliers_file ->
                            tuple(batch_key, outliers_file.text.trim().toInteger())
                        }

                    if (params.run_filter_batch_samples) {
                        filter_batch_samples_input = batch_cohort
                            .join(sites_filtered_depth_vcf)
                            .join(sites_filtered_manta_vcf)
                            .join(sites_filtered_wham_vcf)
                            .join(sites_filtered_scramble_vcf)
                            .join(filter_batch_sites_cutoffs)
                            .map { batch_key, cohort, sites_filtered_depth_vcf, sites_filtered_manta_vcf, sites_filtered_wham_vcf, sites_filtered_scramble_vcf, cutoffs ->
                                tuple(
                                    batch_key,
                                    cohort,
                                    selectSinglePath(sites_filtered_depth_vcf, "sites_filtered_depth_vcf"),
                                    selectSinglePath(sites_filtered_manta_vcf, "sites_filtered_manta_vcf"),
                                    selectSinglePath(sites_filtered_wham_vcf, "sites_filtered_wham_vcf"),
                                    selectSinglePath(sites_filtered_scramble_vcf, "sites_filtered_scramble_vcf"),
                                    selectSinglePath(cutoffs, "cutoffs")
                                )
                            }

                        GATKSV_FILTERBATCHSAMPLES (
                            filter_batch_samples_input
                        )
                        versions = versions.mix(GATKSV_FILTERBATCHSAMPLES.out.versions)
                        outlier_filtered_depth_vcf = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_depth_vcf
                        outlier_filtered_depth_vcf_index = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_depth_vcf_index
                        outlier_filtered_manta_vcf = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_manta_vcf
                        outlier_filtered_manta_vcf_index = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_manta_vcf_index
                        outlier_filtered_scramble_vcf = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_scramble_vcf
                        outlier_filtered_scramble_vcf_index = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_scramble_vcf_index
                        outlier_filtered_wham_vcf = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_wham_vcf
                        outlier_filtered_wham_vcf_index = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_wham_vcf_index
                        outlier_filtered_pesr_vcf = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_pesr_vcf
                        outlier_filtered_pesr_vcf_index = GATKSV_FILTERBATCHSAMPLES.out.outlier_filtered_pesr_vcf_index
                        filtered_depth_vcf = outlier_filtered_depth_vcf
                            .map { batch_key, file_value ->
                                tuple(batch_key, selectSinglePath(file_value, "filtered_depth_vcf"))
                            }
                        filtered_pesr_vcf = outlier_filtered_pesr_vcf
                            .map { batch_key, file_value ->
                                tuple(batch_key, selectSinglePath(file_value, "filtered_pesr_vcf"))
                            }
                        outlier_samples_excluded_file = GATKSV_FILTERBATCHSAMPLES.out.outlier_samples_excluded_file
                            .map { batch_key, file_value ->
                                tuple(batch_key, selectSinglePath(file_value, "outlier_samples_excluded_file"))
                            }
                        filtered_batch_samples_file = GATKSV_FILTERBATCHSAMPLES.out.filtered_batch_samples_file
                            .map { batch_key, file_value ->
                                tuple(batch_key, selectSinglePath(file_value, "filtered_batch_samples_file"))
                            }
                    }
                }
            }
        }
    }

    emit:
    versions
    generate_batch_metrics_metrics
    metrics_file_batchmetrics
    generate_batch_metrics_ploidy_table
    batch_cohort
    merged_BAF
    merged_BAF_index
    merged_PE
    merged_PE_index
    merged_bincov
    merged_bincov_index
    merged_SR
    merged_SR_index
    merged_dels
    merged_dups
    cnmops_del
    cnmops_del_index
    cnmops_dup
    cnmops_dup_index
    cnmops_large_del
    cnmops_large_del_index
    cnmops_large_dup
    cnmops_large_dup_index
    median_cov
    std_manta_vcf_tar
    std_wham_vcf_tar
    std_scramble_vcf_tar
    batch_ploidy_matrix
    batch_ploidy_plots
    BAF_stats
    RD_stats
    PE_stats
    SR_stats
    Matrix_QC_plot
    metrics_file_batchevidence
    manta_tloc
    filter_batch_sites_cutoffs
    filter_batch_sites_scores
    filter_batch_sites_RF_intermediate_files
    filter_batch_sites_sv_counts
    filter_batch_sites_sv_count_plots
    sites_filtered_depth_vcf
    sites_filtered_manta_vcf
    sites_filtered_wham_vcf
    sites_filtered_scramble_vcf
    filter_batch_sites_outlier_samples_preview
    filter_batch_sites_outlier_samples_with_reason
    filter_batch_sites_num_outlier_samples
    outlier_filtered_depth_vcf
    outlier_filtered_depth_vcf_index
    outlier_filtered_manta_vcf
    outlier_filtered_manta_vcf_index
    outlier_filtered_scramble_vcf
    outlier_filtered_scramble_vcf_index
    outlier_filtered_wham_vcf
    outlier_filtered_wham_vcf_index
    outlier_filtered_pesr_vcf
    outlier_filtered_pesr_vcf_index
    clustered_depth_vcf
    clustered_depth_vcf_index
    clustered_manta_vcf
    clustered_manta_vcf_index
    clustered_wham_vcf
    clustered_wham_vcf_index
    clustered_scramble_vcf
    clustered_scramble_vcf_index
    clustered_sv_counts
    clustered_sv_count_plots
    clustered_outlier_samples_preview
    clustered_outlier_samples_with_reason
    cluster_batch_num_outlier_samples
    metrics_file_clusterbatch
    cluster_batch_num_outlier_samples_int
    filter_batch_sites_num_outlier_samples_int
    filtered_depth_vcf
    filtered_pesr_vcf
    outlier_samples_excluded_file
    filtered_batch_samples_file
}
