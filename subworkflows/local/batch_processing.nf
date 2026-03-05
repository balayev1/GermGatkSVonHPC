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
    train_gcnv_input  // channel: [ cohort_or_batch, sample_ids, count_files ]
    batch_evidence_input

    main:
    versions = Channel.empty()
    def gcnv_model_tars = Channel.empty()
    def contig_ploidy_model_tar = Channel.empty()

    if (params.run_train_gcnv) {
        GATK_TRAINGCNV (
            train_gcnv_input
        )
        versions = versions.mix(GATK_TRAINGCNV.out.versions)
        gcnv_model_tars = GATK_TRAINGCNV.out.cohort_gcnv_model_tars
        contig_ploidy_model_tar = GATK_TRAINGCNV.out.cohort_contig_ploidy_model_tar
    }

    def gather_batch_evidence_results = Channel.empty()
    def cluster_batch_results = Channel.empty()
    def generate_batch_metrics_results = Channel.empty()
    def generate_batch_metrics_metrics = Channel.empty()
    def generate_batch_metrics_ploidy_table = Channel.empty()
    def batch_cohort = Channel.empty()
    def merged_PE = Channel.empty()
    def merged_bincov = Channel.empty()
    def merged_SR = Channel.empty()
    def median_cov = Channel.empty()
    def filter_batch_sites_results = Channel.empty()
    def filter_batch_sites_cutoffs = Channel.empty()
    def filter_batch_sites_sv_counts = Channel.empty()
    def filter_batch_sites_sv_count_plots = Channel.empty()
    def filter_batch_samples_results = Channel.empty()
    def filtered_depth_vcf = Channel.empty()
    def filtered_pesr_vcf = Channel.empty()
    def outlier_samples_excluded_file = Channel.empty()
    def filtered_batch_samples_file = Channel.empty()
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
            if (candidates.size() == 1) {
                return candidates[0]
            }
            def merged = candidates.findAll { it.getName().toLowerCase().contains("merged") }
            return merged ? merged[0] : candidates[0]
        }

        batch_cohort = batch_evidence_input
            .map { batch_key, cohort, ped_file, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files ->
                tuple(batch_key, cohort)
            }
            .unique()

        batch_evidence_with_models = batch_evidence_input
            .join(gcnv_model_tars)
            .map { batch_key, cohort, ped_file, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files, gcnv_model_tars ->
                tuple(
                    batch_key,
                    cohort,
                    ped_file,
                    sample_ids,
                    counts_files,
                    pe_files,
                    sr_files,
                    sd_files,
                    manta_files,
                    wham_files,
                    scramble_files,
                    gcnv_model_tars
                )
            }
            .join(contig_ploidy_model_tar)
            .map { batch_key, cohort, ped_file, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files, gcnv_model_tars, contig_ploidy_model_tar ->
                tuple(
                    batch_key,
                    cohort,
                    ped_file,
                    asList(sample_ids).collect { it.toString() },
                    asList(counts_files),
                    asList(pe_files),
                    asList(sr_files),
                    asList(sd_files),
                    asList(manta_files),
                    asList(wham_files),
                    asList(scramble_files),
                    asList(gcnv_model_tars),
                    contig_ploidy_model_tar
                )
            }

        GATKSV_GATHERBATCHEVIDENCE (
            batch_evidence_with_models
        )
        versions = versions.mix(GATKSV_GATHERBATCHEVIDENCE.out.versions)
        gather_batch_evidence_results = GATKSV_GATHERBATCHEVIDENCE.out.gather_batch_evidence_results
        merged_PE = GATKSV_GATHERBATCHEVIDENCE.out.merged_PE
        merged_bincov = GATKSV_GATHERBATCHEVIDENCE.out.merged_bincov
        merged_SR = GATKSV_GATHERBATCHEVIDENCE.out.merged_SR
        median_cov = GATKSV_GATHERBATCHEVIDENCE.out.median_cov

        if (params.run_cluster_batch) {
            ped_by_batch = batch_evidence_input
                .map { batch_key, cohort, ped_file, sample_ids, counts_files, pe_files, sr_files, sd_files, manta_files, wham_files, scramble_files ->
                    tuple(batch_key, cohort, ped_file)
                }

            cluster_batch_input = ped_by_batch
                .join(GATKSV_GATHERBATCHEVIDENCE.out.merged_dels)
                .join(GATKSV_GATHERBATCHEVIDENCE.out.merged_dups)
                .join(GATKSV_GATHERBATCHEVIDENCE.out.std_wham_vcf_tar)
                .join(GATKSV_GATHERBATCHEVIDENCE.out.std_manta_vcf_tar)
                .join(GATKSV_GATHERBATCHEVIDENCE.out.std_scramble_vcf_tar)
                .map { batch_key, cohort, ped_file, merged_dels, merged_dups, std_wham_vcf_tar, std_manta_vcf_tar, std_scramble_vcf_tar ->
                    tuple(
                        batch_key,
                        cohort,
                        ped_file,
                        selectSinglePath(merged_dels, "merged_dels"),
                        selectSinglePath(merged_dups, "merged_dups"),
                        selectSinglePath(std_wham_vcf_tar, "std_wham_vcf_tar"),
                        selectSinglePath(std_manta_vcf_tar, "std_manta_vcf_tar"),
                        selectSinglePath(std_scramble_vcf_tar, "std_scramble_vcf_tar")
                    )
                }

            GATKSV_CLUSTERBATCH (
                cluster_batch_input
            )
            versions = versions.mix(GATKSV_CLUSTERBATCH.out.versions)
            cluster_batch_results = GATKSV_CLUSTERBATCH.out.cluster_batch_results

            if (params.run_generate_batch_metrics) {
                generate_batch_metrics_input = ped_by_batch
                    .join(GATKSV_GATHERBATCHEVIDENCE.out.merged_PE)
                    .join(GATKSV_GATHERBATCHEVIDENCE.out.merged_BAF)
                    .join(GATKSV_GATHERBATCHEVIDENCE.out.merged_bincov)
                    .join(GATKSV_GATHERBATCHEVIDENCE.out.merged_SR)
                    .join(GATKSV_GATHERBATCHEVIDENCE.out.median_cov)
                    .join(GATKSV_CLUSTERBATCH.out.clustered_depth_vcf)
                    .join(GATKSV_CLUSTERBATCH.out.clustered_manta_vcf)
                    .join(GATKSV_CLUSTERBATCH.out.clustered_wham_vcf)
                    .join(GATKSV_CLUSTERBATCH.out.clustered_scramble_vcf)
                    .map { batch_key, cohort, ped_file, merged_pe, merged_baf, merged_bincov, merged_sr, median_cov, clustered_depth_vcf, clustered_manta_vcf, clustered_wham_vcf, clustered_scramble_vcf ->
                        tuple(
                            batch_key,
                            cohort,
                            ped_file,
                            selectSinglePath(merged_pe, "merged_PE"),
                            selectSinglePath(merged_baf, "merged_BAF"),
                            selectSinglePath(merged_bincov, "merged_bincov"),
                            selectSinglePath(merged_sr, "merged_SR"),
                            selectSinglePath(median_cov, "median_cov"),
                            selectSinglePath(clustered_depth_vcf, "clustered_depth_vcf"),
                            selectSinglePath(clustered_manta_vcf, "clustered_manta_vcf"),
                            selectSinglePath(clustered_wham_vcf, "clustered_wham_vcf"),
                            selectSinglePath(clustered_scramble_vcf, "clustered_scramble_vcf")
                        )
                    }

                GATKSV_GENERATEBATCHMETRICS (
                    generate_batch_metrics_input
                )
                versions = versions.mix(GATKSV_GENERATEBATCHMETRICS.out.versions)
                generate_batch_metrics_results = GATKSV_GENERATEBATCHMETRICS.out.generate_batch_metrics_results
                generate_batch_metrics_metrics = GATKSV_GENERATEBATCHMETRICS.out.metrics
                generate_batch_metrics_ploidy_table = GATKSV_GENERATEBATCHMETRICS.out.ploidy_table

                if (params.run_filter_batch_sites) {
                    filter_batch_sites_input = batch_cohort
                        .join(GATKSV_CLUSTERBATCH.out.clustered_depth_vcf)
                        .join(GATKSV_CLUSTERBATCH.out.clustered_manta_vcf)
                        .join(GATKSV_CLUSTERBATCH.out.clustered_wham_vcf)
                        .join(GATKSV_CLUSTERBATCH.out.clustered_scramble_vcf)
                        .join(GATKSV_GENERATEBATCHMETRICS.out.metrics)
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
                    filter_batch_sites_results = GATKSV_FILTERBATCHSITES.out.filter_batch_sites_results
                    filter_batch_sites_cutoffs = GATKSV_FILTERBATCHSITES.out.cutoffs
                    filter_batch_sites_sv_counts = GATKSV_FILTERBATCHSITES.out.sv_counts
                    filter_batch_sites_sv_count_plots = GATKSV_FILTERBATCHSITES.out.sv_count_plots

                    if (params.run_filter_batch_samples) {
                        filter_batch_samples_input = batch_cohort
                            .join(GATKSV_FILTERBATCHSITES.out.sites_filtered_depth_vcf)
                            .join(GATKSV_FILTERBATCHSITES.out.sites_filtered_manta_vcf)
                            .join(GATKSV_FILTERBATCHSITES.out.sites_filtered_wham_vcf)
                            .join(GATKSV_FILTERBATCHSITES.out.sites_filtered_scramble_vcf)
                            .join(GATKSV_FILTERBATCHSITES.out.cutoffs)
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
                        filter_batch_samples_results = GATKSV_FILTERBATCHSAMPLES.out.filter_batch_samples_results
                        filtered_depth_vcf = GATKSV_FILTERBATCHSAMPLES.out.filtered_depth_vcf
                        filtered_pesr_vcf = GATKSV_FILTERBATCHSAMPLES.out.filtered_pesr_vcf
                        outlier_samples_excluded_file = GATKSV_FILTERBATCHSAMPLES.out.outlier_samples_excluded_file
                        filtered_batch_samples_file = GATKSV_FILTERBATCHSAMPLES.out.filtered_batch_samples_file
                    }
                }
            }
        }
    }

    emit:
    versions
    gather_batch_evidence_results
    cluster_batch_results
    generate_batch_metrics_results
    generate_batch_metrics_metrics
    generate_batch_metrics_ploidy_table
    batch_cohort
    merged_PE
    merged_bincov
    merged_SR
    median_cov
    filter_batch_sites_results
    filter_batch_sites_cutoffs
    filter_batch_sites_sv_counts
    filter_batch_sites_sv_count_plots
    filter_batch_samples_results
    filtered_depth_vcf
    filtered_pesr_vcf
    outlier_samples_excluded_file
    filtered_batch_samples_file
}
