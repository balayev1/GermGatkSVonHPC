#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import processes (defined below)

include { GATHER_SAMPLE_EVIDENCE } from './modules/gather_sample_evidence.nf'
include { EVIDENCE_QC     } from './modules/evidence_qc.nf'
include { COLLECT_INSERT_SIZE } from './modules/collect_insert_size.nf'

workflow {
    // 1. Create channel from manifest
    // Structure: [sample_id, bam, bai]
    samples_ch = Channel
        .fromPath(params.manifest)
        .splitCsv(sep: '\t')
        .map { row -> 
            def meta = [ id: row[0], gender: row[3] ]
            
            def bam = file(row[1])
            def bai = file(row[2])

            return [ meta, bam, bai ] }

    // 2. Run Picard CollectInsertSizeMetrics
    picard_results = COLLECT_INSERT_SIZE(samples_ch)

    // 3. Run GatherSampleEvidence
    gse_results = GATHER_SAMPLE_EVIDENCE(samples_ch)

    // 4. Run EvidenceQC (Jointly for all samples)
    // Collect all sample IDs and folders
    sample_ids = gse_results.full_evidence_dir.map{ it[0].id }.collect()
    evidence_folders = gse_results.full_evidence_dir.map{ it[1] }.collect()

    // Execute Evidence QC
    evidence_qc_results = EVIDENCE_QC(sample_ids, evidence_folders)

    // 5. Run Sample QC (Jointly for all samples)
    // Get the manifest metadata file
    meta_file = gse_results.full_evidence_dir
        .map { meta, path -> "${meta.id}\t${meta.gender}" }
        .collectFile(name: 'sample_metadata.tsv', newLine: true)

    // Count the number of samples from the manifest channel
    num_samples_ch = samples_ch.count()

    // Execute Sample QC 
    SAMPLE_QC(meta_file, workflow.workDir, num_samples_ch, picard_results, evidence_qc_results)
}