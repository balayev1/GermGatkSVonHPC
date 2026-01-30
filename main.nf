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
    // This will save files to /scratch.global/balay011/GermlineSV_outs/CollectInsertSizeMetrics
    picard_results = COLLECT_INSERT_SIZE(samples_ch)

    // 3. Run GatherSampleEvidence
    gse_results = GATHER_SAMPLE_EVIDENCE(samples_ch)

    // 4. Run EvidenceQC (Jointly for all samples)
    // .collect() waits for all Step 3 tasks to finish
    EVIDENCE_QC(gse_results.full_evidence_dir.map{ it[1] }.collect())
}