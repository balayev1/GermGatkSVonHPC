//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { PED_CHECK } from '../../modules/local/ped_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    validated_rows = SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .set { rows_ch }

    ped_inputs = rows_ch
        .map { row -> tuple(row.sample_set_set_id, row.cohort_ped, row.sample_id) }
        .groupTuple()
        .map { cohort, ped_files, sample_ids ->
            def unique_peds = ped_files.unique()
            if (unique_peds.size() != 1) {
                exit 1, "ERROR: Please check input samplesheet -> cohort ${cohort} has multiple PED files: ${unique_peds.join(', ')}"
            }
            tuple(cohort, file(unique_peds[0]), sample_ids.unique().sort())
        }

    ped_pass = PED_CHECK(ped_inputs).pass

    rows_ch
        .map { row ->
            def bam_meta = create_bam_channel(row)
            tuple(row.sample_set_set_id, bam_meta)
        }
        .join(ped_pass)
        .map { cohort, bam_meta, _ -> bam_meta }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions.mix(PED_CHECK.out.versions) // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ bam_or_cram ],[ bai_or_crai] ]
def create_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample_id
    meta.batch      = row.sample_set_id
    meta.cohort     = row.sample_set_set_id
    meta.ped        = row.cohort_ped

    // add path(s) of the bam file(s) to the meta map
    def bam_meta = []
    if (!file(row.bam_or_cram).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> BAM/CRAM file does not exist!\n${row.bam_or_cram}"
        }
    if (!file(row.bai_or_crai).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> BAI/CRAI file does not exist!\n${row.bai_or_crai}"
        }
    bam_meta = [ meta, file(row.bam_or_cram), file(row.bai_or_crai) ]

    return bam_meta
}
