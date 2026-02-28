process MAKE_PED_SAMPLE_LIST {
    tag "$cohort"

    input:
    tuple val(cohort), val(sample_ids)

    output:
    tuple val(cohort), path("${cohort}.samples.list"), emit: sample_list

    script:
    def ids = (sample_ids instanceof List ? sample_ids : [sample_ids])
        .collect { it.toString() }
        .unique()
        .sort()
    def quoted_ids = ids.collect { it.inspect() }.join(' ')
    """
    printf "%s\\n" ${quoted_ids} > ${cohort}.samples.list
    """
}
