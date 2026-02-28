process PED_CHECK {
    tag "$cohort"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(cohort), path(ped_file), path(sample_list)

    output:
    tuple val(cohort), path("ped_validation.pass"), emit: pass
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    validate_ped.py \\
        -p ${ped_file} \\
        -s ${sample_list}

    touch ped_validation.pass

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
