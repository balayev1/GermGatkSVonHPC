process PED_CHECK {
    tag "$cohort"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(cohort), path(ped_file), val(sample_list)

    output:
    tuple val(cohort), path("ped_validation.pass"), emit: pass
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample_lines = sample_list.join('\n')
    """
    cat > samples.list << 'EOF'
    ${sample_lines}
    EOF

    python validate_ped.py \\
        -p ${ped_file} \\
        -s samples.list

    touch ped_validation.pass

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
