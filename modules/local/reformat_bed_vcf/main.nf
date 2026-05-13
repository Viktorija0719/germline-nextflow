process REFORMAT_BED_VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python>=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--1' :
        'docker.io/python:3.11-slim' }"

    input:
    tuple val(meta), path(bed)
    path script

    output:
    tuple val(meta), path("${prefix}.vcf"), emit: vcf
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script} -i ${bed} -o ${prefix}.vcf -v

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
