process XCNV2VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python>=3.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://python:3.11-slim' :
    'docker.io/library/python:3.11-slim' }"

    input:
    tuple val(meta), path(xcnv)
    path py_script

    output:
    path "*.xcnv.vcf", emit: vcfs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${py_script} ${xcnv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
