process KNOTANNOTSV {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/knotannotsv:1.1.5--hdfd78af_0'
        : 'biocontainers/knotannotsv:1.1.5--hdfd78af_0'}"

    input:
    tuple val(meta), path(annotsv_tsv)
    path(config_file)            // custom config (pass [] to use container default)

    output:
    tuple val(meta), path("*.xlsm"), emit: xl,   optional: true
    tuple val(meta), path("*.html"), emit: html,  optional: true
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Use user-supplied config if staged, otherwise fall back to container default
    def cfg    = config_file
                 ? "--configFile ${config_file}"
                 : "--configFile \${CONDA_PREFIX:-/usr/local}/share/knotAnnotSV/config_AnnotSV.yaml"
    """
    # knotAnnotSV names the output after the input filename; rename to control prefix
    mv ${annotsv_tsv} ${prefix}.tsv

    knotAnnotSV2XL.pl \\
        --annotSVfile ${prefix}.tsv \\
        --genomeBuild hg38 \\
        ${cfg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        knotAnnotSV: \$(knotAnnotSV2XL.pl --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.xlsm
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        knotAnnotSV: stub
    END_VERSIONS
    """
}
