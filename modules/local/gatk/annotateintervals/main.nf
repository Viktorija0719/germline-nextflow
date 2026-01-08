process GATK_ANNOTATEINTERVALS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://broadinstitute/gatk:4.5.0.0' :
        'docker.io/broadinstitute/gatk:4.5.0.0' }"

    input:
    tuple val(meta), path(intervals)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple val(meta), path("*.annotated_intervals.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def extra_args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: (meta.id ?: 'DATA')

    """
    set -euo pipefail

    ln -sf ${ref_fasta} ref.fa
    ln -sf ${ref_fai}   ref.fa.fai
    ln -sf ${ref_dict}  ref.dict

    gatk AnnotateIntervals \\
      -R ref.fa \\
      -L ${intervals} \\
      -O ${prefix}.annotated_intervals.tsv \\
      ${extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: "\$(gatk --version 2>&1 | head -n 1 || echo unknown)"
    END_VERSIONS
    """
}
