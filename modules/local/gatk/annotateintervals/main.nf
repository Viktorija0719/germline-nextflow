process GATK_ANNOTATEINTERVALS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://broadinstitute/gatk:4.5.0.0' :
        'docker.io/broadinstitute/gatk:4.5.0.0' }"
        
        input:
        tuple val(meta), path(intervals_bed)
        path fasta
        path fasta_fai
        path fasta_dict

        output:
        tuple val(meta), path("${meta.id}.annotated_intervals.tsv"), emit: tsv
        path "versions.yml", emit: versions

        when:
        task.ext.when == null || task.ext.when

        script:
        def args = task.ext.args ?: ''
        def out  = "${meta.id}.annotated_intervals.tsv"

        """
        set -euo pipefail

        # IMPORTANT for docker/non-root + HPC: tmp must exist & be writable
        mkdir -p tmp
        chmod 1777 tmp || true

        ln -sf ${fasta}      ref.fa
        ln -sf ${fasta_fai}  ref.fa.fai
        ln -sf ${fasta_dict} ref.dict

        gatk --java-options "-Djava.io.tmpdir=tmp" AnnotateIntervals \\
        -R ref.fa \\
        -L ${intervals_bed} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --tmp-dir tmp \\
        -O ${out} \\
        ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: "\$(gatk --version 2>&1 | head -n 1 || echo unknown)"
        END_VERSIONS
        """
    }