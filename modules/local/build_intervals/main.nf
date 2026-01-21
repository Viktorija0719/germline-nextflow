process BUILD_INTERVALS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'quay.io/biocontainers/gawk:5.3.1' }"
    

    input:
    tuple val(meta), path(fasta_fai)

    output:
    tuple val(meta), path("${fasta_fai.baseName}.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Window size for Case B. If 0 or unset => one interval per contig.
    def window_bp = (params.genome_window_bp ?: 0) as long

    // Optional regex to drop contigs (decoys etc.). Example: '^(GL|KI|chrUn)'
    def exclude_regex = params.build_intervals_exclude_regex ?: ''

    """
    set -euo pipefail

    awk -v FS='\\t' -v OFS='\\t' \\
        -v W=${window_bp} \\
        -v EXCL='${exclude_regex}' \\
        '{
            contig=\$1; len=\$2
            if (EXCL != "" && contig ~ EXCL) next

            if (W <= 0) {
                print contig, 0, len
            } else {
                for (i=0; i<len; i+=W) {
                    j=i+W; if (j>len) j=len
                    print contig, i, j
                }
            }
        }' ${fasta_fai} > ${fasta_fai.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
