/*
 * modules/local/pad_bed/main.nf
 *
 * Produce a +100 bp padded, clipped-to-chromosome-bounds, sorted, merged BED
 * from a canonical (unpadded) target BED.
 *
 * Output is stored permanently in resources/ so the process is skipped on
 * subsequent runs if the file already exists there.
 */

process PAD_BED {
    tag "${bed.baseName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"

    // Store permanently in resources/ — skipped automatically if file exists
    storeDir "${projectDir}/resources"

    input:
    path bed
    path fasta_fai

    output:
    path "${bed.baseName}.plus100.merged.bed", emit: padded_bed

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    LC_ALL=C sort -k1,1 -k2,2n "${bed}" \
    | bedtools slop -i - -g "${fasta_fai}" -b 100 \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${bed.baseName}.plus100.merged.bed"
    """
}
