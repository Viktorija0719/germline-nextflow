/*
 * Reformat biobambam metrics to look like Picard DuplicationMetrics
 * so MultiQC can parse them without a container (avoids BeeGFS SIGBUS
 * issues with coreutils in Singularity on some HPC setups).
 */
process PICARDLIKE_DUPMETRICS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(metrics)

    output:
    tuple val(meta), path("${meta.id}.duplicate_metrics_picard_like.txt"), emit: metrics

    script:
    """
    {
        echo '##htsjdk.samtools.metrics.StringHeader'
        echo '##METRICS CLASS picard.sam.DuplicationMetrics'
        awk 'NR>3' ${metrics}
    } > ${meta.id}.duplicate_metrics_picard_like.txt
    """
}
