process ADD_READGROUPS {
    tag "${meta.id}"
    label 'process_small'
    container 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.rg.bam"), emit: bam

    script:
    def rgid     = meta.id
    def sample   = meta.sample   ?: meta.id
    def lib      = meta.library  ?: meta.sample ?: meta.id
    def platform = meta.platform ?: 'ILLUMINA'
    def pu       = meta.platform_unit ?: "${meta.sample}.${meta.lane ?: 'L1'}"
    """
    picard AddOrReplaceReadGroups \\
        I=${bam} \\
        O=${meta.id}.rg.bam \\
        RGID=${rgid} \\
        RGLB=${lib} \\
        RGPL=${platform} \\
        RGPU=${pu} \\
        RGSM=${sample}
    """
}
