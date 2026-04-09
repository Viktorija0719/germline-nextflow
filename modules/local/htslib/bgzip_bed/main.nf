process BGZIP_BED {
    tag "${bed.baseName}"
    label 'process_single'
    container 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'

    input:
    path bed

    output:
    path "${bed.baseName}.bed.gz", emit: bedgz

    script:
    """
    bgzip -c ${bed} > ${bed.baseName}.bed.gz
    """
}
