process GATHER_DEEPVARIANT_VCFS {
    tag "${meta.id}"
    label 'process_small'
    container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_2'

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("${meta.id}.deepvariant.vcf.gz"), path("${meta.id}.deepvariant.vcf.gz.tbi"), emit: vcf

    script:
    """
    printf "%s\\n" ${vcfs} | LC_ALL=C sort -V > vcfs.list
    bcftools concat -a -D -Oz -o ${meta.id}.deepvariant.vcf.gz -f vcfs.list
    bcftools index -t ${meta.id}.deepvariant.vcf.gz
    """
}
