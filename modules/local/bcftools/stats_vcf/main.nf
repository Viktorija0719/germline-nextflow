process BCFTOOLS_STATS_VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bcftools>=1.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data' :
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.vcf_stats.txt"), emit: stats
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-f PASS'
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    bcftools stats ${args} ${vcf} > ${prefix}.vcf_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf_stats.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
