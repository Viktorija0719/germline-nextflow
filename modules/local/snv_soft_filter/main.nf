// Soft-tags SpliceAI high-confidence and Mastermind-cited SNVs.
// ALL variants are KEPT — flags are purely informational for downstream prioritization.
process SNV_SOFT_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bcftools>=1.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data' :
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.soft_filt.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}.soft_filt.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # SpliceAI35: CSQ contains "PASS" → SpliceAI score ≥ 0.35 cutoff
    # MM:         CSQ contains "&"    → multiple Mastermind citation entries
    bcftools filter \\
        -e 'INFO/CSQ ~ "PASS"' \\
        -m + -s 'SpliceAI35' \\
        ${vcf} | \\
    bcftools filter \\
        -e 'INFO/CSQ ~ "&"' \\
        -m + -s 'MM' \\
        -Oz -o ${prefix}.soft_filt.vcf.gz

    tabix -p vcf ${prefix}.soft_filt.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.soft_filt.vcf.gz
    touch ${prefix}.soft_filt.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
