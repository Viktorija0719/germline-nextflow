process STRELKA_CHRM_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    /*
     * INPUT:
     *   - meta       : sample metadata
     *   - vcf        : Strelka SNV VCF (variants.vcf.gz)
     *   - vcf_index  : corresponding index (variants.vcf.gz.tbi)
     */
    input:
    tuple val(meta), path(vcf), path(vcf_index)

    /*
     * OUTPUT:
     *   - filtered VCF (haploid GT="1" removed, only PASS)
     *   - its tabix index
     */

    output:
    tuple val(meta), path("*.filtered.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.filtered.vcf.gz.tbi") , emit: tbi
    path "versions.yml"                            , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.strelka2"

    """
    set -euo pipefail

    in_vcf="${vcf}"
    prefix="${prefix}"

    # 1) Remove haploid GT="1" calls (typically chrM continuous-VF calls)
    bcftools view -e 'GT="1"' -Oz -o \${prefix}.no_haploid.vcf.gz "\${in_vcf}"
    tabix -p vcf \${prefix}.no_haploid.vcf.gz

    # 2) Keep only FILTER=PASS variants
    bcftools view -f PASS -Oz -o \${prefix}.filtered.vcf.gz \${prefix}.no_haploid.vcf.gz
    tabix -p vcf \${prefix}.filtered.vcf.gz

    # 3) Clean intermediate
    rm \${prefix}.no_haploid.vcf.gz \${prefix}.no_haploid.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        tabix:    \$(tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.strelka2"
    """
    echo "" | gzip > ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
        tabix:    "stub"
    END_VERSIONS
    """
}
