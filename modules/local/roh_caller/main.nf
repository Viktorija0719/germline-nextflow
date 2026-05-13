process ROH_CALLER {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools>=1.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data' :
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.roh.txt"),      emit: roh
    tuple val(meta), path("${prefix}.roh.filt.txt"), emit: roh_filt
    tuple val(meta), path("${prefix}.roh.bed"),      emit: roh_bed
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--AF-dflt 0.4'
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    # Pre-filter: remove low-depth / low-quality variants before ROH detection
    bcftools view -e 'DP<10 || QUAL<7' ${vcf} | \\
        bcftools roh \\
            ${args} \\
            -o ${prefix}.roh.txt \\
            -

    # Extract only ROH region lines (RG prefix)
    grep '^RG' ${prefix}.roh.txt > ${prefix}.roh.filt.txt || touch ${prefix}.roh.filt.txt

    # Convert RG lines to BED: columns are RG, Sample, Chr, Start, End, Length, nSNPs, Quality
    awk 'BEGIN{OFS="\\t"} /^RG/{print \$3, \$4-1, \$5, \$2"_ROH_"NR, \$8, "."}' \\
        ${prefix}.roh.filt.txt > ${prefix}.roh.bed || touch ${prefix}.roh.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.roh.txt ${prefix}.roh.filt.txt ${prefix}.roh.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
