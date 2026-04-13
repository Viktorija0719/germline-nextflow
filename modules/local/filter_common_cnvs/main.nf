process FILTER_COMMON_CNVS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bcftools>=1.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data' :
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta),
          path(ed_vcf,   stageAs: 'ed_input.vcf'),
          path(xhmm_vcf, stageAs: 'xhmm_input.vcf'),
          path(manta_vcf, stageAs: 'manta_input.vcf')

    output:
    tuple val(meta), path("${meta.id}.frq_filt_ed.vcf"),    emit: ed_vcf
    tuple val(meta), path("${meta.id}.frq_filt_xhmm.vcf"),  emit: xhmm_vcf
    tuple val(meta), path("${meta.id}.frq_filt_manta.vcf"), emit: manta_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def frq_threshold = task.ext.frq_threshold ?: '0.1'
    def svlen_max     = task.ext.svlen_max      ?: '1000000'
    """
    bcftools filter \\
        -e 'INFO/FRQ>${frq_threshold}' \\
        -o ${meta.id}.frq_filt_ed.vcf \\
        ${ed_vcf}

    bcftools filter \\
        -e 'INFO/FRQ>${frq_threshold}' \\
        -o ${meta.id}.frq_filt_xhmm.vcf \\
        ${xhmm_vcf}

    bcftools view -f PASS ${manta_vcf} | \\
        bcftools filter \\
            -e 'INFO/FRQ>${frq_threshold} || INFO/SVLEN>${svlen_max} || INFO/SVLEN<-${svlen_max}' \\
            -o ${meta.id}.frq_filt_manta.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
