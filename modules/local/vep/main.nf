process VEP_ANNOTATE {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::ensembl-vep>=111"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:111.0--pl5321h2a838b0_1' :
        'biocontainers/ensembl-vep:111.0--pl5321h2a838b0_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path cache_dir          // pass [] to use container-bundled cache
    path spliceai_snv       // SpliceAI SNV VCF; pass [] to skip plugin
    path spliceai_snv_tbi
    path spliceai_indel     // SpliceAI indel VCF; pass [] to skip plugin
    path spliceai_indel_tbi
    path mastermind         // Mastermind cited-variants VCF; pass [] to skip plugin
    path mastermind_tbi

    output:
    tuple val(meta), path("${prefix}.vep.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}.vep.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def fields     = task.ext.fields ?: 'SpliceAI_cutoff,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL,SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL,Mastermind_MMID3,Mastermind_counts,Consequence,IMPACT,SYMBOL,Gene,Feature,BIOTYPE,EXON,INTRON,cDNA_position,Protein_position,Amino_acids,Codons,STRAND,AF'
    def cache_opt  = cache_dir    ? "--cache --dir_cache ${cache_dir} --offline" : '--cache --offline'
    def sai_opt    = (spliceai_snv && spliceai_indel) ? "--plugin SpliceAI,snv=${spliceai_snv},indel=${spliceai_indel},cutoff=0.4" : ''
    def mm_opt     = mastermind   ? "--plugin Mastermind,${mastermind}" : ''
    prefix         = task.ext.prefix ?: "${meta.id}"
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${prefix}.vep.vcf.gz \\
        ${cache_opt} \\
        --format vcf \\
        --vcf \\
        --compress_output gzip \\
        --fork ${task.cpus} \\
        --gencode_basic \\
        --no_stats \\
        --force_overwrite \\
        --allow_non_variant \\
        --dont_skip \\
        --distance 100 \\
        --af \\
        --fields "${fields}" \\
        ${sai_opt} \\
        ${mm_opt} \\
        ${args}

    tabix -p vcf ${prefix}.vep.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(vep --help 2>&1 | grep -i 'version' | head -1 | sed 's/.*version //' || echo unknown)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.vep.vcf.gz
    touch ${prefix}.vep.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: stub
    END_VERSIONS
    """
}
