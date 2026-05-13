/*
 * SNV annotation: VEP → soft-tag SpliceAI35 / MM → bcftools stats.
 * Activated when params.vep_enable = true.
 * All variants are preserved throughout; FILTER flags are informational only.
 *
 * Input:  ch_vcf – [meta, vcf.gz, tbi]  (normalized DV+Strelka or DV-only VCF)
 */

include { VEP_ANNOTATE    } from '../../modules/local/vep'
include { SNV_SOFT_FILTER } from '../../modules/local/snv_soft_filter'
include { BCFTOOLS_STATS_VCF } from '../../modules/local/bcftools/stats_vcf'

workflow SNV_ANNOTATE {
    take:
    ch_vcf   // [meta, vcf.gz, tbi]

    main:
    def ch_final_vcf = Channel.empty()

    if (params.vep_enable) {
        if (!params.vep_cache_dir) {
            error "vep_enable=true requires vep_cache_dir to be set"
        }

        def ch_cache          = Channel.value(file(params.vep_cache_dir, checkIfExists: true))
        def ch_spliceai_snv   = params.spliceai_snv_vcf   ? Channel.value(file(params.spliceai_snv_vcf))                       : Channel.value([])
        def ch_spliceai_sntbi = params.spliceai_snv_vcf   ? Channel.value(file("${params.spliceai_snv_vcf}.tbi"))               : Channel.value([])
        def ch_spliceai_indel = params.spliceai_indel_vcf ? Channel.value(file(params.spliceai_indel_vcf))                     : Channel.value([])
        def ch_spliceai_itbi  = params.spliceai_indel_vcf ? Channel.value(file("${params.spliceai_indel_vcf}.tbi"))             : Channel.value([])
        def ch_mastermind     = params.mastermind_vcf     ? Channel.value(file(params.mastermind_vcf))                         : Channel.value([])
        def ch_mastermind_tbi = params.mastermind_vcf     ? Channel.value(file("${params.mastermind_vcf}.tbi"))                 : Channel.value([])

        VEP_ANNOTATE(
            ch_vcf,
            ch_cache,
            ch_spliceai_snv,
            ch_spliceai_sntbi,
            ch_spliceai_indel,
            ch_spliceai_itbi,
            ch_mastermind,
            ch_mastermind_tbi
        )

        // Soft-tag SpliceAI35 and MM variants — all variants kept
        SNV_SOFT_FILTER(VEP_ANNOTATE.out.vcf.join(VEP_ANNOTATE.out.tbi))

        // QC stats on PASS variants
        BCFTOOLS_STATS_VCF(SNV_SOFT_FILTER.out.vcf.join(SNV_SOFT_FILTER.out.tbi))

        ch_final_vcf = SNV_SOFT_FILTER.out.vcf.join(SNV_SOFT_FILTER.out.tbi)
    }

    emit:
    vcf = ch_final_vcf   // [meta, vcf.gz, tbi]  (empty when vep_enable=false)
}
