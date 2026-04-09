/*
 * VCF merging and normalisation (all optional).
 *
 *   combine_dv_strelka_enable  → concat DV + Strelka2, optionally normalise
 *   norm_dv_only_enable        → normalise DeepVariant VCF independently
 *
 * Input:
 *   ch_dv_vcf      – [meta, vcf.gz, tbi]  (empty when DV disabled)
 *   ch_strelka_vcf – [meta, vcf.gz, tbi]  (empty when Strelka disabled)
 *   ch_ref         – val [meta_ref, fasta]
 */

include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF      } from '../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_NORM   as BCFTOOLS_NORM_COMBINED   } from '../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_NORM   as BCFTOOLS_NORM_DVONLY     } from '../../modules/nf-core/bcftools/norm'

workflow VARIANT_MERGE {
    take:
    ch_dv_vcf      // [meta, vcf.gz, tbi]  or empty channel
    ch_strelka_vcf // [meta, vcf.gz, tbi]  or empty channel
    ch_ref         // val [meta_ref, fasta]

    main:
    // ------------------------------------------------------------------
    // DV + Strelka2 concat (and optional normalisation)
    // ------------------------------------------------------------------
    if (params.combine_dv_strelka_enable) {
        if (!(params.deepvariant_enable && params.strelka_enable)) {
            if (params.strict_variant_requirements) {
                error "combine_dv_strelka_enable=true requires deepvariant_enable=true and strelka_enable=true"
            } else {
                log.warn "Skipping DV+Strelka2 combine: both deepvariant_enable and strelka_enable must be true"
            }
        } else {
            def ch_concat_in = ch_dv_vcf
                .map { meta, vcf, tbi -> tuple(meta.id, meta, vcf, tbi) }
                .join(ch_strelka_vcf.map { meta, vcf, tbi -> tuple(meta.id, meta, vcf, tbi) })
                .map { sid, meta_dv, dv_vcf, dv_tbi, meta_st, st_vcf, st_tbi ->
                    def meta    = (meta_st ?: meta_dv).clone()
                    meta.id     = "${sid}.DV_ST2.comb"
                    tuple(meta, [dv_vcf, st_vcf], [dv_tbi, st_tbi])
                }

            BCFTOOLS_CONCAT_VCF(ch_concat_in)

            if (params.norm_combined_enable) {
                def ch_norm_in = BCFTOOLS_CONCAT_VCF.out.vcf
                    .join(BCFTOOLS_CONCAT_VCF.out.tbi)
                    .map { meta, vcf, tbi ->
                        def m = meta.clone()
                        m.id  = "${meta.id}.norm"
                        tuple(m, vcf, tbi)
                    }
                BCFTOOLS_NORM_COMBINED(ch_norm_in, ch_ref)
            }
        }
    }

    // ------------------------------------------------------------------
    // DeepVariant-only normalisation
    // ------------------------------------------------------------------
    if (params.deepvariant_enable && params.norm_dv_only_enable) {
        def ch_dv_norm_in = ch_dv_vcf
            .map { meta, vcf, tbi ->
                def m = meta.clone()
                m.id  = "${meta.id}.DV_only.norm"
                tuple(m, vcf, tbi)
            }
        BCFTOOLS_NORM_DVONLY(ch_dv_norm_in, ch_ref)
    }
}
