/*
 * SNV/indel/SV calling: DeepVariant (scattered), Strelka2, Manta.
 * Each caller is individually toggled via params.
 *
 * Input:
 *   ch_bam_bai        – [meta, bam, bai]
 *   ch_var_bed_raw    – path — raw variant target BED
 *   ch_master_bed_raw – path — master intervals BED (for DV scatter)
 *   ch_ref            – val [meta_ref, fasta]
 *   ch_ref_fai        – val [meta_ref, fai]
 *   ch_ref_fasta_path – val path
 *   ch_ref_fai_path   – val path
 *   ch_deep_gzi       – val [meta_ref, gzi_dummy]
 *   ch_deep_par       – val [meta_ref, par_regions_bed]
 *   ch_manta_config   – val path
 *
 * Emit:
 *   dv_vcf      – [meta, vcf.gz, tbi]  (empty when deepvariant_enable=false)
 *   strelka_vcf – [meta, vcf.gz, tbi]  (empty when strelka_enable=false)
 */

include { DEEPVARIANT_RUNDEEPVARIANT } from '../../modules/nf-core/deepvariant/rundeepvariant'
include { STRELKA_GERMLINE           } from '../../modules/nf-core/strelka/germline'
include { MANTA_GERMLINE             } from '../../modules/nf-core/manta/germline'
include { TABIX_TABIX                } from '../../modules/nf-core/tabix/tabix'
include { STRELKA_CHRM_FILTER        } from '../../modules/local/strelka/chrm_filter'
include { GATHER_DEEPVARIANT_VCFS    } from '../../modules/local/deepvariant/gather_vcfs'
include { BGZIP_BED                  } from '../../modules/local/htslib/bgzip_bed'
include { CREATE_INTERVALS_BED       } from '../../modules/local/create_intervals_bed'
include { PREPARE_BED as PREPARE_BED_MASTER  } from '../../modules/local/prepare_bed'
include { PREPARE_BED as PREPARE_BED_VARIANT } from '../../modules/local/prepare_bed'

workflow VARIANT_CALLING {
    take:
    ch_bam_bai        // [meta, bam, bai]
    ch_var_bed_raw    // path — variant target BED
    ch_master_bed_raw // path — master intervals BED (DeepVariant scatter)
    ch_ref            // val [meta_ref, fasta]
    ch_ref_fai        // val [meta_ref, fai]
    ch_ref_fasta_path // val path
    ch_ref_fai_path   // val path
    ch_deep_gzi       // val [meta_ref, gzi]
    ch_deep_par       // val [meta_ref, par_bed]
    ch_manta_config   // val path

    main:
    ch_dv_vcf      = Channel.empty()
    ch_strelka_vcf = Channel.empty()

    // ------------------------------------------------------------------
    // DeepVariant — scatter by interval chunks, then gather per sample
    // ------------------------------------------------------------------
    if (params.deepvariant_enable) {
        PREPARE_BED_MASTER(ch_master_bed_raw, ch_ref_fai_path)

        def ch_interval_chunks
        if (params.no_intervals) {
            ch_interval_chunks = PREPARE_BED_MASTER.out.out_bed
        } else {
            CREATE_INTERVALS_BED(PREPARE_BED_MASTER.out.out_bed, ch_ref_fai_path)
            ch_interval_chunks = CREATE_INTERVALS_BED.out.bed
        }

        def ch_dv_scatter = ch_bam_bai
            .combine(ch_interval_chunks)
            .map { meta, bam, bai, bed ->
                def m = meta.clone()
                m.sample_id = meta.id
                m.id        = "${meta.id}__${bed.baseName}"
                tuple(m, bam, bai, bed)
            }

        DEEPVARIANT_RUNDEEPVARIANT(ch_dv_scatter, ch_ref, ch_ref_fai, ch_deep_gzi, ch_deep_par)

        DEEPVARIANT_RUNDEEPVARIANT.out.vcf
            .join(DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index)
            .map { meta, vcf, tbi ->
                def sample_meta = [ id: meta.sample_id, sample: meta.sample_id, patient: meta.patient ]
                tuple(sample_meta, vcf, tbi)
            }
            .groupTuple()
            .map { meta, vcfs, tbis -> tuple(meta, vcfs, tbis) }
            | GATHER_DEEPVARIANT_VCFS

        ch_dv_vcf = GATHER_DEEPVARIANT_VCFS.out.vcf
    }

    // ------------------------------------------------------------------
    // Shared bgzipped + indexed BED for Strelka2 and Manta
    // ------------------------------------------------------------------
    def ch_bedgz  = Channel.empty()
    def ch_bedtbi = Channel.empty()

    if (params.strelka_enable || params.manta_enable) {
        PREPARE_BED_VARIANT(ch_var_bed_raw, ch_ref_fai_path)
        BGZIP_BED(PREPARE_BED_VARIANT.out.out_bed)
        ch_bedgz = BGZIP_BED.out.bedgz

        TABIX_TABIX(ch_bedgz.map { bedgz -> tuple([ id: 'variant_bed' ], bedgz) })
        ch_bedtbi = TABIX_TABIX.out.index.map { meta, idx -> idx }
    }

    // ------------------------------------------------------------------
    // Strelka2
    // ------------------------------------------------------------------
    if (params.strelka_enable) {
        def ch_strelka_in = ch_bam_bai
            .combine(ch_bedgz)
            .combine(ch_bedtbi)
            .map { meta, bam, bai, bedgz, bedtbi -> tuple(meta, bam, bai, bedgz, bedtbi) }

        STRELKA_GERMLINE(ch_strelka_in, ch_ref_fasta_path, ch_ref_fai_path)

        STRELKA_CHRM_FILTER(
            STRELKA_GERMLINE.out.vcf.join(STRELKA_GERMLINE.out.vcf_tbi)
        )

        ch_strelka_vcf = STRELKA_CHRM_FILTER.out.vcf.join(STRELKA_CHRM_FILTER.out.tbi)
    }

    // ------------------------------------------------------------------
    // Manta
    // ------------------------------------------------------------------
    if (params.manta_enable) {
        def ch_manta_in = ch_bam_bai
            .combine(ch_bedgz)
            .combine(ch_bedtbi)
            .map { meta, bam, bai, bedgz, bedtbi -> tuple(meta, [bam], bai, bedgz, bedtbi) }

        MANTA_GERMLINE(ch_manta_in, ch_ref, ch_ref_fai, ch_manta_config)
    }

    emit:
    dv_vcf      = ch_dv_vcf       // [meta, vcf.gz, tbi]  or empty
    strelka_vcf = ch_strelka_vcf  // [meta, vcf.gz, tbi]  or empty
}
