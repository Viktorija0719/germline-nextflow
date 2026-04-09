/*
 * Alignment QC — all tools are individually toggled via params.
 *
 * Input:
 *   ch_bam_bai     – [meta, bam, bai]   (for idxstats, verifybamid2)
 *   ch_bam         – [meta, bam]        (for fastqc, qualimap)
 *   ch_qualimap_bed – path              (prepared feature BED)
 *   ch_ref         – val [meta_ref, fasta]
 *   ch_vbid_svd    – val [ud, mu, bed]
 *   ch_vbid_refvcf – val file
 *
 * Emit:
 *   idxstats – [meta, idxstats_file]   (empty channel when disabled)
 */

include { FASTQC                   } from '../../modules/nf-core/fastqc'
include { SAMTOOLS_IDXSTATS        } from '../../modules/nf-core/samtools/idxstats'
include { VERIFYBAMID_VERIFYBAMID2 } from '../../modules/nf-core/verifybamid/verifybamid2'
include { QUALIMAP_BAMQC           } from '../../modules/nf-core/qualimap/bamqc'

workflow BAM_QC {
    take:
    ch_bam_bai      // [meta, bam, bai]
    ch_bam          // [meta, bam]
    ch_qualimap_bed // path (single prepared BED)
    ch_ref          // val [meta_ref, fasta]
    ch_vbid_svd     // val [ud, mu, bed]
    ch_vbid_refvcf  // val file

    main:
    ch_idxstats = Channel.empty()

    if (params.fastqc_enable) {
        FASTQC(ch_bam)
    }

    if (params.idxstats_enable) {
        SAMTOOLS_IDXSTATS(ch_bam_bai)
        ch_idxstats = SAMTOOLS_IDXSTATS.out.idxstats
    }

    if (params.verifybamid2_enable) {
        VERIFYBAMID_VERIFYBAMID2(
            ch_bam_bai,
            ch_vbid_svd,
            ch_vbid_refvcf,
            ch_ref.map { meta, fasta -> fasta }
        )
    }

    if (params.qualimap_enable) {
        QUALIMAP_BAMQC(ch_bam, ch_qualimap_bed)
    }

    emit:
    idxstats = ch_idxstats  // [meta, idxstats_file] — empty when idxstats_enable=false
}
