/*
 * BWA alignment (lane-level) → read-group tagging → per-sample merge
 * → duplicate marking → BAM indexing.
 *
 * Input:
 *   ch_reads  – [meta(lane), [fq1, fq2]]
 *   ch_index  – val [meta_ref, bwa_dir]
 *   ch_ref    – val [meta_ref, fasta]
 *
 * Emit:
 *   bam_bai  – [meta(sample), bam, bai]
 *   bam      – [meta(sample), bam]          (post-dedup)
 *   metrics  – [meta(sample), metrics.txt]  (biobambam dup metrics)
 */

include { BWA_MEM                       } from '../../modules/nf-core/bwa/mem'
include { BIOBAMBAM_BAMMARKDUPLICATES2  } from '../../modules/nf-core/biobambam/bammarkduplicates2'
include { SAMTOOLS_INDEX                } from '../../modules/nf-core/samtools/index'
include { ADD_READGROUPS                } from '../../modules/local/picard/addreadgroups'
include { MERGE_BAMS                    } from '../../modules/local/samtools/merge'

workflow ALIGN_AND_DEDUP {
    take:
    ch_reads    // [meta, [fq1, fq2]]
    ch_index    // val [meta_ref, bwa_dir]
    ch_ref      // val [meta_ref, fasta]

    main:
    BWA_MEM(ch_reads, ch_index, ch_ref, true)
    ADD_READGROUPS(BWA_MEM.out.bam)

    // Group lane BAMs by sample, then merge
    ADD_READGROUPS.out.bam
        .map { meta, bam -> tuple(meta.sample, meta, bam) }
        .groupTuple()
        .map { sample, metas, bams ->
            def merged_meta = [ id: sample, sample: sample, patient: metas[0].patient ]
            tuple(merged_meta, bams)
        }
        | MERGE_BAMS

    BIOBAMBAM_BAMMARKDUPLICATES2(MERGE_BAMS.out.bam)
    SAMTOOLS_INDEX(BIOBAMBAM_BAMMARKDUPLICATES2.out.bam)

    // Join deduped BAM with its index into a 3-tuple
    def ch_bam_bai = BIOBAMBAM_BAMMARKDUPLICATES2.out.bam
        .map { meta, bam -> tuple(meta.id, meta, bam) }
        .join(SAMTOOLS_INDEX.out.bai.map { meta, bai -> tuple(meta.id, bai) })
        .map { id, meta, bam, bai -> tuple(meta, bam, bai) }

    emit:
    bam_bai = ch_bam_bai
    bam     = BIOBAMBAM_BAMMARKDUPLICATES2.out.bam
    metrics = BIOBAMBAM_BAMMARKDUPLICATES2.out.metrics
}
