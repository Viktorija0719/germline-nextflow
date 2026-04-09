/*
 * ExomeDepth cohort-level CNV calling.
 * Requires params.idxstats_enable = true (asserted in germline.nf).
 *
 * Input:
 *   ch_bam_bai     – [meta, bam, bai]  (all samples)
 *   ch_idxstats    – [meta, idxstats]  (from SAMTOOLS_IDXSTATS)
 *   ch_exd_bed_raw – path — raw BED (will be prepared inside)
 *   ch_ref_fasta   – val path
 *   ch_ref_fai     – val path
 *
 * ExomeDepth is cohort-wide — all samples are assembled into one tuple
 * and passed to the R script together.
 */

include { EXOMEDEPTH  } from '../../modules/local/exomedepth'
include { PREPARE_BED } from '../../modules/local/prepare_bed'

workflow CNV_EXOMEDEPTH {
    take:
    ch_bam_bai      // [meta, bam, bai]
    ch_idxstats     // [meta, idxstats_file]
    ch_exd_bed_raw  // path — raw BED
    ch_ref_fasta    // val path
    ch_ref_fai      // val path

    main:
    PREPARE_BED(ch_exd_bed_raw, ch_ref_fai)

    def meta_cohort = [ id: 'cohort', sample: 'cohort', patient: 'cohort' ]

    // Build sample_id → idxstats path map
    def ch_idx_map = ch_idxstats
        .map { meta, idx ->
            def sid = (meta instanceof Map && meta.id != null) ? meta.id.toString() : meta.toString()
            tuple(sid, idx)
        }
        .collect(flat: false)
        .map { pairs ->
            def m = [:]
            pairs.each { p -> m[p[0].toString()] = p[1] }
            m
        }

    // Collect cohort rows, combine with idxmap, bed, and fasta, then assemble the EXOMEDEPTH tuple
    ch_bam_bai
        .collect(flat: false)
        .map  { cohort -> tuple([cohort: cohort]) }   // box to prevent combine flattening
        .combine(ch_idx_map)
        .combine(PREPARE_BED.out.out_bed)
        .combine(ch_ref_fasta)
        .map { box, idxmap, bed, fasta ->
            def rows = box.cohort
                .collect { row ->
                    def m   = row[0]
                    def sid = (m instanceof Map && m.id     != null) ? m.id.toString()
                            : (m instanceof Map && m.sample != null) ? m.sample.toString()
                            : m.toString()
                    [ sid, row[1], row[2] ]
                }
                .sort { a, b -> a[0] <=> b[0] }

            def samples = rows.collect { it[0] }
            def bams    = rows.collect { it[1] }
            def bais    = rows.collect { it[2] }
            def idxs    = samples.collect { sid ->
                def p = idxmap[sid]
                if (!p) error "Missing idxstats for '${sid}' (available: ${idxmap.keySet().sort().join(', ')})"
                p
            }

            tuple(meta_cohort, samples, bams, bais, idxs, fasta, bed)
        }
        | EXOMEDEPTH
}
