/*
 * XHMM cohort CNV discovery from GATK DepthOfCoverage outputs.
 * Requires CNV_DEPTHOFCOVERAGE to have run first (asserted in germline.nf).
 *
 * Input:
 *   ch_doc_coverage – [meta, [coverage_files]] from CNV_DEPTHOFCOVERAGE
 *   ch_xhmm_bed_raw – path — raw XHMM intervals BED (will be prepared)
 *   ch_ref_fasta    – val path
 *   ch_ref_fai      – val path
 *   ch_ref_dict     – val path
 *   ch_xhmm_params  – val path
 *   ch_extreme_gc   – val path (empty file is acceptable)
 *
 * Emit:
 *   xcnv – [meta, xcnv]
 */

include { GATK_ANNOTATEINTERVALS } from '../../modules/local/gatk/annotateintervals'
include { XHMM_PIPELINE          } from '../../modules/local/xhmm'
include { PREPARE_BED            } from '../../modules/local/prepare_bed'

workflow CNV_XHMM {
    take:
    ch_doc_coverage // [meta, [coverage_files]]
    ch_xhmm_bed_raw // path — raw XHMM intervals BED
    ch_ref_fasta    // val path
    ch_ref_fai      // val path
    ch_ref_dict     // val path
    ch_xhmm_params  // val path
    ch_extreme_gc   // val path

    main:
    PREPARE_BED(ch_xhmm_bed_raw, ch_ref_fai)

    def meta_xhmm = [ id: (params.xhmm_prefix ?: 'DATA') ]

    // Extract .sample_interval_summary from per-sample DepthOfCoverage output
    def ch_xhmm_summaries = ch_doc_coverage
        .map { meta, files ->
            def sis = files.find { it.name.endsWith('.sample_interval_summary') }
            if (!sis) error "No .sample_interval_summary in DepthOfCoverage output for '${meta.id}'"
            sis
        }
        .collect()
        .map { summaries -> tuple(meta_xhmm, summaries) }

    GATK_ANNOTATEINTERVALS(
        PREPARE_BED.out.out_bed.map { bed -> tuple(meta_xhmm, bed) },
        ch_ref_fasta,
        ch_ref_fai,
        ch_ref_dict
    )

    XHMM_PIPELINE(
        ch_xhmm_summaries,
        GATK_ANNOTATEINTERVALS.out.tsv.map { meta, tsv -> tsv },
        ch_extreme_gc,
        ch_xhmm_params
    )

    emit:
    xcnv = XHMM_PIPELINE.out.xcnv
}
