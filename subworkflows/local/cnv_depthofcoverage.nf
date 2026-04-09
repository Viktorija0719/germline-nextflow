/*
 * GATK DepthOfCoverage — per-sample coverage depth over target intervals.
 * Output is consumed by CNV_XHMM when xhmm_enable=true, and published
 * independently when depthofcoverage_enable=true.
 *
 * Input:
 *   ch_bam_bai     – [meta, bam, bai]
 *   ch_cov_bed_raw – path — raw coverage target BED (will be prepared)
 *   ch_ref_fasta   – val path
 *   ch_ref_fai     – val path
 *   ch_ref_dict    – val path
 *   ch_gene_list   – val path
 *
 * Emit:
 *   coverage – [meta, [coverage_files]]
 */

include { GATK_DEPTHOFCOVERAGE } from '../../modules/local/gatk/depthofcoverage'
include { PREPARE_BED          } from '../../modules/local/prepare_bed'

workflow CNV_DEPTHOFCOVERAGE {
    take:
    ch_bam_bai     // [meta, bam, bai]
    ch_cov_bed_raw // path — raw BED
    ch_ref_fasta   // val path
    ch_ref_fai     // val path
    ch_ref_dict    // val path
    ch_gene_list   // val path

    main:
    PREPARE_BED(ch_cov_bed_raw, ch_ref_fai)

    GATK_DEPTHOFCOVERAGE(
        ch_bam_bai,
        ch_ref_fasta,
        ch_ref_fai,
        ch_ref_dict,
        PREPARE_BED.out.out_bed,
        ch_gene_list
    )

    emit:
    coverage = GATK_DEPTHOFCOVERAGE.out.coverage  // [meta, [files]]
}
