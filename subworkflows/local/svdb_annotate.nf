/*
 * SVDB cohort-frequency annotation for Manta SV, XHMM CNV, and ExomeDepth CNV.
 *
 * For each tool the pipeline:
 *   1. Converts non-VCF outputs to VCF  (XCNV2VCF, EXOMEDEPTH_CNV2VCF)
 *   2. Collects all cohort VCFs → builds a cohort SVDB SQLite database
 *   3. Queries each sample against that cohort database
 *   4. (Optional) Filters common variants by cohort frequency  [params.frq_enable]
 *   5. (Optional) Annotates each per-sample VCF with AnnotSV   [params.annotsv_enable]
 *
 * No external databases are required. Annotation arms are driven purely by
 * whether the input channels carry data — an empty channel is a no-op.
 * Use multiMap before collect() to avoid splitting a queue channel between
 * two consumers (build and query).
 */

include { XCNV2VCF                       } from '../../modules/local/xcnv2vcf'
include { EXOMEDEPTH_CNV2VCF             } from '../../modules/local/exomedepth_cnv2vcf'
include { SVDB_BUILD as SVDB_BUILD_MANTA } from '../../modules/nf-core/svdb/build'
include { SVDB_BUILD as SVDB_BUILD_XHMM  } from '../../modules/nf-core/svdb/build'
include { SVDB_BUILD as SVDB_BUILD_ED    } from '../../modules/nf-core/svdb/build'
include { SVDB_QUERY as SVDB_QUERY_MANTA } from '../../modules/nf-core/svdb/query'
include { SVDB_QUERY as SVDB_QUERY_XHMM  } from '../../modules/nf-core/svdb/query'
include { SVDB_QUERY as SVDB_QUERY_ED    } from '../../modules/nf-core/svdb/query'
include { FILTER_COMMON_CNVS as FILTER_COMMON_CNVS_ED    } from '../../modules/local/filter_common_cnvs'
include { FILTER_COMMON_CNVS as FILTER_COMMON_CNVS_XHMM  } from '../../modules/local/filter_common_cnvs'
include { FILTER_COMMON_CNVS as FILTER_COMMON_CNVS_MANTA } from '../../modules/local/filter_common_cnvs'
include { ANNOTSV as ANNOTSV_ED               } from '../../modules/local/annotsv'
include { ANNOTSV as ANNOTSV_XHMM             } from '../../modules/local/annotsv'
include { ANNOTSV as ANNOTSV_MANTA            } from '../../modules/local/annotsv'
include { KNOTANNOTSV as KNOTANNOTSV_ED       } from '../../modules/local/knotannotsv'
include { KNOTANNOTSV as KNOTANNOTSV_XHMM     } from '../../modules/local/knotannotsv'
include { KNOTANNOTSV as KNOTANNOTSV_MANTA    } from '../../modules/local/knotannotsv'

workflow SVDB_ANNOTATE {
    take:
    ch_manta_vcf   // [meta, vcf.gz]  — Channel.empty() when no Manta data
    ch_xhmm_xcnv   // [meta, xcnv]    — Channel.empty() when no XHMM data
    ch_exd_calls   // [meta, csv]     — Channel.empty() when no ExomeDepth data

    main:

    // ------------------------------------------------------------------
    // Manta: cohort diploid_sv.vcf.gz → SVDB build → per-sample query
    // ------------------------------------------------------------------
    ch_manta_vcf.multiMap { meta, vcf ->
        for_build: vcf
        for_query: tuple(meta, vcf)
    }.set { ch_manta }

    def ch_manta_build = ch_manta.for_build
        .collect()
        .filter  { it.size() > 0 }
        .map     { vcfs -> tuple([id: 'manta_cohort'], vcfs) }

    SVDB_BUILD_MANTA(ch_manta_build, 'files')

    def ch_manta_db = SVDB_BUILD_MANTA.out.db
        .map   { meta, db -> [db] }
        .first()                      // value channel — broadcasts to all N queries

    SVDB_QUERY_MANTA(ch_manta.for_query, [], [], [], [], ch_manta_db, [])

    // ------------------------------------------------------------------
    // XHMM: DATA.xcnv → per-sample VCF → SVDB build → per-sample query
    // ------------------------------------------------------------------
    XCNV2VCF(
        ch_xhmm_xcnv,
        Channel.value(file("${projectDir}/modules/local/xcnv2vcf/xcnv2vcf.py"))
    )

    def ch_xhmm_vcfs = XCNV2VCF.out.vcfs
        .flatten()
        .map { vcf ->
            def sample = vcf.name.replaceAll(/\.xcnv\.vcf$/, '')
            tuple([id: sample, sample: sample], vcf)
        }

    ch_xhmm_vcfs.multiMap { meta, vcf ->
        for_build: vcf
        for_query: tuple(meta, vcf)
    }.set { ch_xhmm }

    def ch_xhmm_build = ch_xhmm.for_build
        .collect()
        .filter  { it.size() > 0 }
        .map     { vcfs -> tuple([id: 'xhmm_cohort'], vcfs) }

    SVDB_BUILD_XHMM(ch_xhmm_build, 'files')

    def ch_xhmm_db = SVDB_BUILD_XHMM.out.db
        .map   { meta, db -> [db] }
        .first()

    SVDB_QUERY_XHMM(ch_xhmm.for_query, [], [], [], [], ch_xhmm_db, [])

    // ------------------------------------------------------------------
    // ExomeDepth: per-sample CSV → VCF → SVDB build → per-sample query
    // ------------------------------------------------------------------
    EXOMEDEPTH_CNV2VCF(
        ch_exd_calls,
        Channel.value(file("${projectDir}/modules/local/exomedepth_cnv2vcf/exomedepth_cnv2vcf.py"))
    )

    EXOMEDEPTH_CNV2VCF.out.vcf.multiMap { meta, vcf ->
        for_build: vcf
        for_query: tuple(meta, vcf)
    }.set { ch_ed }

    def ch_ed_build = ch_ed.for_build
        .collect()
        .filter  { it.size() > 0 }
        .map     { vcfs -> tuple([id: 'exomedepth_cohort'], vcfs) }

    SVDB_BUILD_ED(ch_ed_build, 'files')

    def ch_ed_db = SVDB_BUILD_ED.out.db
        .map   { meta, db -> [db] }
        .first()

    SVDB_QUERY_ED(ch_ed.for_query, [], [], [], [], ch_ed_db, [])

    // ------------------------------------------------------------------
    // (Optional) Frequency filter — params.frq_enable
    //
    // When enabled:  FILTER_COMMON_CNVS runs and its outputs feed the
    //                next step (AnnotSV).
    // When disabled: SVDB_QUERY outputs are used directly.
    //
    // The filter requires a 4-tuple joining all three tool outputs.
    // If any tool was disabled its SVDB_QUERY arm is empty, so the inner
    // join produces no items and the filter is effectively a no-op.
    // ------------------------------------------------------------------
    def ch_final_ed    = Channel.empty()
    def ch_final_xhmm  = Channel.empty()
    def ch_final_manta = Channel.empty()

    if (params.frq_enable) {
        // Each tool is filtered independently — a missing tool (e.g. XHMM when
        // DepthOfCoverage fails) simply produces an empty channel for that arm
        // without blocking annotation of the other tools.
        FILTER_COMMON_CNVS_ED(SVDB_QUERY_ED.out.vcf)
        FILTER_COMMON_CNVS_XHMM(SVDB_QUERY_XHMM.out.vcf)
        FILTER_COMMON_CNVS_MANTA(SVDB_QUERY_MANTA.out.vcf)

        ch_final_ed    = FILTER_COMMON_CNVS_ED.out.vcf
        ch_final_xhmm  = FILTER_COMMON_CNVS_XHMM.out.vcf
        ch_final_manta = FILTER_COMMON_CNVS_MANTA.out.vcf
    } else {
        ch_final_ed    = SVDB_QUERY_ED.out.vcf
        ch_final_xhmm  = SVDB_QUERY_XHMM.out.vcf
        ch_final_manta = SVDB_QUERY_MANTA.out.vcf
    }

    // ------------------------------------------------------------------
    // (Optional) AnnotSV annotation — params.annotsv_enable
    //
    // Runs on the final per-sample VCF from each tool arm independently.
    // Pass [] for config_file or panels_tsv when the param is not set.
    // ------------------------------------------------------------------
    def ch_annotsv_ed_tsv    = Channel.empty()
    def ch_annotsv_xhmm_tsv  = Channel.empty()
    def ch_annotsv_manta_tsv = Channel.empty()

    if (params.annotsv_enable) {
        if (!params.annotsv_annotations_dir) {
            error "annotsv_enable=true requires annotsv_annotations_dir to be set"
        }
        def ch_annotations = Channel.value(file(params.annotsv_annotations_dir, checkIfExists: true))
        def ch_config      = params.annotsv_config_file
                             ? Channel.value(file(params.annotsv_config_file))
                             : Channel.value([])
        def ch_panels      = params.annotsv_panels_tsv
                             ? Channel.value(file(params.annotsv_panels_tsv))
                             : Channel.value([])

        ANNOTSV_ED(ch_final_ed,    ch_annotations, ch_config, ch_panels)
        ANNOTSV_XHMM(ch_final_xhmm,  ch_annotations, ch_config, ch_panels)
        ANNOTSV_MANTA(ch_final_manta, ch_annotations, ch_config, ch_panels)

        ch_annotsv_ed_tsv    = ANNOTSV_ED.out.tsv
        ch_annotsv_xhmm_tsv  = ANNOTSV_XHMM.out.tsv
        ch_annotsv_manta_tsv = ANNOTSV_MANTA.out.tsv
    }

    // ------------------------------------------------------------------
    // (Optional) knotAnnotSV — params.knotannotsv_enable
    //
    // Converts each AnnotSV TSV to an annotated Excel (.xlsm) workbook.
    // Runs only when annotsv_enable AND knotannotsv_enable are both true.
    // ------------------------------------------------------------------
    def ch_knot_ed_xl    = Channel.empty()
    def ch_knot_xhmm_xl  = Channel.empty()
    def ch_knot_manta_xl = Channel.empty()

    if (params.knotannotsv_enable && params.annotsv_enable) {
        def ch_knot_cfg = params.knotannotsv_config_file
                          ? Channel.value(file(params.knotannotsv_config_file, checkIfExists: true))
                          : Channel.value([])

        KNOTANNOTSV_ED(ch_annotsv_ed_tsv,    ch_knot_cfg)
        KNOTANNOTSV_XHMM(ch_annotsv_xhmm_tsv,  ch_knot_cfg)
        KNOTANNOTSV_MANTA(ch_annotsv_manta_tsv, ch_knot_cfg)

        ch_knot_ed_xl    = KNOTANNOTSV_ED.out.xl
        ch_knot_xhmm_xl  = KNOTANNOTSV_XHMM.out.xl
        ch_knot_manta_xl = KNOTANNOTSV_MANTA.out.xl
    }

    emit:
    manta_annotated  = SVDB_QUERY_MANTA.out.vcf     // [meta, vcf] — raw SVDB-annotated
    xhmm_annotated   = SVDB_QUERY_XHMM.out.vcf      // [meta, vcf]
    ed_annotated     = SVDB_QUERY_ED.out.vcf         // [meta, vcf]
    manta_filtered   = ch_final_manta                // [meta, vcf] — FRQ-filtered or same as above
    xhmm_filtered    = ch_final_xhmm                 // [meta, vcf]
    ed_filtered      = ch_final_ed                   // [meta, vcf]
    manta_annotsv    = ch_annotsv_manta_tsv           // [meta, tsv] — AnnotSV output
    xhmm_annotsv     = ch_annotsv_xhmm_tsv            // [meta, tsv]
    ed_annotsv       = ch_annotsv_ed_tsv              // [meta, tsv]
    manta_knot       = ch_knot_manta_xl               // [meta, xlsm] — knotAnnotSV output
    xhmm_knot        = ch_knot_xhmm_xl                // [meta, xlsm]
    ed_knot          = ch_knot_ed_xl                  // [meta, xlsm]
}
