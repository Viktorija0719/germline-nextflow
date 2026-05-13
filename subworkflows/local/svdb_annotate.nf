/*
 * SVDB cohort-frequency annotation for Manta SV, XHMM CNV, and ExomeDepth CNV.
 *
 * Per-tool arm (always runs when tool is enabled):
 *   1. Convert non-VCF → VCF  (XCNV2VCF, EXOMEDEPTH_CNV2VCF)
 *   2. Build cohort SVDB database from all samples
 *   3. Query each sample against cohort DB
 *   4. Soft-tag common variants [params.frq_enable]   ← soft filter: keeps ALL variants
 *   5. AnnotSV + knotAnnotSV per tool [params.annotsv_enable]
 *
 * Merged arm (params.svdb_merge_enable — requires ≥2 tools with data):
 *   6. SVDB_MERGE: join ED + XHMM + Manta VCFs per sample
 *   7. BEDTOOLS_MERGE_CNVS: split by type, merge overlapping, recombine
 *   8. REFORMAT_BED_VCF: BED → VCF
 *   9. Build second cohort DB from merged VCFs + query
 *  10. Soft-tag merged at FRQ ≥ 0.05 (TooCommon) — all variants preserved
 *  11. AnnotSV + knotAnnotSV on merged [params.annotsv_enable]
 */

include { XCNV2VCF                            } from '../../modules/local/xcnv2vcf'
include { EXOMEDEPTH_CNV2VCF                  } from '../../modules/local/exomedepth_cnv2vcf'
include { SVDB_BUILD as SVDB_BUILD_MANTA      } from '../../modules/nf-core/svdb/build'
include { SVDB_BUILD as SVDB_BUILD_XHMM       } from '../../modules/nf-core/svdb/build'
include { SVDB_BUILD as SVDB_BUILD_ED         } from '../../modules/nf-core/svdb/build'
include { SVDB_BUILD as SVDB_BUILD_MERGED     } from '../../modules/nf-core/svdb/build'
include { SVDB_QUERY as SVDB_QUERY_MANTA      } from '../../modules/nf-core/svdb/query'
include { SVDB_QUERY as SVDB_QUERY_XHMM       } from '../../modules/nf-core/svdb/query'
include { SVDB_QUERY as SVDB_QUERY_ED         } from '../../modules/nf-core/svdb/query'
include { SVDB_QUERY as SVDB_QUERY_MERGED     } from '../../modules/nf-core/svdb/query'
include { SVDB_MERGE as SVDB_MERGE_CALLERS    } from '../../modules/nf-core/svdb/merge'
include { FILTER_COMMON_CNVS as FILTER_COMMON_CNVS_ED    } from '../../modules/local/filter_common_cnvs'
include { FILTER_COMMON_CNVS as FILTER_COMMON_CNVS_XHMM  } from '../../modules/local/filter_common_cnvs'
include { FILTER_COMMON_CNVS as FILTER_COMMON_CNVS_MANTA } from '../../modules/local/filter_common_cnvs'
include { BCFTOOLS_FILTER as SOFT_FILTER_MERGED           } from '../../modules/nf-core/bcftools/filter'
include { BEDTOOLS_MERGE_CNVS                 } from '../../modules/local/bedtools/merge_cnvs'
include { REFORMAT_BED_VCF                    } from '../../modules/local/reformat_bed_vcf'
include { ANNOTSV as ANNOTSV_ED               } from '../../modules/local/annotsv'
include { ANNOTSV as ANNOTSV_XHMM             } from '../../modules/local/annotsv'
include { ANNOTSV as ANNOTSV_MANTA            } from '../../modules/local/annotsv'
include { ANNOTSV as ANNOTSV_MERGED           } from '../../modules/local/annotsv'
include { KNOTANNOTSV as KNOTANNOTSV_ED       } from '../../modules/local/knotannotsv'
include { KNOTANNOTSV as KNOTANNOTSV_XHMM     } from '../../modules/local/knotannotsv'
include { KNOTANNOTSV as KNOTANNOTSV_MANTA    } from '../../modules/local/knotannotsv'
include { KNOTANNOTSV as KNOTANNOTSV_MERGED   } from '../../modules/local/knotannotsv'

workflow SVDB_ANNOTATE {
    take:
    ch_manta_vcf    // [meta, vcf.gz]  — Channel.empty() when no Manta data
    ch_xhmm_xcnv    // [meta, xcnv]    — Channel.empty() when no XHMM data
    ch_exd_calls    // [meta, csv]     — Channel.empty() when no ExomeDepth data
    ch_ref_fai_path // val path        — reference FAI (for bedtools sort in merge arm)

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
        .first()

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
    // Soft frequency filter — params.frq_enable
    //
    // SOFT filter: adds FILTER=TooCommon (and TooLarge for Manta) tags.
    // ALL variants are preserved — downstream tools see FILTER flags.
    // When disabled: SVDB_QUERY outputs pass through unchanged.
    // ------------------------------------------------------------------
    def ch_final_ed    = Channel.empty()
    def ch_final_xhmm  = Channel.empty()
    def ch_final_manta = Channel.empty()

    if (params.frq_enable) {
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
    // Per-tool AnnotSV + knotAnnotSV — params.annotsv_enable
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
                             ? Channel.value(file(params.annotsv_config_file, checkIfExists: true))
                             : Channel.value([])
        def ch_panels      = params.annotsv_panels_tsv
                             ? Channel.value(file(params.annotsv_panels_tsv, checkIfExists: true))
                             : Channel.value([])

        ANNOTSV_ED(ch_final_ed,    ch_annotations, ch_config, ch_panels)
        ANNOTSV_XHMM(ch_final_xhmm,  ch_annotations, ch_config, ch_panels)
        ANNOTSV_MANTA(ch_final_manta, ch_annotations, ch_config, ch_panels)

        ch_annotsv_ed_tsv    = ANNOTSV_ED.out.tsv
        ch_annotsv_xhmm_tsv  = ANNOTSV_XHMM.out.tsv
        ch_annotsv_manta_tsv = ANNOTSV_MANTA.out.tsv
    }

    if (params.knotannotsv_enable && params.annotsv_enable) {
        def ch_knot_cfg = params.knotannotsv_config_file
                          ? Channel.value(file(params.knotannotsv_config_file, checkIfExists: true))
                          : Channel.value([])

        KNOTANNOTSV_ED(ch_annotsv_ed_tsv,    ch_knot_cfg)
        KNOTANNOTSV_XHMM(ch_annotsv_xhmm_tsv,  ch_knot_cfg)
        KNOTANNOTSV_MANTA(ch_annotsv_manta_tsv, ch_knot_cfg)
    }

    // ------------------------------------------------------------------
    // Merged arm — params.svdb_merge_enable
    //
    // Joins ED + XHMM + Manta VCFs per sample (needs ≥2 tools).
    // Runs bedtools merge per SVTYPE, reformats BED → VCF,
    // annotates with a second cohort-level SVDB, then soft-flags common.
    // ------------------------------------------------------------------
    if (params.svdb_merge_enable) {

        // Group per-sample VCFs from all enabled tools by sample ID.
        // Alphabetical sort ensures consistent priority: ed < manta < xhmm
        // (matches pipeline_nextflow.txt: --priority ed,xhmm,manta)
        def ch_for_merge = ch_final_manta
            .map { meta, vcf -> tuple(meta.id, 'manta', vcf) }
            .mix(ch_final_xhmm.map  { meta, vcf -> tuple(meta.id, 'xhmm',  vcf) })
            .mix(ch_final_ed.map    { meta, vcf -> tuple(meta.id, 'ed',    vcf) })
            .groupTuple(by: 0)
            .filter { sid, callers, vcfs -> vcfs.size() >= 2 }
            .map    { sid, callers, vcfs ->
                def pairs          = [callers, vcfs].transpose().sort { a, b -> a[0] <=> b[0] }
                def sorted_callers = pairs.collect { it[0] }
                def sorted_vcfs    = pairs.collect { it[1] }
                tuple([id: sid], sorted_vcfs, sorted_callers)
            }

        ch_for_merge.multiMap { meta, vcfs, callers ->
            meta_vcfs: tuple(meta, vcfs)
            priority:  callers
        }.set { ch_svdb_merge_in }

        SVDB_MERGE_CALLERS(ch_svdb_merge_in.meta_vcfs, ch_svdb_merge_in.priority, false)

        BEDTOOLS_MERGE_CNVS(SVDB_MERGE_CALLERS.out.vcf, ch_ref_fai_path)

        REFORMAT_BED_VCF(
            BEDTOOLS_MERGE_CNVS.out.bed,
            Channel.value(file("${projectDir}/modules/local/reformat_bed_vcf/reformat_bed_vcf.py"))
        )

        // Second SVDB: build cohort DB from merged VCFs, query each sample
        REFORMAT_BED_VCF.out.vcf.multiMap { meta, vcf ->
            for_build: vcf
            for_query: tuple(meta, vcf)
        }.set { ch_merged }

        def ch_merged_build = ch_merged.for_build
            .collect()
            .filter { it.size() > 0 }
            .map    { vcfs -> tuple([id: 'merged_cohort'], vcfs) }

        SVDB_BUILD_MERGED(ch_merged_build, 'files')

        def ch_merged_db = SVDB_BUILD_MERGED.out.db
            .map   { meta, db -> [db] }
            .first()

        SVDB_QUERY_MERGED(ch_merged.for_query, [], [], [], [], ch_merged_db, [])

        // Soft-tag: FRQ ≥ 0.05 → FILTER=TooCommon; all variants kept
        SOFT_FILTER_MERGED(
            SVDB_QUERY_MERGED.out.vcf.map { meta, vcf -> tuple(meta, vcf, []) }
        )

        // AnnotSV on the merged VCF (if enabled)
        if (params.annotsv_enable) {
            def ch_annotations_m = Channel.value(file(params.annotsv_annotations_dir, checkIfExists: true))
            def ch_config_m      = params.annotsv_config_file
                                   ? Channel.value(file(params.annotsv_config_file, checkIfExists: true))
                                   : Channel.value([])
            def ch_panels_m      = params.annotsv_panels_tsv
                                   ? Channel.value(file(params.annotsv_panels_tsv, checkIfExists: true))
                                   : Channel.value([])

            ANNOTSV_MERGED(SOFT_FILTER_MERGED.out.vcf, ch_annotations_m, ch_config_m, ch_panels_m)

            if (params.knotannotsv_enable) {
                def ch_knot_cfg_m = params.knotannotsv_config_file
                                    ? Channel.value(file(params.knotannotsv_config_file, checkIfExists: true))
                                    : Channel.value([])
                KNOTANNOTSV_MERGED(ANNOTSV_MERGED.out.tsv, ch_knot_cfg_m)
            }
        }
    }

    emit:
    manta_filtered   = ch_final_manta   // [meta, vcf] — soft-filtered or raw SVDB query
    xhmm_filtered    = ch_final_xhmm
    ed_filtered      = ch_final_ed
}
