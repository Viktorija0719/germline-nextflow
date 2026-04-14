/*
 * Top-level germline variant calling workflow.
 *
 * Orchestrates:
 *   1. Reference channel setup
 *   2. BAM strategy (reuse existing or align from FASTQ)
 *   3. QC
 *   4. CNV calling (ExomeDepth, XHMM)
 *   5. Variant calling (DeepVariant, Strelka2, Manta)
 *   6. VCF merging / normalisation
 */

include { BUILD_INTERVALS       } from '../modules/local/build_intervals'
include { PICARDLIKE_DUPMETRICS } from '../modules/local/utils/picardlike_dupmetrics'

include { ALIGN_AND_DEDUP             } from '../subworkflows/local/align_and_dedup'
include { BAM_QC                      } from '../subworkflows/local/bam_qc'
include { CNV_EXOMEDEPTH              } from '../subworkflows/local/cnv_exomedepth'
include { CNV_DEPTHOFCOVERAGE         } from '../subworkflows/local/cnv_depthofcoverage'
include { CNV_XHMM                    } from '../subworkflows/local/cnv_xhmm'
include { VARIANT_CALLING             } from '../subworkflows/local/variant_calling'
include { VARIANT_MERGE               } from '../subworkflows/local/variant_merge'
include { PREPARE_BED as PREPARE_BED_QUALIMAP } from '../modules/local/prepare_bed'
include { PAD_BED                            } from '../modules/local/pad_bed'
include { SVDB_ANNOTATE                      } from '../subworkflows/local/svdb_annotate'

workflow GERMLINE {

    // ------------------------------------------------------------------
    // Validate required params
    // ------------------------------------------------------------------
    if (!params.genome || params.igenomes_ignore) {
        error "Provide --genome (e.g. GATK.GRCh38) and do not set --igenomes_ignore."
    }

    // ------------------------------------------------------------------
    // Reference channels
    // ------------------------------------------------------------------
    def meta_ref    = [ id: params.genome ]
    def fasta_path  = file(Utils.getGenomeAttr(params, 'fasta'))
    def bwa_path    = file(Utils.getGenomeAttr(params, 'bwa'))
    def fai_path    = file(Utils.getGenomeAttr(params, 'fasta_fai'))
    def dict_path   = file(Utils.getGenomeAttr(params, 'dict'))

    def ch_ref      = Channel.value(tuple(meta_ref, fasta_path))
    def ch_ref_fai  = Channel.value(tuple(meta_ref, fai_path))
    def ch_index    = Channel.value(tuple(meta_ref, bwa_path))

    def ch_ref_fasta_path = Channel.value(fasta_path)
    def ch_ref_fai_path   = Channel.value(fai_path)
    def ch_ref_dict_path  = Channel.value(dict_path)

    // ------------------------------------------------------------------
    // Build genome-wide intervals from FAI (used when no BED provided)
    // ------------------------------------------------------------------
    BUILD_INTERVALS(ch_ref_fai)
    def ch_genome_bed = BUILD_INTERVALS.out.bed.map { meta, bed -> bed }

    // ------------------------------------------------------------------
    // BAM strategy:
    //   1. BAM samplesheet  → load paths from CSV directly (skip alignment)
    //   2. Existing BAM dir → reuse precomputed BAMs from outdir/bam
    //   3. FASTQ samplesheet → align + dedup
    // ------------------------------------------------------------------
    def ssType = Utils.detectSampleSheetType(params.input)

    def ch_bam_final
    def ch_bai
    def ch_dup_metrics = null
    boolean run_picardlike = false

    if (ssType == 'bam') {
        // ------------------------------------------------------------------
        // BAM samplesheet: paths are declared in the CSV — skip alignment
        // ------------------------------------------------------------------
        log.info "BAM samplesheet detected — skipping alignment and deduplication."

        def ch_ss = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def sample   = row.sample.toString()
                def patient  = (row.patient ?: sample).toString()
                def meta     = [ id: sample, sample: sample, patient: patient ]
                def bam      = file(row.bam, checkIfExists: true)
                def bai_path = row.bai ? row.bai.toString() : row.bam.toString() + '.bai'
                def bai      = file(bai_path, checkIfExists: true)
                tuple(meta, bam, bai)
            }
            .multiMap { meta, bam, bai ->
                bams: tuple(meta, bam)
                bais: tuple(meta, bai)
            }

        ch_bam_final = ch_ss.bams
        ch_bai       = ch_ss.bais

    } else {
        // ------------------------------------------------------------------
        // FASTQ samplesheet: check for precomputed BAMs or align
        // ------------------------------------------------------------------
        def ss             = Utils.readSampleSheetInfo(params.input)
        def sampleIds      = ss.sampleIds
        def sample2patient = ss.sample2patient
        def bamDir         = params.existing_bam_dir.toString()

        def missingPairs    = Utils.findMissingBamBai(sampleIds, bamDir)
        boolean useExisting = false

        if (params.use_existing_bams as boolean) {
            if (missingPairs) {
                def msg = "Missing BAM/BAI in ${bamDir} for: ${missingPairs.join(', ')}"
                if (params.existing_bams_strict as boolean) error msg
                log.warn "${msg} — falling back to FASTQ alignment."
            } else {
                useExisting = true
                log.info "Reusing precomputed BAM/BAI from ${bamDir}."
            }
        }

        if (useExisting) {
            ch_bam_final = Channel.fromList(sampleIds).map { sid ->
                def meta = [ id: sid, sample: sid, patient: sample2patient[sid] ]
                tuple(meta, file("${bamDir}/${sid}.bam", checkIfExists: true))
            }
            ch_bai = Channel.fromList(sampleIds).map { sid ->
                def meta = [ id: sid, sample: sid, patient: sample2patient[sid] ]
                tuple(meta, file("${bamDir}/${sid}.bam.bai", checkIfExists: true))
            }

            def missingMetrics = Utils.findMissingMetrics(sampleIds, bamDir)
            if (missingMetrics) {
                log.warn "Missing *.metrics.txt for ${missingMetrics.size()} sample(s) — PICARDLIKE_DUPMETRICS skipped."
            } else {
                ch_dup_metrics = Channel.fromList(sampleIds).map { sid ->
                    def meta = [ id: sid, sample: sid, patient: sample2patient[sid] ]
                    tuple(meta, file("${bamDir}/${sid}.metrics.txt", checkIfExists: true))
                }
                run_picardlike = true
            }

        } else {
            def ch_reads = Channel
                .fromPath(params.input, checkIfExists: true)
                .splitCsv(header: true)
                .map { row ->
                    def sample = row.sample.toString()
                    def lane   = (row.lane ?: 'L1').toString()
                    def meta   = [ id: "${sample}_${lane}", sample: sample, patient: row.patient, lane: lane ]
                    tuple(meta, [ file(row.fastq_1), file(row.fastq_2) ])
                }

            ALIGN_AND_DEDUP(ch_reads, ch_index, ch_ref)

            ch_bam_final   = ALIGN_AND_DEDUP.out.bam
            ch_dup_metrics = ALIGN_AND_DEDUP.out.metrics
            run_picardlike = true

            ch_bai = ALIGN_AND_DEDUP.out.bam_bai.map { meta, bam, bai -> tuple(meta, bai) }
        }
    }

    // Convert biobambam metrics to Picard-like format for MultiQC
    if ((params.picardlike_enable == null || params.picardlike_enable) && run_picardlike && ch_dup_metrics != null) {
        PICARDLIKE_DUPMETRICS(ch_dup_metrics)
    }

    // Robust BAM + BAI join by meta.id → 3-tuple [meta, bam, bai]
    def ch_bam_bai = ch_bam_final
        .map { meta, bam -> tuple(meta.id, meta, bam) }
        .join(ch_bai.map { meta, bai -> tuple(meta.id, bai) })
        .map { id, meta, bam, bai -> tuple(meta, bam, bai) }

    // ------------------------------------------------------------------
    // QC
    // ------------------------------------------------------------------
    def ch_qualimap_bed = Channel.empty()

    if (params.qualimap_enable) {
        def ch_qualimap_bed_raw = params.qualimap_feature_bed ? Channel.value(file(params.qualimap_feature_bed))
                                : params.master_bed            ? Channel.value(file(params.master_bed))
                                : ch_genome_bed

        PREPARE_BED_QUALIMAP(ch_qualimap_bed_raw, ch_ref_fai_path)
        ch_qualimap_bed = PREPARE_BED_QUALIMAP.out.out_bed
    }

    def ch_vbid_svd    = Channel.value(tuple(
                            file(params.verifybamid2_ud),
                            file(params.verifybamid2_mu),
                            file(params.verifybamid2_bed)
                         ))
    def ch_vbid_refvcf = Channel.value(file(params.verifybamid2_dummy_refvcf))

    BAM_QC(
        ch_bam_bai,
        ch_qualimap_bed,
        ch_ref,
        ch_vbid_svd,
        ch_vbid_refvcf
    )

    // ------------------------------------------------------------------
    // SVDB annotation input channels (populated by each enabled tool below)
    // ------------------------------------------------------------------
    def ch_svdb_manta = Channel.empty()
    def ch_svdb_xhmm  = Channel.empty()
    def ch_svdb_exd   = Channel.empty()

    // ------------------------------------------------------------------
    // CNV: ExomeDepth (cohort-level)
    // ------------------------------------------------------------------
    if (params.exomedepth_enable) {
        if (!params.idxstats_enable) {
            error "exomedepth_enable=true requires idxstats_enable=true"
        }

        def ch_exd_bed_raw = params.exomedepth_target_bed ? Channel.value(file(params.exomedepth_target_bed))
                           : params.master_bed             ? Channel.value(file(params.master_bed))
                           : ch_genome_bed

        CNV_EXOMEDEPTH(
            ch_bam_bai,
            BAM_QC.out.idxstats,
            ch_exd_bed_raw,
            ch_ref_fasta_path,
            ch_ref_fai_path
        )
        ch_svdb_exd = CNV_EXOMEDEPTH.out.calls
            .flatten()
            .map { csv ->
                def sample = csv.name.replaceAll(/\.exomedepth\.cnv\.csv$/, '')
                tuple([id: sample, sample: sample], csv)
            }
    }

    // ------------------------------------------------------------------
    // CNV: DepthOfCoverage (standalone or as XHMM pre-requisite)
    // ------------------------------------------------------------------
    def ch_doc_coverage = Channel.empty()

    if (params.depthofcoverage_enable) {
        def ch_cov_bed_raw = params.coverage_target_bed ? Channel.value(file(params.coverage_target_bed))
                           : params.master_bed           ? Channel.value(file(params.master_bed))
                           : ch_genome_bed

        CNV_DEPTHOFCOVERAGE(
            ch_bam_bai,
            ch_cov_bed_raw,
            ch_ref_fasta_path,
            ch_ref_fai_path,
            ch_ref_dict_path,
            Channel.value(file(params.gatk_gene_list))
        )

        ch_doc_coverage = CNV_DEPTHOFCOVERAGE.out.coverage
    }

    // ------------------------------------------------------------------
    // CNV: XHMM (requires DepthOfCoverage)
    // ------------------------------------------------------------------
    if (params.xhmm_enable) {
        if (!params.depthofcoverage_enable) {
            error "xhmm_enable=true requires depthofcoverage_enable=true"
        }

        def ch_xhmm_bed_raw = params.xhmm_intervals_bed ? Channel.value(file(params.xhmm_intervals_bed))
                            : params.master_bed          ? Channel.value(file(params.master_bed))
                            : ch_genome_bed

        def ch_extreme_gc = Channel.value(
            params.xhmm_extreme_gc_targets
                ? file(params.xhmm_extreme_gc_targets)
                : file("${projectDir}/resources/xhmm/empty_targets.txt")
        )

        CNV_XHMM(
            ch_doc_coverage,
            ch_xhmm_bed_raw,
            ch_ref_fasta_path,
            ch_ref_fai_path,
            ch_ref_dict_path,
            Channel.value(file(params.xhmm_param_file)),
            ch_extreme_gc
        )
        ch_svdb_xhmm = CNV_XHMM.out.xcnv
    }

    // ------------------------------------------------------------------
    // Variant calling
    // ------------------------------------------------------------------
    if (params.deepvariant_enable || params.strelka_enable || params.manta_enable) {

        // Build +100 bp padded BED from the canonical target BED.
        // Stored in resources/ via storeDir — skipped if already present.
        def ch_canonical_bed = params.master_bed          ? Channel.value(file(params.master_bed))
                             : params.coverage_target_bed ? Channel.value(file(params.coverage_target_bed))
                             : ch_genome_bed

        PAD_BED(ch_canonical_bed, ch_ref_fai_path)
        def ch_padded_bed = PAD_BED.out.padded_bed

        VARIANT_CALLING(
            ch_bam_bai,
            ch_padded_bed,
            ch_padded_bed,
            ch_ref,
            ch_ref_fai,
            ch_ref_fasta_path,
            ch_ref_fai_path,
            Channel.value(tuple(meta_ref, file(params.deepvariant_gzi_dummy))),
            Channel.value(tuple(meta_ref, file(params.par_regions_bed))),
            Channel.value(file(params.manta_config))
        )

        if (params.manta_enable) {
            ch_svdb_manta = VARIANT_CALLING.out.manta_vcf
        }

        // ------------------------------------------------------------------
        // VCF merging / normalisation
        // ------------------------------------------------------------------
        VARIANT_MERGE(
            VARIANT_CALLING.out.dv_vcf,
            VARIANT_CALLING.out.strelka_vcf,
            ch_ref
        )
    }

    // ------------------------------------------------------------------
    // SVDB annotation (Manta SV, XHMM CNV, ExomeDepth CNV)
    //
    // Each input channel is built by mixing:
    //   a) live outputs from this run  (populated above when tool is enabled)
    //   b) Channel.fromPath scan of existing results/ (fallback when tool is
    //      disabled but a previous run already produced the files)
    //
    // If neither source has data the arm is a silent no-op.
    // ------------------------------------------------------------------
    if (params.svdb_enable) {
        def ch_m = ch_svdb_manta.mix(
            params.manta_enable ? Channel.empty()
            : Channel.fromPath("${params.outdir}/cnv/manta/**/*.diploid_sv.vcf.gz", checkIfExists: false)
                  .map { vcf -> tuple([id: vcf.parent.name, sample: vcf.parent.name], vcf) }
        )

        def ch_x = ch_svdb_xhmm.mix(
            params.xhmm_enable ? Channel.empty()
            : Channel.fromPath("${params.outdir}/cnv/xhmm/DATA.xcnv", checkIfExists: false)
                  .map { xcnv -> tuple([id: 'DATA'], xcnv) }
        )

        def ch_e = ch_svdb_exd.mix(
            params.exomedepth_enable ? Channel.empty()
            : Channel.fromPath("${params.outdir}/cnv/exomedepth/**/calls/*.exomedepth.cnv.csv", checkIfExists: false)
                  .map { csv ->
                      def sample = csv.name.replaceAll(/\.exomedepth\.cnv\.csv$/, '')
                      tuple([id: sample, sample: sample], csv)
                  }
        )

        SVDB_ANNOTATE(ch_m, ch_x, ch_e)
    }
}
