#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * -------------------------
 * Module includes
 * -------------------------
 */
include { BWA_MEM                       } from './modules/nf-core/bwa/mem'
include { BIOBAMBAM_BAMMARKDUPLICATES2 } from './modules/nf-core/biobambam/bammarkduplicates2'
include { SAMTOOLS_INDEX               } from './modules/nf-core/samtools/index'
include { FASTQC                       } from './modules/nf-core/fastqc'
include { SAMTOOLS_IDXSTATS            } from './modules/nf-core/samtools/idxstats'
include { VERIFYBAMID_VERIFYBAMID2     } from './modules/nf-core/verifybamid/verifybamid2'
include { QUALIMAP_BAMQC               } from './modules/nf-core/qualimap/bamqc'

include { GATK_DEPTHOFCOVERAGE         } from './modules/local/gatk/depthofcoverage'
include { GATK_ANNOTATEINTERVALS       } from './modules/local/gatk/annotateintervals'
include { XHMM_PIPELINE                } from './modules/local/xhmm'

include { DEEPVARIANT_RUNDEEPVARIANT   } from './modules/nf-core/deepvariant/rundeepvariant'
include { STRELKA_GERMLINE             } from './modules/nf-core/strelka/germline'
include { MANTA_GERMLINE               } from './modules/nf-core/manta/germline'

include { TABIX_TABIX                  } from './modules/nf-core/tabix/tabix'
include { STRELKA_CHRM_FILTER          } from './modules/local/strelka/chrm_filter'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF } from './modules/nf-core/bcftools/concat'
include { BCFTOOLS_NORM  as BCFTOOLS_NORM_COMBINED } from './modules/nf-core/bcftools/norm'
include { BCFTOOLS_NORM  as BCFTOOLS_NORM_DVONLY } from './modules/nf-core/bcftools/norm'

/*
 * Your requested interval modules
 */
include { BUILD_INTERVALS      } from './modules/local/build_intervals'
include { CREATE_INTERVALS_BED } from './modules/local/create_intervals_bed'

/*
 * BED preparation module (imported multiple times with aliases because DSL2
 * does not allow reusing the same included process name in one workflow scope).
 */
include { PREPARE_BED as PREPARE_BED_QUALIMAP } from './modules/local/prepare_bed'
include { PREPARE_BED as PREPARE_BED_COVERAGE } from './modules/local/prepare_bed'
include { PREPARE_BED as PREPARE_BED_XHMM } from './modules/local/prepare_bed'
include { PREPARE_BED as PREPARE_BED_MASTER } from './modules/local/prepare_bed'
include { PREPARE_BED as PREPARE_BED_VARIANT } from './modules/local/prepare_bed'


/*
 * -------------------------
 * Helpers
 * -------------------------
 */
def getGenomeAttr(String attr) {
  if (!params.genome) return null
  if (!params.genomes || !params.genomes.containsKey(params.genome)) {
    error "Unknown --genome '${params.genome}'. Add it to conf/igenomes.config"
  }
  def g = params.genomes[params.genome]
  if (!g.containsKey(attr)) {
    error "Genome '${params.genome}' does not define '${attr}' in conf/igenomes.config"
  }
  return g[attr]
}


/*
 * -------------------------
 * Local processes
 * -------------------------
 */

// Add read groups at lane-level
process ADD_READGROUPS {
  tag "${meta.id}"
  label 'process_small'
  container 'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6'

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("${meta.id}.rg.bam"), emit: bam

  script:
  def rgid     = meta.id
  def sample   = meta.sample ?: meta.id
  def lib      = meta.library ?: meta.sample ?: meta.id
  def platform = meta.platform ?: 'ILLUMINA'
  def pu       = meta.platform_unit ?: "${meta.sample}.${meta.lane ?: 'L1'}"

  """
  set -euo pipefail
  picard AddOrReplaceReadGroups \
    I=${bam} \
    O=${meta.id}.rg.bam \
    RGID=${rgid} \
    RGLB=${lib} \
    RGPL=${platform} \
    RGPU=${pu} \
    RGSM=${sample}
  """
}


// Merge lane BAMs per sample, then output a single per-sample BAM
process MERGE_BAMS {
  tag "${meta.id}"
  label 'process_medium'
  container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

  input:
  tuple val(meta), path(bams)

  output:
  tuple val(meta), path("${meta.id}.merged.bam"), emit: bam

  script:
  def nbams = bams.size()
  def bam_list = bams.collect { "\"${it}\"" }.join(' ')
  def first_bam = "\"${bams[0]}\""

  """
  set -euo pipefail

  if [ ${nbams} -eq 1 ]; then
    ln -s ${first_bam} ${meta.id}.merged.bam
  else
    samtools merge -@ ${task.cpus} -O BAM - ${bam_list} | \
      samtools sort  -@ ${task.cpus} -o ${meta.id}.merged.bam -
  fi
  """
}


// Make biobambam metrics “Picard-like” for MultiQC
process PICARDLIKE_DUPMETRICS {
  tag "${meta.id}"
  label 'process_single'
  container 'quay.io/biocontainers/coreutils:9.5--a9c29d0e2be5c7b6'  // just needs tail/cat

  input:
  tuple val(meta), path(metrics)

  output:
  tuple val(meta), path("${meta.id}.duplicate_metrics_picard_like.txt"), emit: metrics

  script:
  """
  set -euo pipefail
  tail -n +4 ${metrics} > body.txt

  cat > header.txt <<'EOF'
##htsjdk.samtools.metrics.StringHeader
##METRICS CLASS picard.sam.DuplicationMetrics
EOF

  cat header.txt body.txt > ${meta.id}.duplicate_metrics_picard_like.txt
  """
}


// bgzip a BED
process BGZIP_BED {
  tag "${bed.baseName}"
  label 'process_single'
  container 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'

  input:
  path bed

  output:
  path "${bed.baseName}.bed.gz", emit: bedgz

  script:
  """
  set -euo pipefail
  bgzip -c ${bed} > ${bed.baseName}.bed.gz
  """
}




/*
 * Gather scattered DV chunk VCFs into one per sample
 */
process GATHER_DEEPVARIANT_VCFS {
  tag "${sample_id}"
  label 'process_small'
  container 'quay.io/biocontainers/bcftools:1.22--h3a4d415_2'

  input:
  val sample_id
  path vcfs
  path tbis

  output:
  tuple val([id: sample_id]), path("${sample_id}.deepvariant.vcf.gz"), path("${sample_id}.deepvariant.vcf.gz.tbi"), emit: vcf

  script:
  """
  set -euo pipefail
  printf "%s\\n" ${vcfs} | LC_ALL=C sort -V > vcfs.list
  bcftools concat -a -D -Oz -o ${sample_id}.deepvariant.vcf.gz -f vcfs.list
  bcftools index -t ${sample_id}.deepvariant.vcf.gz
  """
}


/*
 * -------------------------
 * Workflow
 * -------------------------
 */
workflow {

  if (!params.genome || params.igenomes_ignore) {
    error "This pipeline is iGenomes-only. Provide --genome (e.g. GATK.GRCh38) and do not set --igenomes_ignore."
  }

  def meta_ref = [ id: params.genome ]

  // Reference paths
  def fasta_path = file(getGenomeAttr('fasta'))
  def bwa_path   = file(getGenomeAttr('bwa'))
  def fai_path   = file(getGenomeAttr('fasta_fai'))
  def dict_path  = file(getGenomeAttr('dict'))

  // Reference channels
  Channel.value(tuple(meta_ref, fasta_path)).set { ref_val }
  Channel.value(tuple(meta_ref, bwa_path)).set   { index_val }
  Channel.value(tuple(meta_ref, fai_path)).set   { fai_val }

  Channel.value(fasta_path).set { ch_ref_fasta_path }
  Channel.value(fai_path).set   { ch_ref_fai_path }
  Channel.value(dict_path).set  { ch_ref_dict_path }

  /*
   * Case B base: build genome.bed from .fai using your module
   */
    BUILD_INTERVALS(fai_val)
    BUILD_INTERVALS.out.bed
    .map { meta, bed -> bed }
    .set { ch_genome_bed }


  /*
   * -------------------------
   * Input samplesheet: patient,sample,lane,fastq_1,fastq_2
   * meta.id unique per row (lane); meta.sample preserved for biological sample
   * -------------------------
   */
  Channel
    .fromPath(params.input, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
      def sample = row.sample.toString()
      def lane   = (row.lane ?: 'L1').toString()

      def meta = [
        id      : "${sample}_${lane}",
        sample  : sample,
        patient : row.patient,
        lane    : lane
      ]

      tuple(meta, [ file(row.fastq_1), file(row.fastq_2) ])
    }
    .set { ch_reads }

  /*
   * -------------------------
   * Alignment → lane BAM
   * -------------------------
   */
  BWA_MEM(ch_reads, index_val, ref_val, true)
  BWA_MEM.out.bam.set { ch_bam_lane }

  /*
   * Add read groups per lane
   */
  ADD_READGROUPS(ch_bam_lane)
  ADD_READGROUPS.out.bam.set { ch_bam_rg_lane }

  /*
   * Merge per sample (group by meta.sample), then mark duplicates once
   */
  ch_bam_rg_lane
    .map { meta, bam -> tuple(meta.sample, meta, bam) }
    .groupTuple()
    .map { sample, metas, bams ->
      def m0 = metas[0]
      def merged_meta = [ id: sample, sample: sample, patient: m0.patient ]
      tuple(merged_meta, bams)
    }
    .set { ch_merge_in }

  MERGE_BAMS(ch_merge_in)
  MERGE_BAMS.out.bam.set { ch_bam_merged }

  BIOBAMBAM_BAMMARKDUPLICATES2(ch_bam_merged)
  BIOBAMBAM_BAMMARKDUPLICATES2.out.bam.set     { ch_bam_final }
  BIOBAMBAM_BAMMARKDUPLICATES2.out.metrics.set { ch_dup_metrics }

  SAMTOOLS_INDEX(ch_bam_final)
  SAMTOOLS_INDEX.out.bai.set { ch_bai }

  PICARDLIKE_DUPMETRICS(ch_dup_metrics)

  /*
   * Robust BAM/BAI join by meta.id
   */
  ch_bam_final.map { meta, bam -> tuple(meta.id, meta, bam) }.set { ch_bam_keyed }
  ch_bai.map      { meta, bai -> tuple(meta.id, bai) }.set        { ch_bai_keyed }

  ch_bam_keyed
    .join(ch_bai_keyed)
    .map { id, meta, bam, bai -> tuple(meta, bam, bai) }
    .set { ch_bam_bai }

  /*
   * -------------------------
   * QC toggles
   * -------------------------
   */
  if (params.fastqc_enable) {
    FASTQC(ch_bam_final)
  }

  if (params.idxstats_enable) {
    SAMTOOLS_IDXSTATS(ch_bam_bai)
  }

  if (params.verifybamid2_enable) {
    Channel.value(tuple(file(params.verifybamid2_ud), file(params.verifybamid2_mu), file(params.verifybamid2_bed))).set { ch_vbid_svd }
    Channel.value(file(params.verifybamid2_dummy_refvcf)).set { ch_vbid_refvcf }
    ref_val.map { meta, fasta -> fasta }.set { ch_vbid_ref }
    VERIFYBAMID_VERIFYBAMID2(ch_bam_bai, ch_vbid_svd, ch_vbid_refvcf, ch_vbid_ref)
  }

  if (params.qualimap_enable) {
    // Qualimap feature BED: if not provided, fall back to master_bed else genome bed
    def role_bed_raw = params.master_bed ? Channel.value(file(params.master_bed)) : ch_genome_bed
    def qualimap_raw = params.qualimap_feature_bed ? Channel.value(file(params.qualimap_feature_bed)) : role_bed_raw

    PREPARE_BED_QUALIMAP(qualimap_raw, ch_ref_fai_path)
    def ch_qualimap_bed = PREPARE_BED_QUALIMAP.out.out_bed

    QUALIMAP_BAMQC(ch_bam_final, ch_qualimap_bed)
  }

  /*
   * -------------------------
   * DepthOfCoverage toggle
   * -------------------------
   */
  if (params.depthofcoverage_enable) {
    def role_bed_raw = params.master_bed ? Channel.value(file(params.master_bed)) : ch_genome_bed
    def cov_raw      = params.coverage_target_bed ? Channel.value(file(params.coverage_target_bed)) : role_bed_raw

    PREPARE_BED_COVERAGE(cov_raw, ch_ref_fai_path)
    def ch_cov_bed = PREPARE_BED_COVERAGE.out.out_bed

    Channel.value(file(params.gatk_gene_list)).set { ch_gene_list }

    GATK_DEPTHOFCOVERAGE(
      ch_bam_bai,
      ch_ref_fasta_path,
      ch_ref_fai_path,
      ch_ref_dict_path,
      ch_cov_bed,
      ch_gene_list
    )
  }

  /*
   * -------------------------
   * XHMM toggle (requires DepthOfCoverage outputs)
   * -------------------------
   */
  if (params.xhmm_enable) {
    if (!params.depthofcoverage_enable) {
      error "xhmm_enable=true requires depthofcoverage_enable=true (XHMM needs DepthOfCoverage outputs)."
    }

    def meta_xhmm = [ id: (params.xhmm_prefix ?: 'DATA') ]

    def ch_doc_coverage = GATK_DEPTHOFCOVERAGE.out.coverage

    def ch_doc_sample_interval_summary = ch_doc_coverage
      .map { meta, files ->
        def sis = files.find { it.name.endsWith('.sample_interval_summary') }
        if (!sis) error "DepthOfCoverage outputs for '${meta.id}' do not contain a .sample_interval_summary file"
        tuple(meta, sis)
      }

    def ch_xhmm_gatk_summaries = ch_doc_sample_interval_summary
      .map { meta, sis -> sis }
      .collect()
      .map { summaries -> tuple(meta_xhmm, summaries) }

    def xhmm_raw = params.xhmm_intervals_bed ? Channel.value(file(params.xhmm_intervals_bed))
                                            : (params.master_bed ? Channel.value(file(params.master_bed)) : ch_genome_bed)

    PREPARE_BED_XHMM(xhmm_raw, ch_ref_fai_path)
    def ch_xhmm_intervals_bed = PREPARE_BED_XHMM.out.out_bed

    def ch_xhmm_intervals_for_annot = ch_xhmm_intervals_bed.map { bed -> tuple(meta_xhmm, bed) }

    GATK_ANNOTATEINTERVALS(ch_xhmm_intervals_for_annot, ch_ref_fasta_path, ch_ref_fai_path, ch_ref_dict_path)

    def ch_xhmm_annotated = GATK_ANNOTATEINTERVALS.out.tsv.map { meta, tsv -> tsv }

    def ch_xhmm_params = Channel.value(file(params.xhmm_param_file))
    def ch_xhmm_extreme_gc = Channel.value(
      params.xhmm_extreme_gc_targets ? file(params.xhmm_extreme_gc_targets) : file("${baseDir}/resources/xhmm/empty_targets.txt")
    )

    XHMM_PIPELINE(ch_xhmm_gatk_summaries, ch_xhmm_annotated, ch_xhmm_extreme_gc, ch_xhmm_params)
  }

  /*
   * -------------------------
   * Variant calling toggles
   * -------------------------
   */
  def run_dv      = params.deepvariant_enable as boolean
  def run_strelka = params.strelka_enable as boolean
  def run_manta   = params.manta_enable as boolean

  def ch_dv_vcf      = null
  def ch_strelka_vcf = null

  // ---- DeepVariant (Case A/B intervals + scatter/gather) ----
  if (run_dv) {

    // Case A/B master intervals selection
    def ch_master_intervals_raw
    if (params.intervals) {
      ch_master_intervals_raw = Channel.value(file(params.intervals))
    } else if (!params.no_intervals) {
      ch_master_intervals_raw = ch_genome_bed
    } else {
      ch_master_intervals_raw = ch_genome_bed
    }

    PREPARE_BED_MASTER(ch_master_intervals_raw, ch_ref_fai_path)
    def ch_master_intervals = PREPARE_BED_MASTER.out.out_bed


    def ch_interval_chunks
    if (!params.no_intervals) {
    CREATE_INTERVALS_BED(ch_master_intervals, ch_ref_fai_path)
    ch_interval_chunks = CREATE_INTERVALS_BED.out.bed
    } else {
    ch_interval_chunks = ch_master_intervals
    }


    Channel.value(tuple(meta_ref, file(params.deepvariant_gzi_dummy))).set { ch_deep_gzi }
    Channel.value(tuple(meta_ref, file(params.par_regions_bed))).set       { ch_deep_par }

    ch_bam_bai
      .combine(ch_interval_chunks)
      .map { meta, bam, bai, bed ->
        def m = meta.clone()
        m.sample_id = meta.id
        m.id = "${meta.id}__${bed.baseName}"
        tuple(m, bam, bai, bed)
      }
      .set { ch_dv_scatter_in }

    DEEPVARIANT_RUNDEEPVARIANT(ch_dv_scatter_in, ref_val, fai_val, ch_deep_gzi, ch_deep_par)

    DEEPVARIANT_RUNDEEPVARIANT.out.vcf
      .join(DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index)
      .map { meta, vcf, tbi -> tuple(meta.sample_id, vcf, tbi) }
      .groupTuple()
      .map { sample_id, vcfs, tbis -> tuple(sample_id, vcfs, tbis) }
      .set { ch_dv_gather_in }

    GATHER_DEEPVARIANT_VCFS(ch_dv_gather_in)
    ch_dv_vcf = GATHER_DEEPVARIANT_VCFS.out.vcf
  }

  // ---- Prepare target BED.gz + .tbi once for Strelka/Manta ----
  def ch_target_bed_gz = null
  def ch_target_bed_tbi = null

  if (run_strelka || run_manta) {
    def role_bed_raw = params.master_bed ? Channel.value(file(params.master_bed)) : ch_genome_bed
    def var_raw      = params.variant_target_bed ? Channel.value(file(params.variant_target_bed)) : role_bed_raw

    PREPARE_BED_VARIANT(var_raw, ch_ref_fai_path)
    def ch_variant_bed = PREPARE_BED_VARIANT.out.out_bed

    BGZIP_BED(ch_variant_bed)
    ch_target_bed_gz = BGZIP_BED.out.bedgz

    ch_target_bed_gz
      .map { bedgz -> tuple([id: 'variant_bed'], bedgz) }
      .set { ch_for_tabix }

    TABIX_TABIX(ch_for_tabix)
    ch_target_bed_tbi = TABIX_TABIX.out.index.map { meta, idx -> idx }
  }

  // ---- Strelka2 ----
  if (run_strelka) {
    ch_bam_bai
      .combine(ch_target_bed_gz)
      .combine(ch_target_bed_tbi)
      .map { meta, bam, bai, bedgz, bedidx -> tuple(meta, bam, bai, bedgz, bedidx) }
      .set { ch_strelka_in }

    STRELKA_GERMLINE(ch_strelka_in, ch_ref_fasta_path, ch_ref_fai_path)

    STRELKA_GERMLINE.out.vcf
      .join(STRELKA_GERMLINE.out.vcf_tbi)
      .set { ch_strelka_with_index }

    STRELKA_CHRM_FILTER(ch_strelka_with_index)
    ch_strelka_vcf = STRELKA_CHRM_FILTER.out.vcf.join(STRELKA_CHRM_FILTER.out.tbi)
  }

  // ---- Manta ----
  if (run_manta) {
    def ch_manta_samples = ch_bam_bai
      .combine(ch_target_bed_gz)
      .combine(ch_target_bed_tbi)
      .map { meta, bam, bai, bedgz, bedidx -> tuple(meta, [bam], bai, bedgz, bedidx) }

    Channel.value(file(params.manta_config)).set { ch_manta_config }
    MANTA_GERMLINE(ch_manta_samples, ref_val, fai_val, ch_manta_config)
  }

  // ---- Combine DV + Strelka ----
  if (params.combine_dv_strelka_enable) {
    if (!(run_dv && run_strelka)) {
      if (params.strict_variant_requirements) {
        error "combine_dv_strelka_enable=true requires deepvariant_enable=true and strelka_enable=true"
      } else {
        log.warn "Skipping DV+Strelka combine: requires deepvariant_enable && strelka_enable"
      }
    } else {
      def ch_concat_in = ch_dv_vcf
        .join(ch_strelka_vcf)
        .map { meta, dv_vcf, dv_tbi, st_vcf, st_tbi ->
          def m = meta.clone()
          m.id = "${meta.id}.DV_ST2.comb"
          tuple(m, [dv_vcf, st_vcf], [dv_tbi, st_tbi])
        }

      BCFTOOLS_CONCAT_VCF(ch_concat_in)

      def ch_combined_vcf = BCFTOOLS_CONCAT_VCF.out.vcf.join(BCFTOOLS_CONCAT_VCF.out.tbi)

      if (params.norm_combined_enable) {
        def ch_norm_in = ch_combined_vcf
          .map { meta, vcf, tbi ->
            def m = meta.clone()
            m.id = "${meta.id}.DV_ST2.comb.norm"
            tuple(m, vcf, tbi)
          }
        BCFTOOLS_NORM_COMBINED(ch_norm_in, ref_val)
      }
    }
  }

  // ---- DV-only normalization ----
  if (run_dv && params.norm_dv_only_enable) {
    def ch_dv_norm_in = ch_dv_vcf
      .map { meta, vcf, tbi ->
        def m = meta.clone()
        m.id = "${meta.id}.DV_only.norm"
        tuple(m, vcf, tbi)
      }

    BCFTOOLS_NORM_DVONLY(ch_dv_norm_in, ref_val)
  }
}