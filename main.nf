#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Include nf-core modules
 */
include { BWA_INDEX                        } from './modules/nf-core/bwa/index/main'
include { BWA_MEM                          } from './modules/nf-core/bwa/mem/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from './modules/nf-core/picard/createsequencedictionary/main'
include { BIOBAMBAM_BAMMARKDUPLICATES2    } from './modules/nf-core/biobambam/bammarkduplicates2/main'
include { SAMTOOLS_INDEX                  } from './modules/nf-core/samtools/index/main'
include { FASTQC                          } from './modules/nf-core/fastqc/main'
include { SAMTOOLS_IDXSTATS               } from './modules/nf-core/samtools/idxstats/main'
include { VERIFYBAMID_VERIFYBAMID2        } from './modules/nf-core/verifybamid/verifybamid2/main'
include { QUALIMAP_BAMQC                  } from './modules/nf-core/qualimap/bamqc/main'
include { GATK_DEPTHOFCOVERAGE            } from './modules/local/gatk/depthofcoverage/main'
include { DEEPVARIANT_RUNDEEPVARIANT } from './modules/nf-core/deepvariant/rundeepvariant/main'
include { STRELKA_GERMLINE } from './modules/nf-core/strelka/germline/main'
include { TABIX_TABIX      } from './modules/nf-core/tabix/tabix/main'
include { STRELKA_CHRM_FILTER } from './modules/local/strelka/chrm_filter/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF      } from './modules/nf-core/bcftools/concat'
include { BCFTOOLS_NORM  as BCFTOOLS_NORM_COMBINED   } from './modules/nf-core/bcftools/norm'
include { BCFTOOLS_NORM  as BCFTOOLS_NORM_DVONLY     } from './modules/nf-core/bcftools/norm'
include { MANTA_GERMLINE } from './modules/nf-core/manta/germline/main'
include { GATK_ANNOTATEINTERVALS } from './modules/local/gatk/annotateintervals/main'
include { XHMM_PIPELINE         } from './modules/local/xhmm/main'




def getGenomeAttr(String attr) {
  if ( !params.genome ) return null
  if ( !params.genomes || !params.genomes.containsKey(params.genome) )
    error "Unknown --genome '${params.genome}'. Add it to conf/igenomes.config"
  def g = params.genomes[params.genome]
  if ( !g.containsKey(attr) )
    error "Genome '${params.genome}' does not define '${attr}' in conf/igenomes.config"
  return g[attr]
}




/*
 * SAMTOOLS_FAIDX_SIMPLE
 * (simple faidx wrapper for the reference FASTA)
 */
process SAMTOOLS_FAIDX_SIMPLE {

    tag "${meta.id}"
    publishDir "${params.outdir}/ref", mode: 'copy'
    container 'quay.io/biocontainers/samtools:1.3.1--h0cf4675_11'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta}.fai")

    script:
    """
    set -euo pipefail
    samtools faidx ${fasta}
    """
}


/*
 * MAKE_GENOME_BED
 * Build a simple whole-genome BED from the FASTA index (.fai)
 * One line per contig: chrom, 0, length
 */
process MAKE_GENOME_BED {

    tag "${meta.id}"
    publishDir "${params.outdir}/ref", mode: 'copy'

    input:
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path("genome.autobed")

    script:
    """
    set -euo pipefail
    awk 'BEGIN{OFS="\\t"} {print \$1,0,\$2}' ${fai} > genome.autobed
    """
}


/*
 * ADD_READGROUPS
 */
process ADD_READGROUPS {

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.rg.bam")

    script:
    def rgid     = meta.id
    def sample   = meta.sample   ?: meta.id
    def lib      = meta.library  ?: meta.patient ?: meta.id
    def platform = meta.platform ?: 'ILLUMINA'
    def pu       = meta.platform_unit ?: "${meta.id}.${meta.lane ?: 'L1'}"

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


/*
 * PICARDLIKE_DUPMETRICS
 * Wrap biobambam2 metrics into Picard-like format for MultiQC.
 */
process PICARDLIKE_DUPMETRICS {

    tag "${meta.id}"

    input:
    tuple val(meta), path(metrics)

    output:
    tuple val(meta), path("*.duplicate_metrics_picard_like.txt")

    script:
    """
    set -euo pipefail

    in_metrics="${metrics}"
    prefix="${meta.id}"

    # Drop the first 3 lines and keep the tabular metrics
    tail -n +4 "\${in_metrics}" > tmp.txt

    # Build a Picard-like header
    cat <<EOF > header.txt
##htsjdk.samtools.metrics.StringHeader
# MarkDuplicates INPUT=\${prefix}.bam OUTPUT=\${prefix}.bam METRICS_FILE=\${prefix}.duplication_metrics.txt ...
##htsjdk.samtools.metrics.StringHeader
##METRICS CLASS picard.sam.DuplicationMetrics
EOF

    cat header.txt tmp.txt > "\${prefix}.duplicate_metrics_picard_like.txt"
    rm tmp.txt header.txt
    """
}



/*
 * BGZIP_BED
 * Compress a BED file to BED.GZ (needed for tabix indexing)
 */
process BGZIP_BED {

    tag "bgzip_bed"
    container 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'

    input:
    path bed

    output:
    path "*.bed.gz"

    """
    bgzip -c ${bed} > \$(basename ${bed}).gz
    """
}



/*
 * Main workflow
 */
workflow {


    /*
    * 1) Reference
    */
    def meta_ref = [ id: params.genome ?: params.ref_id ]

    if ( params.genome && !params.igenomes_ignore ) {

        // Use iGenomes assets directly (no download, no indexing)
        def fasta_path = file(getGenomeAttr('fasta'))
        def bwa_path   = file(getGenomeAttr('bwa'))
        def fai_path   = file(getGenomeAttr('fasta_fai'))
        def dict_path  = file(getGenomeAttr('dict'))

        Channel.value( tuple(meta_ref, fasta_path) ).set { ref_val }     // [meta, fasta]
        Channel.value( tuple(meta_ref, bwa_path)   ).set { index_val }   // [meta, bwa index dir]
        Channel.value( tuple(meta_ref, fai_path)   ).set { fai_val }     // [meta, fai]

        // Path-only channels for DepthOfCoverage module inputs:
        Channel.value( fai_path  ).set { ch_doc_fai }
        Channel.value( dict_path ).set { ch_doc_dict }

    } else {
    error "This pipeline is iGenomes-only. Use --genome (e.g. GATK.GRCh38) and do not set --igenomes_ignore."
    }

    /*
    * 1b) Qualimap regions (must run for both iGenomes and fallback)
    */
    def ch_qualimap_regions

    if ( params.qualimap_bed ) {
        ch_qualimap_regions = Channel.value( file(params.qualimap_bed) )
    } else {
        MAKE_GENOME_BED(fai_val)  // uses (meta, fai) tuple
        ch_qualimap_regions = MAKE_GENOME_BED.out.map { meta, bed -> bed }
    }





    /*
     * 2) Samples from samplesheet.csv
     *    patient,sample,lane,fastq_1,fastq_2
     */
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id      : row.sample,
                sample  : row.sample,
                patient : row.patient,
                lane    : row.lane
            ]
            def reads = [ file(row.fastq_1), file(row.fastq_2) ]
            tuple(meta, reads)
        }
        .set { ch_reads }


    /*
     * 3) BWA-MEM alignment (sorted BAM)
     *    nf-core BWA_MEM internally does: bwa mem | samtools sort
     */
    BWA_MEM(
        ch_reads,
        index_val,  // value channel
        ref_val,    // value channel
        true        // sort_bam = true
    )

    BWA_MEM.out.bam.set { ch_bam_sorted }   // [meta_sample, sorted.bam]


    /*
     * 4) Add read groups
     */
    ADD_READGROUPS(ch_bam_sorted)
    ADD_READGROUPS.out.set { ch_bam_rg }    // [meta_sample, bam_with_rg]


    /*
     * 5) Mark duplicates (biobambam2, remove dups)
     */
    ch_bam_rg
        .map { meta, bam ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_rmdup"  // prefix for final BAM
            tuple(new_meta, bam)
        }
        .set { ch_bam_for_markdup }

    BIOBAMBAM_BAMMARKDUPLICATES2(ch_bam_for_markdup)

    BIOBAMBAM_BAMMARKDUPLICATES2.out.bam.set     { ch_bam_dedup }     // [meta_dedup, dedup.bam]
    BIOBAMBAM_BAMMARKDUPLICATES2.out.metrics.set { ch_dup_metrics }   // [meta_dedup, *.metrics.txt]


    /*
     * 6) Index final BAMs
     */
    SAMTOOLS_INDEX(ch_bam_dedup)
    SAMTOOLS_INDEX.out.bai.set { ch_bai }   // [meta_dedup, *.bai]


    /*
     * 7) Picard-like duplication metrics (for MultiQC)
     */
    PICARDLIKE_DUPMETRICS(ch_dup_metrics)


    /*
     * 8a) BAM-level FastQC on deduplicated BAMs
     */
    FASTQC(ch_bam_dedup)


    /*
     * 8b) samtools idxstats (BAM + BAI)
     */
    ch_bam_dedup
        .join(ch_bai)                       // [meta, bam, bai]
        .set { ch_bam_bai }

    SAMTOOLS_IDXSTATS(ch_bam_bai)


    /*
     * 9) VerifyBamID2 contamination estimates
     *    Using precomputed SVD (UD / mu / bed) + reference FASTA.
     *    No RefVCF (dummy placeholder).
     */

    // 9.1 SVD triple as a value channel
    Channel.value(
        tuple(
            file(params.verifybamid2_ud),
            file(params.verifybamid2_mu),
            file(params.verifybamid2_bed)
        )
    ).set { ch_vbid_svd }

    // 9.2 Dummy "refvcf" path that is NOT a .vcf
    Channel.value(
        file(params.verifybamid2_dummy_refvcf)
    ).set { ch_vbid_refvcf }

    // 9.3 Reference FASTA path as value channel
    ref_val
        .map { meta, fasta -> fasta }
        .set { ch_vbid_ref }

    // 9.4 Use BAM+BAI channel directly
    VERIFYBAMID_VERIFYBAMID2(
        ch_bam_bai,     // (meta, bam, bai)
        ch_vbid_svd,    // (UD, mu, bed)
        ch_vbid_refvcf, // dummy refvcf (no --RefVCF)
        ch_vbid_ref     // reference FASTA
    )

    // Capture selfSM (FREEMIX etc.)
    VERIFYBAMID_VERIFYBAMID2.out.self_sm.set { ch_verifybamid_selfSM }


    /*
     * 10) Qualimap2 (bamqc)
     *     - If WGS: uses auto-generated genome.autobed
     *     - If WES/targeted: uses params.qualimap_bed (BED)
     */
    QUALIMAP_BAMQC(
        ch_bam_dedup,       // (meta, bam)
        ch_qualimap_regions // single value channel with BED / regions file
    )

    QUALIMAP_BAMQC.out.results.set { ch_qualimap_results }


    /*
     * 11) GATK DepthOfCoverage (local module, GATK4)
     *     Using deduplicated BAM + BAI, reference FASTA/FAI/DICT, BED, gene list.
     */

    // Reference FASTA as value channel (just the path, not meta)
    ref_val
        .map { meta, fasta -> fasta }
        .set { ch_doc_ref }

    // ch_doc_fai and ch_doc_dict are defined above from fai_val and PICARD output


    // BED for coverage: prefer dedicated GATK BED, otherwise reuse qualimap_bed
    def bed_path = params.gatk_target_bed ?: params.qualimap_bed
    if( !bed_path ) {
        MAKE_GENOME_BED(fai_val)
        MAKE_GENOME_BED.out
            .first()                   // convert to value channel
            .map { meta, bed -> bed }
            .set { ch_doc_bed }
    } else {
        Channel.value( file(bed_path) ) // <-- value channel, broadcasts to all samples
            .set { ch_doc_bed }
    }

    // Gene list for DepthOfCoverage (e.g. refGene.hg38.txt)
    Channel.value( file(params.gatk_gene_list) ) // <-- value channel
        .set { ch_doc_gene_list }


    // Finally call the module
    GATK_DEPTHOFCOVERAGE(
        ch_bam_bai,      // (meta, bam, bai)
        ch_doc_ref,      // FASTA
        ch_doc_fai,      // FASTA.fai
        ch_doc_dict,     // FASTA.dict
        ch_doc_bed,      // BED
        ch_doc_gene_list // gene list
    )
    

    GATK_DEPTHOFCOVERAGE.out.coverage.set { ch_doc_coverage }


    /*
    * 11b) XHMM CNV pipeline (cohort-level)
    *      Requires DepthOfCoverage to output TABLE (tab/space delimited)
    */
if ( params.xhmm_enable ) {

    def meta_xhmm = [ id: (params.xhmm_prefix ?: 'DATA') ]

    // Extract per-sample .sample_interval_summary
    def ch_doc_sample_interval_summary = ch_doc_coverage
        .map { meta, files ->
            def sis = files.find { it.name.endsWith('.sample_interval_summary') }
            if( !sis ) error "DepthOfCoverage outputs for '${meta.id}' do not contain a .sample_interval_summary file"
            tuple(meta, sis)
        }

    // Collect into one list (cohort-level run)
    def ch_xhmm_gatk_summaries = ch_doc_sample_interval_summary
        .map { meta, sis -> sis }
        .collect()
        .map { summaries -> tuple(meta_xhmm, summaries) }

    def ch_xhmm_params = Channel.value( file(params.xhmm_param_file) )

    // Always pass a file (either user-provided or empty placeholder)
    def ch_xhmm_extreme_gc = Channel.value(
        params.xhmm_extreme_gc_targets
            ? file(params.xhmm_extreme_gc_targets)
            : file("${baseDir}/resources/xhmm/empty_targets.txt")
    )

    // Always run AnnotateIntervals (simple + avoids “optional input” problems)
    def ch_xhmm_intervals_for_annot = Channel.value(
        tuple(meta_xhmm, file(params.xhmm_intervals))
    )

    GATK_ANNOTATEINTERVALS(
        ch_xhmm_intervals_for_annot,
        ch_doc_ref,
        ch_doc_fai,
        ch_doc_dict
    )

    def ch_xhmm_annotated = GATK_ANNOTATEINTERVALS.out.tsv.map { meta, tsv -> tsv }

    XHMM_PIPELINE(
        ch_xhmm_gatk_summaries,
        ch_xhmm_annotated,
        ch_xhmm_extreme_gc,
        ch_xhmm_params
    )
}






    /*
     * 12) DeepVariant (variant calling)
     *     Uses final deduplicated BAM + BAI, reference FASTA/FAI, target regions BED,
     *     dummy GZI and PAR BED.
     */

    // 12.1 Target regions for DeepVariant: use params.target_regions_bed
    Channel.value( file(params.target_regions_bed) )
        .set { ch_target_regions_bed }

    // 12.2 Combine BAM+BAI with target intervals -> (meta, bam, bai, intervals)
    ch_bam_bai
        .combine(ch_target_regions_bed)
        .map { meta, bam, bai, intervals ->
            tuple(meta, bam, bai, intervals)
        }
        .set { ch_deep_bam }

    // 12.3 Reference FASTA and FAI as value channels
    def ch_deep_fasta = ref_val   // (meta_ref, fasta)
    def ch_deep_fai   = fai_val   // (meta_ref, fasta.fai)

    // 12.4 Dummy GZI for DeepVariant (must NOT be the .fai file)
    Channel.value(tuple( meta_ref, file(params.deepvariant_gzi_dummy) )).set { ch_deep_gzi }

    // 12.5 PAR BED (can be empty, no effect if empty)
    Channel.value(tuple( meta_ref, file(params.par_regions_bed) )).set { ch_deep_par }

    // 12.6 Run DeepVariant with all 5 required input tuples
    DEEPVARIANT_RUNDEEPVARIANT(
        ch_deep_bam,    // (meta, bam, bai, intervals)
        ch_deep_fasta,  // (meta_ref, fasta)
        ch_deep_fai,    // (meta_ref, fasta.fai)
        ch_deep_gzi,    // (meta_ref, dummy_fasta.gzi)
        ch_deep_par     // (meta_ref, par_regions.bed)
    )

    // 12.7 Capture DeepVariant outputs
    DEEPVARIANT_RUNDEEPVARIANT.out.vcf       .set { ch_dv_vcf_raw }   // (meta, vcf)
    DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index .set { ch_dv_vcf_tbi }   // (meta, tbi)
    DEEPVARIANT_RUNDEEPVARIANT.out.gvcf      .set { ch_dv_gvcf }
    DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index.set { ch_dv_gvcf_index }

    // Join DV VCF + index -> (meta, vcf, tbi)
    def ch_dv_vcf = ch_dv_vcf_raw.join(ch_dv_vcf_tbi)



    /*
     * 12b) Prepare compressed BED + tabix index for Strelka
     *      Starting from the plain target BED (params.target_regions_bed)
     */

    // Plain BED as a value channel for Strelka
    Channel.value( file(params.target_regions_bed) )
        .set { ch_strelka_bed_plain }

    // Compress BED -> *.bed.gz
    BGZIP_BED(ch_strelka_bed_plain)
    BGZIP_BED.out.set { ch_strelka_bed_gz }   // path to .bed.gz

    // Index compressed BED with tabix (module expects (meta, bedgz))
    ch_strelka_bed_gz
        .map { bedgz -> tuple( [ id: 'strelka_bed' ], bedgz ) }
        .set { ch_strelka_bed_for_tabix }

    TABIX_TABIX(ch_strelka_bed_for_tabix)

    TABIX_TABIX.out.index
        .map { meta, idx -> idx }
        .set { ch_strelka_bed_index }         // path to .bed.gz.tbi



    /*
     * 13) Strelka2 germline (alternative variant caller)
     */

    // STRELKA_GERMLINE expects:
    // tuple(meta, bam, bai, target_bed, target_bed_index), fasta, fai
    ch_bam_bai
        .combine(ch_strelka_bed_gz)
        .combine(ch_strelka_bed_index)
        .map { meta, bam, bai, bedgz, bedidx ->
            tuple(meta, bam, bai, bedgz, bedidx)
        }
        .set { ch_strelka_input }

    STRELKA_GERMLINE(
        ch_strelka_input,
        ch_doc_ref,   // FASTA (path-only channel; matches module API)
        ch_doc_fai    // FASTA.fai (path-only channel)
    )

    // Capture Strelka outputs
    STRELKA_GERMLINE.out.vcf            .set { ch_strelka_variants_vcf }      // (meta, vcf)
    STRELKA_GERMLINE.out.vcf_tbi        .set { ch_strelka_variants_vcf_tbi }  // (meta, tbi)
    STRELKA_GERMLINE.out.genome_vcf     .set { ch_strelka_genome_vcf }
    STRELKA_GERMLINE.out.genome_vcf_tbi .set { ch_strelka_genome_vcf_tbi }



    /*
     * Strelka2 chrM / haploid filtering + PASS-only
     */

    // Pair variants VCF with its index -> (meta, vcf, tbi)
    ch_strelka_variants_vcf
        .join(ch_strelka_variants_vcf_tbi)
        .set { ch_strelka_variants_with_index }

    // Run local filter module
    STRELKA_CHRM_FILTER(ch_strelka_variants_with_index)

    // Capture filtered outputs
    STRELKA_CHRM_FILTER.out.vcf.set { ch_st2_vcf_raw }   // (meta, vcf)
    STRELKA_CHRM_FILTER.out.tbi.set { ch_st2_vcf_tbi }   // (meta, tbi)

    // Join filtered Strelka VCF + index -> (meta, vcf, tbi)
    def ch_st2_vcf = ch_st2_vcf_raw.join(ch_st2_vcf_tbi)



    /*
     * 14) bcftools concat DeepVariant + Strelka2
     */

    // Join DeepVariant and Strelka2 on meta.id
    def ch_concat_input = ch_dv_vcf
        .join(ch_st2_vcf)
        .map { meta, dv_vcf, dv_tbi, st_vcf, st_tbi ->
            def m = meta.clone()
            m.id = "${meta.id}.DV_ST2.comb"  // will appear in filenames

            def vcfs = [ dv_vcf, st_vcf ]
            def tbis = [ dv_tbi, st_tbi ]
            tuple(m, vcfs, tbis)             // matches: tuple val(meta), path(vcfs), path(tbi)
        }

    // bcftools concat; bcftools options are in nextflow.config (ext.args)
    BCFTOOLS_CONCAT_VCF(ch_concat_input)

    // Join concat VCF + index -> (meta, vcf, tbi)
    def ch_concat_vcf = BCFTOOLS_CONCAT_VCF.out.vcf
        .join(BCFTOOLS_CONCAT_VCF.out.tbi)

    // Build input for bcftools norm on combined calls
    def ch_norm_input_vcf = ch_concat_vcf
        .map { meta, vcf, tbi ->
            def m = meta.clone()
            m.id = "${meta.id}.DV_ST2.comb.norm"
            tuple(m, vcf, tbi)
        }



    /*
     * 15) bcftools norm on:
     *     a) combined DV+Strelka VCF
     *     b) DeepVariant-only VCF
     */

    // Reference FASTA for bcftools norm (tuple (meta_ref, fasta))
    def ch_ref_fasta = ref_val

    // (a) Run bcftools norm on combined DV+Strelka VCF
    BCFTOOLS_NORM_COMBINED(ch_norm_input_vcf, ch_ref_fasta)

    // (b) DeepVariant-only normalization
    def ch_dv_norm_input = ch_dv_vcf
        .map { meta, vcf, tbi ->
            def m = meta.clone()
            m.id = "${meta.id}.DV_only.norm"
            tuple(m, vcf, tbi)
        }

    BCFTOOLS_NORM_DVONLY(ch_dv_norm_input, ch_ref_fasta)

        
    /*
     * 16) Manta germline SV/CNV calling
     *     Uses the same deduplicated BAMs, reference, and target BED as Strelka.
     */

    // 16.1 First Manta input:
    // tuple(meta, input_bams, bai, target_bed, target_bed_tbi)
    //
    // Note: 'input' must be a LIST so the module can do input.collect{"--bam ${it}"}.
    def ch_manta_samples = ch_bam_bai
        .combine(ch_strelka_bed_gz)      // BED.GZ (value channel)
        .combine(ch_strelka_bed_index)   // BED.GZ.TBI (value channel)
        .map { meta, bam, bai, bedgz, bedidx ->
            tuple(meta, [bam], bai, bedgz, bedidx)   // [bam] => list with one BAM
        }

    // 16.2 Second input: reference FASTA as tuple(meta_ref, fasta)
    def ch_manta_fasta = ref_val        // already (meta_ref, fasta)

    // 16.3 Third input: FAI as tuple(meta_ref, fai)
    def ch_manta_fai = fai_val          // already (meta_ref, fai)

    // 16.4 Fourth input: Manta config file (single path)
    Channel.value( file(params.manta_config) )
        .set { ch_manta_config }

    // 16.5 Run Manta
    MANTA_GERMLINE(
        ch_manta_samples,   // (meta, [bam], bai, target_bed, target_bed_tbi)
        ch_manta_fasta,     // (meta_ref, fasta)
        ch_manta_fai,       // (meta_ref, fai)
        ch_manta_config     // path(config)
    )

    // 16.6 Capture outputs (optional channel names if you want them downstream)
    MANTA_GERMLINE.out.diploid_sv_vcf     .set { ch_manta_diploid_sv_vcf }
    MANTA_GERMLINE.out.diploid_sv_vcf_tbi .set { ch_manta_diploid_sv_vcf_tbi }
    MANTA_GERMLINE.out.candidate_sv_vcf   .set { ch_manta_candidate_sv_vcf }
    MANTA_GERMLINE.out.candidate_sv_vcf_tbi.set { ch_manta_candidate_sv_vcf_tbi }
    MANTA_GERMLINE.out.candidate_small_indels_vcf   .set { ch_manta_small_indels_vcf }
    MANTA_GERMLINE.out.candidate_small_indels_vcf_tbi.set { ch_manta_small_indels_vcf_tbi }





}