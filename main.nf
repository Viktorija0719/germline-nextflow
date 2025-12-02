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



/*
 * DOWNLOAD_REFERENCE
 */
process DOWNLOAD_REFERENCE {

    tag "${meta.id}"
    publishDir params.ref_dir, mode: 'copy'

    input:
    tuple val(meta), val(ref_url)

    output:
    tuple val(meta), path("${params.ref_name}.fasta")

    script:
    """
    set -euo pipefail

    mkdir -p "${params.ref_dir}"

    if [ -s "${params.ref_dir}/${params.ref_name}.fasta" ]; then
        echo ">>> Found existing FASTA in ${params.ref_dir}, reusing it."
        ln -s "${params.ref_dir}/${params.ref_name}.fasta" "${params.ref_name}.fasta"
    else
        echo ">>> Downloading reference from: ${ref_url}"
        wget -O "${params.ref_name}.fasta" "${ref_url}"
        cp "${params.ref_name}.fasta" "${params.ref_dir}/"
    fi
    """
}


/*
 * SAMTOOLS_FAIDX_SIMPLE
 * (simple faidx wrapper for the reference FASTA)
 */
process SAMTOOLS_FAIDX_SIMPLE {

    tag "${meta.id}"
    publishDir params.ref_dir, mode: 'copy'
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
    publishDir params.ref_dir, mode: 'copy'

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
 * Main workflow
 */
workflow {


        /*
     * 1) Reference
     */
    def meta_ref = [ id: params.ref_id ]

    Channel.of( [ meta_ref, params.ref_url ] )
        | DOWNLOAD_REFERENCE
        | set { ref_ch }   // [meta_ref, fasta]

    PICARD_CREATESEQUENCEDICTIONARY(ref_ch)
    BWA_INDEX(ref_ch)
    SAMTOOLS_FAIDX_SIMPLE(ref_ch)

    // Turn reference into VALUE channels
    ref_ch.first().set                    { ref_val }    // [meta_ref, fasta]
    BWA_INDEX.out.index.first().set       { index_val }  // [meta_ref, bwa_dir]
    SAMTOOLS_FAIDX_SIMPLE.out.first().set { fai_val }    // [meta_ref, fasta.fai]

    // For DepthOfCoverage: just the FAI path (no meta)
    fai_val
        .map { meta, fai -> fai }
        .set { ch_doc_fai }



    // For DepthOfCoverage: the .dict path from Picard
    PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
        .first()
        .map { meta, dict -> dict }
        .set { ch_doc_dict }
    


    /*
     * 1b) Qualimap regions:
     *     - If params.qualimap_bed is given -> use that file (e.g. WES targets)
     *     - Else -> auto-generate a whole-genome BED from the .fai
     */
    def ch_qualimap_regions

    if ( params.qualimap_bed ) {
        ch_qualimap_regions = Channel.value( file(params.qualimap_bed) )
    } else {
        MAKE_GENOME_BED(fai_val)
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
     * 12) DeepVariant (variant calling)
     *     Uses final deduplicated BAM + BAI, reference FASTA/FAI, target regions BED,
     *     dummy GZI and dummy PAR BED.
     */

    // 12.1 Target regions: use params.target_regions_bed by default
    Channel.fromPath(params.target_regions_bed, checkIfExists: true)
           .set { ch_target_regions_bed }   // value channel with a single BED file

    // 12.2 Combine BAM+BAI with target intervals -> (meta, bam, bai, intervals)
    ch_bam_bai
        .combine(ch_target_regions_bed)
        .map { meta, bam, bai, intervals ->
            tuple(meta, bam, bai, intervals)
        }
        .set { ch_deep_bam }

    // 12.3 Reference FASTA and FAI as-is
    def ch_deep_fasta = ref_val   // (meta_ref, fasta)
    def ch_deep_fai   = fai_val   // (meta_ref, fasta.fai)

    // 12.4 Dummy GZI: must be a *different* filename than the .fai
    Channel.value(
        tuple( [ id: params.ref_id ], file(params.deepvariant_gzi_dummy) )
    ).set { ch_deep_gzi }   // (meta_ref, dummy_fasta.gzi)

    // 12.5 Dummy PAR BED (empty file, no effect)
    Channel.value(
        tuple( [ id: params.ref_id ], file(params.par_regions_bed) )
    ).set { ch_deep_par }   // (meta_ref, empty_par_regions.bed)

    // 12.6 Run DeepVariant with all 5 required input tuples
    DEEPVARIANT_RUNDEEPVARIANT(
        ch_deep_bam,    // (meta, bam, bai, intervals)
        ch_deep_fasta,  // (meta_ref, fasta)
        ch_deep_fai,    // (meta_ref, fasta.fai)
        ch_deep_gzi,    // (meta_ref, dummy_fasta.gzi)
        ch_deep_par     // (meta_ref, empty_par_regions.bed)
    )

    // 12.7 Capture outputs
    DEEPVARIANT_RUNDEEPVARIANT.out.vcf        .set { ch_deep_vcf }
    DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index  .set { ch_deep_vcf_index }
    DEEPVARIANT_RUNDEEPVARIANT.out.gvcf       .set { ch_deep_gvcf }
    DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index .set { ch_deep_gvcf_index }




}
