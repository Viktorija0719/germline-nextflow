#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Include nf-core modules
 */
include { BWA_INDEX                        } from './modules/nf-core/bwa/index/main'
include { BWA_MEM                          } from './modules/nf-core/bwa/mem/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from './modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_MARKDUPLICATES           } from './modules/nf-core/picard/markduplicates/main'


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
    def sample   = meta.sample ?: meta.id
    def lib      = meta.library ?: meta.patient ?: meta.id
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
 *
 * Take Picard MarkDuplicates metrics and wrap them into a
 * Picard-like header + body file, as you described.
 */
process PICARDLIKE_DUPMETRICS {

    tag "${meta.id}"

    input:
    tuple val(meta), path(metrics)

    // one output per metrics file
    output:
    tuple val(meta), path("*.duplicate_metrics_picard_like.txt")

    script:
    """
    set -euo pipefail

    in_metrics="${metrics}"
    prefix="${meta.id}"

    # Drop the first 3 lines and keep the tabular metrics
    tail -n +4 "\${in_metrics}" > tmp.txt

    # Build a Picard-like header; adjust text as you like
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
    BWA_INDEX.out.index.first().set       { index_val }  // [meta_ref, index_dir]
    SAMTOOLS_FAIDX_SIMPLE.out.first().set { fai_val }    // [meta_ref, fasta.fai]


    /*
     * 2) Samples from samplesheet.csv
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
     * 3) BWA-MEM alignment (per sample)
     */
    BWA_MEM(
        ch_reads,
        index_val,  // value
        ref_val,    // value
        true
    )

    BWA_MEM.out.bam.set { ch_bam_raw }   // [meta_sample, bam]


    /*
     * 4) Add read groups (per sample)
     */
    ADD_READGROUPS(ch_bam_raw)
    ADD_READGROUPS.out.set { ch_bam_rg }   // [meta_sample, bam_with_rg]


    /*
     * 5) Mark duplicates (per sample)
     */
    ch_bam_rg
        .map { meta, bam ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}.dedup"
            tuple(new_meta, bam)
        }
        .set { ch_bam_for_markdup }

    PICARD_MARKDUPLICATES(
        ch_bam_for_markdup,
        ref_val,   // value
        fai_val    // value
    )

    // Capture metrics output channel
    PICARD_MARKDUPLICATES.out.metrics.set { ch_dup_metrics }

    /*
     * 6) Picard-like duplication metrics (per sample)
     */
    PICARDLIKE_DUPMETRICS(ch_dup_metrics)
}



// workflow {

//     // meta info for this reference
//     def meta_ref = [ id: params.ref_id ]

//     /*
//      * 1) Download or reuse the FASTA
//      */
//     Channel.of( [ meta_ref, params.ref_url ] )
//         | DOWNLOAD_REFERENCE
//         | set { ref_ch }        // tuple: [meta_ref, fasta]


//     /*
//      * 2) Reference indexes
//      *
//      * SAMTOOLS_FAIDX needs:
//      *   - fasta channel        : ref_ch
//      *   - existing .fai file   : empty (none yet)
//      *   - get_sizes flag       : false (no *.sizes file)
//      */
//     def fai_existing_ch = Channel.empty()
//     def get_sizes_ch    = Channel.of(false)

//     SAMTOOLS_FAIDX(ref_ch, fai_existing_ch, get_sizes_ch)

//     // Keep named channels for reference FASTA + BWA index
//     def ch_fasta = ref_ch

//     BWA_INDEX(ch_fasta)
//     def ch_index = BWA_INDEX.out.index

//     // Picard sequence dictionary
//     PICARD_CREATESEQUENCEDICTIONARY(ch_fasta)


//     /*
//      * 3) Build a channel with reads from the samplesheet
//      *
//      * samplesheet.csv:
//      *   patient,sample,lane,fastq_1,fastq_2
//      *
//      * Channel structure:
//      *   [ meta_sample, [fastq1, fastq2] ]
//      */
//     Channel
//         .fromPath(params.input, checkIfExists: true)
//         .splitCsv(header: true)
//         .map { row ->
//             def meta_sample = [
//                 id         : row.sample,
//                 sample     : row.sample,
//                 patient    : row.patient,
//                 lane       : row.lane,
//                 single_end : false
//             ]
//             tuple(meta_sample, [ file(row.fastq_1), file(row.fastq_2) ])
//         }
//         .set { ch_samples }


//     /*
//      * 4) Alignment with BWA_MEM (nf-core module)
//      *
//      * BWA_MEM(
//      *   ch_samples : [meta_sample, reads]
//      *   ch_index   : [meta_ref, bwa_index_dir]
//      *   ch_fasta   : [meta_ref, fasta]
//      *   true       : sort_bam -> samtools sort in-module
//      * )
//      */
//     BWA_MEM(
//         ch_samples,
//         ch_index,
//         ch_fasta,
//         true       // sort_bam
//     )

//     // BAMs from BWA_MEM (one per sample)
//     def ch_bam = BWA_MEM.out.bam


//     /*
//      * 5) Mark duplicates with Picard MarkDuplicates (nf-core module)
//      *
//      * Input:  [meta_sample, bam]
//      * Output: bam (dedup), metrics, optional bai
//      * (CREATE_INDEX=true is set in nextflow.config)
//      */
//     PICARD_MARKDUPLICATES(ch_bam, ref_ch, SAMTOOLS_FAIDX.out.fai)

//     def ch_dedup_bam     = PICARD_MARKDUPLICATES.out.bam
//     def ch_dedup_metrics = PICARD_MARKDUPLICATES.out.metrics

//     // At this point you have:
//     //  - reference prepared (FASTA, FAI, dict, BWA index)
//     //  - per-sample deduplicated BAMs + metrics
//     // Next steps: add variant calling modules (bcftools / GATK etc.)
// }


// workflow {

//     // meta info for this reference
//     def meta = [ id: params.ref_id ]

//     /*
//      * 1) Download or reuse the FASTA
//      */
//     Channel.of( [ meta, params.ref_url ] )
//         | DOWNLOAD_REFERENCE
//         | set { ref_ch }        // tuple: [meta, fasta]

//     /*
//      * 2) Prepare channels required by SAMTOOLS_FAIDX
//      *
//      *    - fa_ch       : our ref_ch (meta + fasta)
//      *    - fai_ch      : empty for now (no pre-existing .fai to reuse)
//      *    - get_sizes_ch: boolean; false = don't emit *.sizes file
//      */
//     def fai_existing_ch = Channel.empty()
//     def get_sizes_ch    = Channel.of(false)

//     // This now matches the nf-core module signature: 3 input channels
//     SAMTOOLS_FAIDX(ref_ch, fai_existing_ch, get_sizes_ch)

//     /*
//      * 3) Other indexing modules: these only expect a single channel
//      *    containing [meta, fasta]
//      */
//     PICARD_CREATESEQUENCEDICTIONARY(ref_ch)
//     BWA_INDEX(ref_ch)
// }
