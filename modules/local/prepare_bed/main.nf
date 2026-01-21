/*
 * modules/local/prepare_bed/main.nf
 *
 * Purpose:
 *   - Clean BED-like intervals (drop comments/track, require >=3 cols)
 *   - Sort
 *   - Validate contigs vs reference FASTA .fai
 *   - Optional chr-prefix auto-fix (add/remove "chr")
 *
 * Notes:
 *   - Designed to be called multiple times in main.nf. Always capture the returned handle:
 *       def pb = PREPARE_BED(ch_some_bed, ch_fai)
 *       def bed_prepared = pb.out.out_bed
 *
 * Params used:
 *   params.bed_validate_strict  (default true)
 *   params.bed_contig_autofix   (default true)
 *   params.bed_cache_dir        (optional storeDir base)
 */

process PREPARE_BED {
    tag "${bed.baseName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"

    /*
     * Optional HPC cache of prepared BEDs:
     * - Only effective if bed_cache_dir is set.
     * - Keeps outputs deduplicated across runs on shared filesystem.
     */
    storeDir { params.bed_cache_dir ? "${params.bed_cache_dir}/prepared/${bed.baseName}" : null }

    input:
    path bed
    path fasta_fai

    output:
    path "${bed.baseName}.prepared.bed", emit: out_bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def strict  = (params.bed_validate_strict == null ? true : params.bed_validate_strict) as boolean
    def autofix = (params.bed_contig_autofix  == null ? true : params.bed_contig_autofix)  as boolean

    """
    set -euo pipefail

    in="${bed}"
    fai="${fasta_fai}"

    # 1) Clean + keep >=3 cols, keep extra cols intact, then sort
    awk 'BEGIN{FS=OFS="\\t"}
         \$0 ~ /^#/     {next}
         \$0 ~ /^track/ {next}
         NF < 3         {next}
         {print \$0}' "\$in" | bedtools sort -i - > cleaned.sorted.bed

    if [ ! -s cleaned.sorted.bed ]; then
      echo "ERROR: BED is empty after cleaning: \$in" >&2
      exit 1
    fi

    # 2) Contig validation
    cut -f1 "\$fai" | sort -u > fai.contigs
    awk 'BEGIN{FS="\\t"}{print \$1}' cleaned.sorted.bed | sort -u > bed.contigs

    comm -12 fai.contigs bed.contigs > overlap.contigs || true

    if [ ! -s overlap.contigs ]; then
      if ${autofix}; then
        # Try add/remove "chr" prefix
        if grep -q '^chr' bed.contigs && ! grep -q '^chr' fai.contigs; then
          awk 'BEGIN{FS=OFS="\\t"}{sub(/^chr/,"",\$1); print \$0}' cleaned.sorted.bed | bedtools sort -i - > fixed.sorted.bed
        elif ! grep -q '^chr' bed.contigs && grep -q '^chr' fai.contigs; then
          awk 'BEGIN{FS=OFS="\\t"}{\$1=(\$1 ~ /^chr/ ? \$1 : "chr"\$1); print \$0}' cleaned.sorted.bed | bedtools sort -i - > fixed.sorted.bed
        else
          cp cleaned.sorted.bed fixed.sorted.bed
        fi

        awk 'BEGIN{FS="\\t"}{print \$1}' fixed.sorted.bed | sort -u > bed.fixed.contigs
        comm -12 fai.contigs bed.fixed.contigs > overlap.fixed.contigs || true

        if [ -s overlap.fixed.contigs ]; then
          mv fixed.sorted.bed cleaned.sorted.bed
        fi
      fi
    fi

    # 3) Re-check overlap after optional fix
    awk 'BEGIN{FS="\\t"}{print \$1}' cleaned.sorted.bed | sort -u > bed.final.contigs
    comm -12 fai.contigs bed.final.contigs > overlap.final.contigs || true

    if [ ! -s overlap.final.contigs ]; then
      echo "ERROR: BED contigs do not match reference contigs from FASTA .fai" >&2
      echo "BED example contigs:" >&2
      head -n 10 bed.final.contigs >&2 || true
      echo "FAI example contigs:" >&2
      head -n 10 fai.contigs >&2 || true
      if ${strict}; then
        exit 1
      else
        echo "WARNING: continuing despite contig mismatch (bed_validate_strict=false)" >&2
      fi
    fi

    mv cleaned.sorted.bed ${bed.baseName}.prepared.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
