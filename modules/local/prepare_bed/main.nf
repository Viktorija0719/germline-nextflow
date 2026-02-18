/*
 * modules/local/prepare_bed/main.nf
 * Portable, Groovy-safe BED preparation
 */

process PREPARE_BED {
    tag "${bed.baseName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"

    // Optional cache of prepared BEDs
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
    def strict = (params.containsKey('bed_validate_strict')
        ? params.bed_validate_strict
        : true) as boolean

    // backward-compatible: new param preferred, old param accepted
    def autofix = (params.containsKey('bed_contig_autofix')
        ? params.bed_contig_autofix
        : (params.containsKey('bed_autocorrect_chr_prefix')
            ? params.bed_autocorrect_chr_prefix
            : true)) as boolean

    """
    set -euo pipefail

    in="${bed}"
    fai="${fasta_fai}"

    # 1) Clean input
    # - remove CRLF safely with tr
    # - drop comments/track lines
    # - require >=3 cols and valid numeric coordinates
    tr -d '\\r' < "\$in" | \
    awk 'BEGIN{FS=OFS="\\t"}
         \$0 ~ /^#/     {next}
         \$0 ~ /^track/ {next}
         NF < 3                               {next}
         \$2 !~ /^[0-9]+$/ || \$3 !~ /^[0-9]+$/ {next}
         \$2 >= \$3                            {next}
         {print \$0}' > cleaned.raw.bed

    if [ ! -s cleaned.raw.bed ]; then
      echo "ERROR: BED is empty after cleaning: \$in" >&2
      exit 1
    fi

    # 2) Portable sorting (no bedtools sort, no sort -V)
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n cleaned.raw.bed > cleaned.sorted.bed

    # 3) Contig validation
    cut -f1 "\$fai" | LC_ALL=C sort -u > fai.contigs
    awk 'BEGIN{FS="\\t"}{print \$1}' cleaned.sorted.bed | LC_ALL=C sort -u > bed.contigs

    LC_ALL=C comm -12 fai.contigs bed.contigs > overlap.contigs || true

    if [ ! -s overlap.contigs ]; then
      if ${autofix}; then
        # Try add/remove chr prefix
        if grep -q '^chr' bed.contigs && ! grep -q '^chr' fai.contigs; then
          awk 'BEGIN{FS=OFS="\\t"}{sub(/^chr/,"",\$1); print \$0}' cleaned.sorted.bed > fixed.raw.bed
        elif ! grep -q '^chr' bed.contigs && grep -q '^chr' fai.contigs; then
          awk 'BEGIN{FS=OFS="\\t"}{\$1=(\$1 ~ /^chr/ ? \$1 : "chr"\$1); print \$0}' cleaned.sorted.bed > fixed.raw.bed
        else
          cp cleaned.sorted.bed fixed.raw.bed
        fi

        LC_ALL=C sort -k1,1 -k2,2n -k3,3n fixed.raw.bed > fixed.sorted.bed

        awk 'BEGIN{FS="\\t"}{print \$1}' fixed.sorted.bed | LC_ALL=C sort -u > bed.fixed.contigs
        LC_ALL=C comm -12 fai.contigs bed.fixed.contigs > overlap.fixed.contigs || true

        if [ -s overlap.fixed.contigs ]; then
          mv fixed.sorted.bed cleaned.sorted.bed
        fi
      fi
    fi

    # 4) Final overlap check
    awk 'BEGIN{FS="\\t"}{print \$1}' cleaned.sorted.bed | LC_ALL=C sort -u > bed.final.contigs
    LC_ALL=C comm -12 fai.contigs bed.final.contigs > overlap.final.contigs || true

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

    btv=\$(bedtools --version 2>/dev/null | sed 's/bedtools v//' || true)
    [ -n "\$btv" ] || btv="unknown"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$btv
    END_VERSIONS
    """
}
