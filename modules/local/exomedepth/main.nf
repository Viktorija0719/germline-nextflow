process EXOMEDEPTH {
  label 'process_medium'
  tag "${meta.id}"
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/r-exomedepth:1.1.18--r44hb2a3317_0' :
      'quay.io/biocontainers/r-exomedepth:1.1.18--r44hb2a3317_0' }"

  input:
    tuple val(meta), val(samples), path(bams), path(bais), path(idxstats), path(fasta), path(bed)

  output:
    path "${prefix}/calls/*.exomedepth.cnv.csv",      emit: calls
    path "${prefix}/exomedepth_summary.tsv",          emit: summary
    path "${prefix}/reference_sets.tsv",              emit: refsets
    path "${prefix}/inferred_sex_from_idxstats.tsv",  emit: sex
    path "${prefix}/run_config.tsv",                  emit: run_config
    path "${prefix}/bam_manifest.tsv",                emit: bam_manifest

  when:
    task.ext.when == null || task.ext.when

  script:
    prefix      = task.ext.prefix ?: 'exomedepth'
    extra_args  = task.ext.args   ?: ''
    include_x   = (task.ext.include_x != null ? task.ext.include_x : false)
    bed_coord   = task.ext.bed_coord ?: 'bed0'
    include_chr = task.ext.include_chr ?: 'auto'

    // Build manifests for robust linking
    def bam_lines = (0..<samples.size()).collect { i ->
      "${samples[i]}\t${bams[i]}\t${bais[i]}"
    }.join('\n')

    def idx_lines = idxstats ? (0..<samples.size()).collect { i ->
      "${samples[i]}\t${idxstats[i]}"
    }.join('\n') : ""

    """
    set -euo pipefail

    mkdir -p bams idxstats ${prefix}

    # BAM/BAI symlinks as <sample>.bam + <sample>.bam.bai
    cat > bam_links.tsv <<'EOF'
sample\tbam\tbai
${bam_lines}
EOF

    while IFS=\$'\\t' read -r sample bam bai; do
      [[ "\$sample" == "sample" ]] && continue

      # Resolve staged names to absolute source paths
      bam_src="\$bam"; [[ "\$bam_src" = /* ]] || bam_src="\$PWD/\$bam_src"
      bai_src="\$bai"; [[ "\$bai_src" = /* ]] || bai_src="\$PWD/\$bai_src"

      [[ -e "\$bam_src" ]] || { echo "Missing staged BAM source: \$bam_src" >&2; exit 1; }
      [[ -e "\$bai_src" ]] || { echo "Missing staged BAI source: \$bai_src" >&2; exit 1; }

      ln -sfn "\$bam_src" "bams/\${sample}.bam"
      ln -sfn "\$bai_src" "bams/\${sample}.bam.bai"
    done < bam_links.tsv

    # idxstats optional
    if [[ "${idxstats ? '1' : '0'}" == "1" ]]; then
      cat > idx_links.tsv <<'EOF'
sample\tidx
${idx_lines}
EOF
      while IFS=\$'\\t' read -r sample idx; do
        [[ "\$sample" == "sample" ]] && continue

        idx_src="\$idx"; [[ "\$idx_src" = /* ]] || idx_src="\$PWD/\$idx_src"
        [[ -e "\$idx_src" ]] || { echo "Missing staged idxstats source: \$idx_src" >&2; exit 1; }

        ln -sfn "\$idx_src" "idxstats/\${sample}.idxstats"
      done < idx_links.tsv
      IDX_ARG="--idxstats_dir idxstats"
    else
      IDX_ARG=""
    fi

    # Minimal samplesheet with required 'sample' column
    {
      echo "sample"
      ${samples.collect{ "echo '${it}'" }.join('\n')}
    } > samplesheet_exomedepth.csv

    # Quick link sanity
    awk 'NR>1{print \$1}' bam_links.tsv | while read -r s; do
      [[ -e "bams/\${s}.bam" ]]     || { echo "Broken link: bams/\${s}.bam" >&2; exit 1; }
      [[ -e "bams/\${s}.bam.bai" ]] || { echo "Broken link: bams/\${s}.bam.bai" >&2; exit 1; }
    done

    Rscript ${projectDir}/modules/local/exomedepth/run_exomedepth.R \\
      --fasta ${fasta} \\
      --bed ${bed} \\
      --bed_coord ${bed_coord} \\
      --bam_dir bams \\
      --samplesheet samplesheet_exomedepth.csv \\
      --outdir ${prefix} \\
      --include_chr ${include_chr} \\
      --include_x ${include_x} \\
      \$IDX_ARG \\
      ${extra_args}
    """