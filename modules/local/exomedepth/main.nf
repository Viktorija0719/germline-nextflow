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
    path "${prefix}/calls/*.exomedepth.cnv.csv", emit: calls
    path "${prefix}/exomedepth_summary.tsv",        emit: summary
    path "${prefix}/reference_sets.tsv",           emit: refsets
    path "${prefix}/inferred_sex_from_idxstats.tsv", emit: sex
    path "${prefix}/run_config.tsv",               emit: run_config
    path "${prefix}/bam_manifest.tsv",             emit: bam_manifest

  when:
    task.ext.when == null || task.ext.when

  script:
    prefix     = task.ext.prefix ?: 'exomedepth'
    extra_args = task.ext.args   ?: ''
    include_x  = (task.ext.include_x != null ? task.ext.include_x : false)
    bed_coord  = task.ext.bed_coord ?: 'bed0'
    include_chr = task.ext.include_chr ?: 'auto'   // strongly consider setting TRUE/FALSE explicitly

    // build manifests for robust linking
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
      ln -sf "\$bam" "bams/\${sample}.bam"
      ln -sf "\$bai" "bams/\${sample}.bam.bai"
    done < bam_links.tsv

    # idxstats (optional, but needed for include_x TRUE to match your script's design)
    if [[ "${idxstats ? '1' : '0'}" == "1" ]]; then
      cat > idx_links.tsv <<'EOF'
sample\tidx
${idx_lines}
EOF
      while IFS=\$'\\t' read -r sample idx; do
        [[ "\$sample" == "sample" ]] && continue
        ln -sf "\$idx" "idxstats/\${sample}.idxstats"
      done < idx_links.tsv
      IDX_ARG="--idxstats_dir idxstats"
    else
      IDX_ARG=""
    fi

    # minimal samplesheet with required 'sample' column
    {
      echo "sample"
      ${samples.collect{ "echo '${it}'" }.join('\n')}
    } > samplesheet_exomedepth.csv

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
}
