process XHMM_PIPELINE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/xhmm:0.0.0.2016_01_04.cc14e52--hedee03e_3' :
        'quay.io/biocontainers/xhmm:0.0.0.2016_01_04.cc14e52--hedee03e_3' }"

    input:
    tuple val(meta), path(gatk_summaries)     // list of *.sample_interval_summary
    path annotated_intervals                  // may be empty placeholder file
    path extreme_gc_targets                   // may be empty placeholder file
    path xhmm_params

    output:
    // Use closures so meta/task.ext are available (prevents "null.*" output expectations)
    tuple val(meta), path { "${task.ext.prefix ?: meta.id ?: 'DATA'}.xcnv" }     , emit: xcnv
    tuple val(meta), path { "${task.ext.prefix ?: meta.id ?: 'DATA'}.aux_xcnv" } , emit: aux_xcnv
    tuple val(meta), path { "${task.ext.prefix ?: meta.id ?: 'DATA'}*" }         , emit: all
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Output prefix used inside the script
    def prefix = task.ext.prefix ?: (meta.id ?: 'DATA')

    // GC thresholds
    def gc_low  = task.ext.gc_low  ?: '0.1'
    def gc_high = task.ext.gc_high ?: '0.9'

    // XHMM args (optional)
    def merge_args      = (task.ext.merge_args      ?: '').trim()
    def matrix1_args    = (task.ext.matrix1_args    ?: '').trim()
    def normalize_args  = (task.ext.normalize_args  ?: '').trim()
    def matrix2_args    = (task.ext.matrix2_args    ?: '').trim()
    def discover_args   = (task.ext.discover_args   ?: '').trim()

    // Filter thresholds (make sane defaults; override via task.ext.* from nextflow.config)
    def minTargetSize    = task.ext.minTargetSize    ?: 10
    def maxTargetSize    = task.ext.maxTargetSize    ?: 10000
    def minMeanTargetRD  = task.ext.minMeanTargetRD  ?: 10
    def maxMeanTargetRD  = task.ext.maxMeanTargetRD  ?: 1000000
    def minMeanSampleRD  = task.ext.minMeanSampleRD  ?: 0
    def maxMeanSampleRD  = task.ext.maxMeanSampleRD  ?: 5000
    def maxSdSampleRD    = task.ext.maxSdSampleRD    ?: 1000000
    def maxSdTargetRD    = task.ext.maxSdTargetRD    ?: 30

    // Guardrails
    def min_samples_hard = task.ext.min_samples_hard ?: 2   // PCA cannot work with <2
    def min_samples_soft = task.ext.min_samples_soft ?: 10  // XHMM suggestion

    """
    set -euo pipefail

    echo "[XHMM] GATK summaries count: ${gatk_summaries.size()}"
    echo "[XHMM] prefix: ${prefix}"

    # ------------------------------------------------------------
    # 0) Ensure tab-delimited summaries (XHMM fails on CSV headers)
    # ------------------------------------------------------------
    mkdir -p summaries_tsv
    for f in ${gatk_summaries.join(' ')}; do
      bn=\$(basename "\$f")
      out="summaries_tsv/\$bn"
      if head -n 1 "\$f" | grep -q "," ; then
        tr ',' '\\t' < "\$f" > "\$out"
      else
        cp "\$f" "\$out"
      fi
    done
    ls -1 summaries_tsv/* > gatkdepths.list

    # ------------------------------------------------------------------
    # 1) mergeGATKdepths -> read-depth matrix (samples x targets)
    # ------------------------------------------------------------------
    xhmm --mergeGATKdepths \\
      -o ${prefix}.RD.txt \\
      --GATKdepthsList gatkdepths.list \\
      ${merge_args}

    # ------------------------------------------------------------------
    # 2) Extreme GC targets list
    # ------------------------------------------------------------------
    if [[ -s "${extreme_gc_targets}" ]]; then
        echo "[XHMM] Using provided extreme GC targets: ${extreme_gc_targets}"
        cp "${extreme_gc_targets}" ${prefix}.extreme_gc_targets.txt
        : > ${prefix}.locus_GC.txt || true
    elif [[ -s "${annotated_intervals}" ]]; then
        echo "[XHMM] Deriving extreme GC targets from: ${annotated_intervals}"

        awk -v OFS="\\t" '
          \$0 ~ /^@/ { next }
          !hdr {
            for (i=1;i<=NF;i++) {
              if (\$i=="CONTIG" || \$i=="contig") c=i
              if (\$i=="START"  || \$i=="start")  s=i
              if (\$i=="END"    || \$i=="end")    e=i
              if (\$i=="GC_CONTENT" || \$i=="GCContent" || \$i ~ /^GC/) gc=i
            }
            if (!c || !s || !e || !gc) { print "ERROR: missing CONTIG/START/END/GC column in header" > "/dev/stderr"; exit 1 }
            hdr=1
            next
          }
          { print \$c ":" \$s "-" \$e, \$gc }
        ' "${annotated_intervals}" > ${prefix}.locus_GC.txt

        awk -v lo=${gc_low} -v hi=${gc_high} '\$2 < lo || \$2 > hi {print \$1}' \\
          ${prefix}.locus_GC.txt > ${prefix}.extreme_gc_targets.txt
    else
        echo "[XHMM] No GC inputs provided; using empty exclude list"
        : > ${prefix}.extreme_gc_targets.txt
        : > ${prefix}.locus_GC.txt
    fi

    # ------------------------------------------------------------------
    # 3) matrix: filter + mean-center by TARGET
    # ------------------------------------------------------------------
    xhmm --matrix \\
      -r ${prefix}.RD.txt \\
      --centerData --centerType target \\
      -o ${prefix}.filtered_centered.RD.txt \\
      --outputExcludedTargets ${prefix}.filtered_centered.RD.txt.filtered_targets.txt \\
      --outputExcludedSamples ${prefix}.filtered_centered.RD.txt.filtered_samples.txt \\
      --excludeTargets ${prefix}.extreme_gc_targets.txt \\
      --minTargetSize ${minTargetSize} --maxTargetSize ${maxTargetSize} \\
      --minMeanTargetRD ${minMeanTargetRD} --maxMeanTargetRD ${maxMeanTargetRD} \\
      --minMeanSampleRD ${minMeanSampleRD} --maxMeanSampleRD ${maxMeanSampleRD} \\
      --maxSdSampleRD ${maxSdSampleRD} \\
      ${matrix1_args}

    kept=\$(awk 'NR>1{c++}END{print c+0}' ${prefix}.filtered_centered.RD.txt)
    echo "[XHMM] Samples kept after matrix1: \$kept"
    if [[ -s ${prefix}.filtered_centered.RD.txt.filtered_samples.txt ]]; then
      echo "[XHMM] Excluded samples (matrix1):"
      cat ${prefix}.filtered_centered.RD.txt.filtered_samples.txt || true
    fi

    if [[ "\$kept" -lt ${min_samples_hard} ]]; then
      echo "[XHMM][ERROR] Too few samples (\$kept) left after filtering for PCA normalization. Relax thresholds (e.g. maxMeanSampleRD/maxSdSampleRD) or use more samples." >&2
      exit 1
    fi

    if [[ "\$kept" -lt ${min_samples_soft} ]]; then
      echo "[XHMM][WARN] Only \$kept samples. XHMM is typically run with ~10+ samples; results may be unstable." >&2
    fi

    # ------------------------------------------------------------------
    # 4) PCA
    # ------------------------------------------------------------------
    xhmm --PCA \\
      -r ${prefix}.filtered_centered.RD.txt \\
      --PCAfiles ${prefix}.RD_PCA

    # ------------------------------------------------------------------
    # 5) normalize (PCA-based)
    # ------------------------------------------------------------------
    xhmm --normalize \\
      -r ${prefix}.filtered_centered.RD.txt \\
      --PCAfiles ${prefix}.RD_PCA \\
      --normalizeOutput ${prefix}.PCA_normalized.txt \\
      --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7 \\
      ${normalize_args}

    # ------------------------------------------------------------------
    # 6) matrix: filter + z-score by SAMPLE
    # ------------------------------------------------------------------
    xhmm --matrix \\
      -r ${prefix}.PCA_normalized.txt \\
      --centerData --centerType sample --zScoreData \\
      -o ${prefix}.PCA_normalized.filtered.sample_zscores.RD.txt \\
      --outputExcludedTargets ${prefix}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \\
      --outputExcludedSamples ${prefix}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \\
      --maxSdTargetRD ${maxSdTargetRD} \\
      ${matrix2_args}

    # ------------------------------------------------------------------
    # 7) matrix: filter original RD to match the final filtered set
    # ------------------------------------------------------------------
    xhmm --matrix \\
      -r ${prefix}.RD.txt \\
      --excludeTargets ${prefix}.filtered_centered.RD.txt.filtered_targets.txt \\
      --excludeTargets ${prefix}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \\
      --excludeSamples ${prefix}.filtered_centered.RD.txt.filtered_samples.txt \\
      --excludeSamples ${prefix}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \\
      -o ${prefix}.same_filtered.RD.txt

    # ------------------------------------------------------------------
    # 8) discover CNVs
    # ------------------------------------------------------------------
    xhmm --discover \\
      -p "${xhmm_params}" \\
      -r ${prefix}.PCA_normalized.filtered.sample_zscores.RD.txt \\
      -R ${prefix}.same_filtered.RD.txt \\
      -c ${prefix}.xcnv \\
      -a ${prefix}.aux_xcnv \\
      -s ${prefix} \\
      ${discover_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xhmm: "\$(xhmm 2>&1 | head -n 1 || echo unknown)"
    END_VERSIONS
    """

    stub:
    def p = task.ext.prefix ?: (meta.id ?: 'DATA')
    """
    touch ${p}.RD.txt
    touch ${p}.filtered_centered.RD.txt
    touch ${p}.filtered_centered.RD.txt.filtered_targets.txt
    touch ${p}.filtered_centered.RD.txt.filtered_samples.txt
    touch ${p}.PCA_normalized.txt
    touch ${p}.PCA_normalized.filtered.sample_zscores.RD.txt
    touch ${p}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt
    touch ${p}.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt
    touch ${p}.same_filtered.RD.txt
    touch ${p}.xcnv
    touch ${p}.aux_xcnv
    touch ${p}.locus_GC.txt
    touch ${p}.extreme_gc_targets.txt
    touch versions.yml
    """
}
