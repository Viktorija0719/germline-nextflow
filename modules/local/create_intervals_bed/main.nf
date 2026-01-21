process CREATE_INTERVALS_BED {
    tag "$intervals"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'quay.io/biocontainers/gawk:5.3.1' }"

    input:
    path(intervals)
    path(fasta_fai)

    output:
    path("*.bed")       , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def nps = (params.nucleotides_per_second ?: 10000000) as double
    def strict = (params.bed_validate_strict ?: true).toString()
    def autofix = (params.bed_contig_autofix ?: true).toString()
    def chunk_target = (params.chunk_target_seconds ?: 600) as double
    def longest_factor = (params.chunk_longest_factor ?: 1.05) as double

    """
    set -euo pipefail

    IN="${intervals}"
    FAI="${fasta_fai}"

    MASTER="master.intervals"
    CLEAN="clean.intervals"
    PREP="prepared.intervals"

    # 1) Convert input into BED-like (tab, >=3 cols) master intervals file
    shopt -s nocasematch
    if [[ "\$IN" == *.bed ]]; then
        cat "\$IN" > "\$MASTER"
    elif [[ "\$IN" == *.interval_list ]]; then
        # Picard interval_list: chr  start  end  ... (start is 1-based)
        awk 'BEGIN{FS=OFS="\\t"} \$0 ~ /^@/ {next} NF>=3 {print \$1, \$2-1, \$3}' "\$IN" > "\$MASTER"
    else
        # region strings: chr:start-end
        awk -vFS="[:-]" -vOFS="\\t" 'NF>=3 {print \$1, \$2-1, \$3}' "\$IN" > "\$MASTER"
    fi
    shopt -u nocasematch

    # 2) Clean (drop comments/track, require >=3 cols)
    awk 'BEGIN{FS=OFS="\\t"}
         \$0 ~ /^#/     {next}
         \$0 ~ /^track/ {next}
         NF < 3         {next}
         {print \$0}' "\$MASTER" > "\$CLEAN"

    if [ ! -s "\$CLEAN" ]; then
      echo "ERROR: intervals contain no usable rows after cleaning" >&2
      exit 1
    fi

    # 3) Validate contigs against .fai (no external cut/comm needed)
    #    First pass: check if ANY contig matches
    awk -v FS="\\t" -v STRICT="${strict}" '
      FNR==NR { fai[\$1]=1; if(\$1 ~ /^chr/) fai_chr=1; next }
      { if(\$1 in fai) { ok=1 } if(\$1 ~ /^chr/) bed_chr=1 }
      END { if(ok) exit 0; exit 2 }
    ' "\$FAI" "\$CLEAN"
    status=\$?

    if [ "\$status" -eq 2 ]; then
      if [ "${autofix}" = "true" ]; then
        # detect chr style
        fai_chr=\$(awk 'NR==1{print (\$1 ~ /^chr/)?1:0; exit}' "\$FAI")
        bed_chr=\$(awk 'NF>=3{print (\$1 ~ /^chr/)?1:0; exit}' "\$CLEAN")

        if [ "\$bed_chr" -eq 1 ] && [ "\$fai_chr" -eq 0 ]; then
          awk 'BEGIN{FS=OFS="\\t"}{sub(/^chr/,"",\$1); print}' "\$CLEAN" > "\$PREP"
        elif [ "\$bed_chr" -eq 0 ] && [ "\$fai_chr" -eq 1 ]; then
          awk 'BEGIN{FS=OFS="\\t"}{\$1=(\$1 ~ /^chr/ ? \$1 : "chr"\$1); print}' "\$CLEAN" > "\$PREP"
        else
          cp "\$CLEAN" "\$PREP"
        fi

        # re-check overlap after fix
        awk -v FS="\\t" '
          FNR==NR { fai[\$1]=1; next }
          { if(\$1 in fai) ok=1 }
          END { if(ok) exit 0; exit 2 }
        ' "\$FAI" "\$PREP"
        status2=\$?

        if [ "\$status2" -eq 2 ]; then
          if [ "${strict}" = "true" ]; then
            echo "ERROR: intervals contigs do not match reference (.fai) even after chr autofix" >&2
            exit 1
          else
            echo "WARNING: contig mismatch vs .fai; continuing (bed_validate_strict=false)" >&2
          fi
        fi
      else
        if [ "${strict}" = "true" ]; then
          echo "ERROR: intervals contigs do not match reference (.fai) (autofix disabled)" >&2
          exit 1
        else
          echo "WARNING: contig mismatch vs .fai; continuing (bed_validate_strict=false)" >&2
          cp "\$CLEAN" "\$PREP"
        fi
      fi
    else
      # status 0 => ok, just use clean as prepared
      cp "\$CLEAN" "\$PREP"
    fi

    # 4) Chunking
    awk -v FS="\\t" \\
        -v nps=${nps} \\
        -v chunk_target=${chunk_target} \\
        -v longest_factor=${longest_factor} \\
        '{
            # t: numeric col5 if present, else length/nps
            if (NF>=5 && \$5 ~ /^[0-9]+(\\.[0-9]+)?\$/) t=\$5;
            else t=(\$3-\$2)/nps;

            if (name == "" || (chunk > chunk_target && (chunk + t) > longest * longest_factor)) {
                name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3);
                chunk = 0;
                longest = 0;
            }
            if (t > longest) longest = t;
            chunk += t;
            print \$0 > name;
        }' "\$PREP"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
