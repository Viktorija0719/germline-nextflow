process ANNOTSV {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/annotsv:3.3.4--py311hdfd78af_1' :
        'quay.io/biocontainers/annotsv:3.3.4--py311hdfd78af_1' }"

    // Needed for Singularity: AnnotSV may write into its annotation directory
    containerOptions "${ workflow.containerEngine == 'singularity' ? '--writable-tmpfs' : ''}"

    input:
    tuple val(meta), path(sv_vcf)
    path(annotations_dir)          // -annotationsDir  (pass [] to use container-bundled annotations)
    path(config_file)              // user config       (pass [] to skip)
    path(panels_tsv)               // -candidateGenesFile (pass [] to skip)

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.annotated.tsv"), emit: tsv
    path "versions.yml",                                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    // AnnotSV 3.x has no -configFile CLI flag; config is read from (in order):
    //   1. $ANNOTSV/etc/AnnotSV/configfile  (read-only in container)
    //   2. ~/.config/AnnotSV/configfile      (user config, writable — overrides system)
    // Copy the user's file to the writable user-config location.
    def cfg_copy    = config_file    ? """if [ -e '${config_file}' ]; then
        mkdir -p ~/.config/AnnotSV && cp '${config_file}' ~/.config/AnnotSV/configfile
    else
        echo '[WARN] AnnotSV config file not accessible in container; using defaults' >&2
    fi""" : ''
    def annot_opt   = annotations_dir ? "-annotationsDir ${annotations_dir}" : ''
    def panels_opt  = panels_tsv     ? "-candidateGenesFile ${panels_tsv}"   : ''
    """
    # Install custom config before running AnnotSV
    ${cfg_copy}

    # AnnotSV 3.x silently skips variants whose FILTER is not '.' or 'PASS'.
    # Reset every FILTER to '.' so all variants (TooCommon, lowQC, …) are
    # annotated, preserving the original value in INFO/ORIG_FILTER so it
    # survives into the AnnotSV TSV and the knotAnnotSV Excel.
    # Also add dummy FORMAT+SAMPLE columns when absent (required by AnnotSV).
    awk 'BEGIN{OFS="\\t"}
         /^##FILTER/              { print; next }
         /^##INFO=<ID=ORIG_FILTER/ { next }
         /^##/                    { print; next }
         /^#CHROM/ {
             print "##INFO=<ID=ORIG_FILTER,Number=1,Type=String,Description=\\"Original FILTER value before AnnotSV pre-processing\\">"
             if (\$0 !~ /FORMAT/) print \$0"\\tFORMAT\\tSAMPLE"
             else print
             next
         }
         {
             orig = \$7
             \$7 = "."
             if (NF < 9) \$0 = \$0"\\tGT\\t./."
             n = split(\$8, info_arr, ";")
             \$8 = "ORIG_FILTER=" orig
             for (i=1; i<=n; i++) \$8 = \$8 ";" info_arr[i]
             print
         }' "${sv_vcf}" > input_fixed.vcf
    input_vcf=input_fixed.vcf

    AnnotSV \\
        ${annot_opt} \\
        -SVinputFile    \${input_vcf} \\
        -outputFile     ${prefix}.annotated \\
        ${panels_opt} \\
        ${args}

    # AnnotSV may write output into a *_AnnotSV sub-directory; flatten it
    if compgen -G '*_AnnotSV/*.tsv' > /dev/null 2>&1; then
        mv *_AnnotSV/*.tsv .
    fi

    # Ensure the expected output name exists.
    # AnnotSV exits 0 with no output when the input VCF is empty;
    # create an empty TSV in that case so Nextflow does not fail.
    if [ ! -f ${prefix}.annotated.tsv ]; then
        mv ${prefix}.annotated*.tsv ${prefix}.annotated.tsv 2>/dev/null || \
        mv *.tsv                    ${prefix}.annotated.tsv 2>/dev/null || \
        touch                       ${prefix}.annotated.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AnnotSV: \$(AnnotSV --version 2>&1 | sed 's/AnnotSV //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotated.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AnnotSV: stub
    END_VERSIONS
    """
}
