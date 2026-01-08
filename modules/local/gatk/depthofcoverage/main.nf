// process GATK_DEPTHOFCOVERAGE {
//     tag "${meta.id}"
//     label 'process_medium'

//     conda "${moduleDir}/environment.yml"
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'docker://broadinstitute/gatk:4.5.0.0' :
//         'docker.io/broadinstitute/gatk:4.5.0.0' }"

process GATK_DEPTHOFCOVERAGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://broadinstitute/gatk:4.5.0.0' :
        'docker.io/broadinstitute/gatk:4.5.0.0' }"

    input:
    // BAM + BAI (deduped BAM + index, same meta as in other modules)
    tuple val(meta), path(bam), path(bai)

    // Reference FASTA, FAI, DICT
    path ref_fasta
    path ref_fai
    path ref_dict

    // Target regions (BED)
    path target_bed

    // Gene list for --calculate-coverage-over-genes
    path gene_list

    output:
    tuple val(meta), path("${prefix}*"), emit: coverage
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.DepthOfCoverage"

    """
    set -euo pipefail

    # Normalize reference aux file names for GATK expectations
    ref_base=\$(basename ${ref_fasta})          # Homo_sapiens_assembly38.fasta
    ref_nosuf=\${ref_base%.fasta}              # Homo_sapiens_assembly38

    # Ensure FAI has the expected name <ref_base>.fai in this work dir
    if [ ! -e "\${ref_base}.fai" ]; then
        ln -s ${ref_fai} "\${ref_base}.fai"
    fi

    # Ensure DICT has the expected name <ref_nosuf>.dict in this work dir
    if [ ! -e "\${ref_nosuf}.dict" ]; then
        ln -s ${ref_dict} "\${ref_nosuf}.dict"
    fi

    # GATK4 DepthOfCoverage (BETA)
    gatk \\
        DepthOfCoverage \\
        --output-format TABLE \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        -R ${ref_fasta} \\
        -O ${prefix} \\
        -I ${bam} \\
        --calculate-coverage-over-genes ${gene_list} \\
        --omit-depth-output-at-each-base \\
        --omit-locus-table \\
        --start 1 \\
        --stop 5000 \\
        --nBins 200 \\
        --include-ref-n-sites \\
        --read-filter MappingQualityReadFilter \\
        --minimum-mapping-quality 5 \\
        -L ${target_bed} \\
        ${extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.DepthOfCoverage"

    """
    touch ${prefix}.summary
    touch versions.yml
    """
}
