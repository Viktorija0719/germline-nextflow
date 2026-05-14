process BEDTOOLS_MERGE_CNVS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_1' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_1' }"

    input:
    tuple val(meta), path(vcf)
    path fai

    output:
    tuple val(meta), path("${prefix}.merged.bed"), emit: bed
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail

    # Split VCF by SVTYPE using awk (header lines preserved for bedtools)
    awk '/^#/ || /SVTYPE=DUP/' ${vcf} > dup.vcf
    awk '/^#/ || /SVTYPE=DEL/' ${vcf} > del.vcf
    awk '/^#/ || /SVTYPE=INS/' ${vcf} > ins.vcf

    for typ in dup del ins; do
        nvar=\$(grep -vc '^#' \${typ}.vcf 2>/dev/null || echo 0)
        if [ "\${nvar}" -gt 0 ]; then
            bedtools merge -i \${typ}.vcf -c 6,8 -o first,collapse > \${typ}_collapsed.bed
        else
            touch \${typ}_collapsed.bed
        fi
    done

    # Combine and sort by reference chromosome order
    cat dup_collapsed.bed del_collapsed.bed ins_collapsed.bed | \\
        bedtools sort -faidx ${fai} > ${prefix}.merged.bed

    rm -f dup.vcf del.vcf ins.vcf dup_collapsed.bed del_collapsed.bed ins_collapsed.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: stub
    END_VERSIONS
    """
}
