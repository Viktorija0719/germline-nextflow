process MERGE_BAMS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.22.1--h96c455f_0'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.id}.merged.bam"), emit: bam

    script:
    def nbams    = bams instanceof List ? bams.size() : 1
    def bam_list = (bams instanceof List ? bams : [bams]).collect { "\"${it}\"" }.join(' ')
    def first    = "\"${bams instanceof List ? bams[0] : bams}\""
    """
    if [ ${nbams} -eq 1 ]; then
        ln -s ${first} ${meta.id}.merged.bam
    else
        samtools merge -@ ${task.cpus} -O BAM - ${bam_list} | \\
            samtools sort -@ ${task.cpus} -o ${meta.id}.merged.bam -
    fi
    """
}
