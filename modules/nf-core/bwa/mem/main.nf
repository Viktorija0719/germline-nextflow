process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam")  , emit: bam,    optional: true
    tuple val(meta), path("*.cram") , emit: cram,   optional: true
    tuple val(meta), path("*.csi")  , emit: csi,    optional: true
    tuple val(meta), path("*.crai") , emit: crai,   optional: true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def sort_mem = task.ext.sort_mem ?: '768M'

    def extension = args2.contains("--output-fmt sam")   ? "sam" :
                    args2.contains("--output-fmt cram")  ? "cram":
                    sort_bam && args2.contains("-O cram")? "cram":
                    !sort_bam && args2.contains("-C")    ? "cram":
                    "bam"

    def reference = fasta && extension=="cram" ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    // Temp handling for samtools sort (fixes tmp.0001.bam errors)
    def tmp_setup = sort_bam ? """
      mkdir -p samtools_sort_tmp
      export TMPDIR="\$PWD/samtools_sort_tmp"
    """ : ""

    def tmp_args = sort_bam ? "-T samtools_sort_tmp/${prefix}.tmp -m ${sort_mem}" : ""

    """
    INDEX=\$(find -L ./ -name "*.amb" -print -quit | sed 's/\\.amb\$//')
    if [[ -z "\$INDEX" ]]; then
      echo "[ERROR] BWA index (.amb) not found in task directory" >&2
      exit 1
    fi

    ${tmp_setup}

    bwa mem \\
        ${args} \\
        -t ${task.cpus} \\
        "\$INDEX" \\
        ${reads} \\
      | samtools ${samtools_command} ${args2} ${reference} ${tmp_args} --threads ${task.cpus} -o ${prefix}.${extension} -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args2_stub = task.ext.args2 ?: ''
    def prefix_stub = task.ext.prefix ?: "${meta.id}"
    def extension_stub = args2_stub.contains("--output-fmt sam")   ? "sam" :
                         args2_stub.contains("--output-fmt cram")  ? "cram":
                         sort_bam && args2_stub.contains("-O cram")? "cram":
                         !sort_bam && args2_stub.contains("-C")    ? "cram":
                         "bam"
    """
    touch ${prefix_stub}.${extension_stub}
    touch ${prefix_stub}.csi
    touch ${prefix_stub}.crai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
