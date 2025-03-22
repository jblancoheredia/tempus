process TELOMEREHUNTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'blancojmskcc/telomerehunter:1.1.0' }"

    input:
    tuple val(meta) , path(tumour_bam), path(tumour_bai)
    tuple val(meta2), path(normal_bam), path(normal_bai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    telomerehunter \\
        -ibt ${tumour_bam} \\
        -ibc ${normal_bam} \\
        -o . \\
        -p ${prefix} \\
        -rt 4 \\
        -pl \\
        -r \"TTAGGG TGAGGG TCAGGG TTCGGG TTGGGG TTTGGG ATAGGG CATGGG CTAGGG GTAGGG TAAGGG\" \\
        -rc \"TTAGGG TGAGGG TCAGGG TTCGGG TTGGGG TTTGGG ATAGGG CATGGG CTAGGG GTAGGG TAAGGG\" \\
        -bp 6 \\
        -mqt 6 \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: "1.1.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: "1.1.0"
    END_VERSIONS
    """
}
