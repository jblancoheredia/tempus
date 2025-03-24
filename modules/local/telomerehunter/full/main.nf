process TELOMEREHUNTER_FULL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/telomerehunter:1.1.0' :
        'blancojmskcc/telomerehunter:1.1.0' }"

    input:
    tuple val(meta) , path(tumour_bam), path(tumour_bai)
    tuple val(meta2), path(normal_bam), path(normal_bai)
    path(banding)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outputDir = "${prefix}_telomerehunter"
    """
    mkdir -p ${outputDir}

    telomerehunter -o ${outputDir}/ \\
    -p ${prefix} \\
    -ibt ${tumour_bam} \\
    -ibc ${normal_bam} \\
    -rt 4 \\
    -pl \\
    -r TTAGGG TGAGGG TCAGGG TTCGGG TTGGGG TTTGGG ATAGGG CATGGG CTAGGG GTAGGG TAAGGG \\
    -b ${banding} \\
    -mqt 6

    cp -r ${outputDir}/${prefix}/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: "1.1.0"
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: "1.1.0"
    END_VERSIONS
    """
}
