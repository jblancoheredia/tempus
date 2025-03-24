process TELOMEREHUNTER_FULL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/telomerehunter:1.1.0' :
        'blancojmskcc/telomerehunter:1.1.0' }"

    input:
    tuple val(meta) , path(tumour_bam)
    tuple val(meta1), path(tumour_bai)
    path(normal_bam)
    path(normal_bai)
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

    telomerehunter --outPath ${outputDir}/ \\
    --pid ${prefix} \\
    --inputBamTumor ${tumour_bam} \\
    --inputBamControl ${normal_bam} \\
    --bandingFile ${banding} \\
    --repeatThreshold 4 \\
    --mappingQualityThreshold 4 \\
    --repeats TTAGGG TGAGGG TCAGGG TTCGGG TTGGGG TTTGGG ATAGGG CATGGG CTAGGG GTAGGG TAAGGG \\
    --repeatsContext TTAGGG TGAGGG TCAGGG TTCGGG TTGGGG TTTGGG ATAGGG CATGGG CTAGGG GTAGGG TAAGGG \\
    --parallel \\
    -p1 -p2 -p3 -p4 -p5 -p6 -p7

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
