process BEDTOOLS_BAMTOFASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3':
        'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*R1.fastq"), path("*R2.fastq"), emit: fastq
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools \\
        bamtofastq \\
        -i ${bam} \\
        -fq ${prefix}_R1.fastq \\
        -fq2 ${prefix}_R2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version |& sed '1!d ; s/bedtools v//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1.fastq
    touch ${prefix}_R2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version |& sed '1!d ; s/bedtools v//')
    END_VERSIONS
    """
}
