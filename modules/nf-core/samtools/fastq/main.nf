process SAMTOOLS_FASTQ {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: fastq
    path "versions.yml"                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools fastq -@ ${task.cpus ?: 16} ${args} ${bam} | gzip > ${sample_id}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n 1 | sed -e "s/samtools //g" )
    END_VERSIONS
    """
}
