process SAMTOOLS_SORT {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam_bai
    path "versions.yml"                                                                        , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools sort -@ ${task.cpus ?: 4} ${args} -o ${sample_id}.sorted.bam ${sam}
    samtools index ${sample_id}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n 1 | sed -e "s/samtools //g" )
    END_VERSIONS
    """
}
