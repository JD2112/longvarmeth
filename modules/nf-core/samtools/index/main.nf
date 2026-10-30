process SAMTOOLS_INDEX {
    tag "$sample_id"
    label 'process_single'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("*.bai"), emit: bam_bai
    path "versions.yml"                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools index -@ ${task.cpus ?: 2} ${args} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n 1 | sed -e "s/samtools //g" )
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """
}

