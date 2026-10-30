process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai
    path "versions.yml", emit: versions

    script:
    """
    samtools faidx ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n 1 | sed -e "s/samtools //g" )
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """
}

