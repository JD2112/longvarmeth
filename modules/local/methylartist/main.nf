process METHYLARTIST_LOCUS {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/methylartist:1.3.0--pyhdfd78af_1"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    val locus

    output:
    tuple val(sample_id), path("${sample_id}_locus_*.png"), emit: plot
    path "versions.yml"                                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    methylartist locus \\
        -b ${bam} \\
        -i "${locus}" \\
        --ref ${ref} \\
        -o ${sample_id}_locus \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        methylartist: \$( methylartist --version 2>&1 | sed -e "s/methylartist //g" )
    END_VERSIONS
    """
}
