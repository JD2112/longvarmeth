process PBCPGTOOLS {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pb-cpg-tools:2.4.0--h9ee0642_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.combined.bed"), emit: bed
    tuple val(sample_id), path("${sample_id}.combined.bw") , optional: true, emit: bigwig
    path "versions.yml"                                    , emit: versions

    script:
    def args = task.ext.args ?: '--min-coverage 4 --model-type consensus'
    """
    aligned_bam_to_cpg_scores \\
        --bam ${bam} \\
        --ref ${reference} \\
        --output-prefix ${sample_id} \\
        --threads ${task.cpus ?: 4} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbcpgtools: \$( aligned_bam_to_cpg_scores --version 2>&1 | head -n1 | sed -e 's/aligned_bam_to_cpg_scores //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.combined.bed
    touch ${sample_id}.combined.bw
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbcpgtools: 2.4.0
    END_VERSIONS
    """
}

