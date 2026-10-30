process WHATSHAP_STATS {
    tag "$sample_id"
    label 'process_single'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/whatshap:2.8--py39h2de1943_0"

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}.phasing_stats.tsv"), emit: stats_tsv
    tuple val(sample_id), path("${sample_id}.blocks.tsv")       , emit: blocks_tsv
    path "versions.yml"                                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    whatshap stats \\
        --tsv=${sample_id}.phasing_stats.tsv \\
        --block-list=${sample_id}.blocks.tsv \\
        ${args} \\
        ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.phasing_stats.tsv
    touch ${sample_id}.blocks.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: 2.8
    END_VERSIONS
    """
}

