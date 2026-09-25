process MINIMAP2 {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_2"

    input:
    tuple val(sample_id), path(reads)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam
    path "versions.yml"                           , emit: versions

    script:
    def args = task.ext.args ?: '-ax map-ont'
    """
    minimap2 ${args} -t ${task.cpus ?: 8} ${reference} ${reads} > ${sample_id}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$( minimap2 --version )
    END_VERSIONS
    """
}
