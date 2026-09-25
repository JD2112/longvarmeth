process CRAMINO {
    tag "$sample_id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/cramino:0.13.0--h03027b5_0"

    input:
    tuple val(sample_id), path(bam)

    output:
    path "${sample_id}_cramino.txt", emit: report
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    cramino --threads ${task.cpus ?: 4} ${args} ${bam} > ${sample_id}_cramino.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cramino: \$( cramino --version | sed -e "s/cramino //g" )
    END_VERSIONS
    """
}
