process NANOPLOT {
    tag "$sample_id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/nanoplot:1.42.0--pyhdfd78af_0"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "nanoplot_${sample_id}", emit: report
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    NanoPlot \\
        --threads ${task.cpus ?: 4} \\
        --fastq ${reads} \\
        -o nanoplot_${sample_id} \\
        --N50 \\
        --drop_outliers \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$( NanoPlot --version | sed -e "s/NanoPlot //g" )
    END_VERSIONS
    """
}
