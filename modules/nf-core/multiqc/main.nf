process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/multiqc:1.18--pyhdfd78af_0"

    input:
    path(qc_files)
    path(multiqc_config)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "versions.yml"        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def config_arg = multiqc_config ? "--config ${multiqc_config}" : ""
    """
    multiqc . ${config_arg} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
