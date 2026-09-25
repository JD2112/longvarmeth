process POD5_CONVERT {
    tag "$run_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pod5:0.3.14--py310h0db2ef8_0"

    input:
    tuple val(run_id), path(fast5_dir)

    output:
    tuple val(run_id), path("converted_pod5"), emit: pod5_dir
    path "versions.yml"                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p converted_pod5
    pod5 convert fast5 ${fast5_dir} --output converted_pod5/ --threads ${task.cpus ?: params.threads} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pod5: \$( pod5 --version 2>&1 | sed -e "s/pod5 //g" )
    END_VERSIONS
    """
}
