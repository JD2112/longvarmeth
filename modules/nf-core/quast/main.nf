process QUAST {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/quast:5.2.0--py39pl5321h4e6d42e_3"

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "quast_${sample_id}", emit: report
    path "versions.yml"      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    quast.py ${assembly} -o quast_${sample_id} --threads ${task.cpus ?: 4} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$( quast.py --version | sed -e "s/QUAST v//g" )
    END_VERSIONS
    """
}
