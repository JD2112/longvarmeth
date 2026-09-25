process BUSCO {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/busco:5.5.0--py310h3289130_0"

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "busco_${sample_id}", emit: report
    path "versions.yml"      , emit: versions

    script:
    def lineage = params.busco_lineage ?: 'eukaryota'
    def args = task.ext.args ?: ''
    """
    busco \\
        -i ${assembly} \\
        -o busco_${sample_id} \\
        -m genome \\
        -l ${lineage} \\
        --cpu ${task.cpus ?: 4} \\
        --offline \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | head -n 1 | sed -e "s/BUSCO //g" )
    END_VERSIONS
    """
}
