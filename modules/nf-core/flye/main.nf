process FLYE {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/flye:2.9.2--py310h2b65f2a_3"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}/assembly.fasta"), emit: assembly
    path "${sample_id}/*"                                    , emit: flye_files
    path "versions.yml"                                      , emit: versions

    script:
    def args = task.ext.args ?: '--nano-raw'
    """
    flye ${args} ${reads} --out-dir ${sample_id} --threads ${task.cpus ?: 8}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """
}
