process MEDAKA {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/medaka:1.8.0--py39hcd45472_0"

    input:
    tuple val(sample_id), path(reads), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fasta"), emit: polished_assembly
    path "versions.yml"                                      , emit: versions

    script:
    def model = params.medaka_model ?: 'r1041_e82_400bps_hac_v4.0.0'
    def args = task.ext.args ?: ''
    """
    medaka_consensus \\
        -i ${reads} \\
        -d ${assembly} \\
        -o polished_dir \\
        -t ${task.cpus ?: 8} \\
        -m ${model} \\
        ${args}
    
    mv polished_dir/consensus.fasta ${sample_id}_polished.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version | sed -e "s/medaka //g" )
    END_VERSIONS
    """
}
