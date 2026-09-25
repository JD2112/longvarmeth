process MODKIT_PILEUP {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/ont-modkit:0.6.3--h7f49ad2_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}_mods.bed"), emit: bed
    path "versions.yml"                                , emit: versions

    script:
    def args = task.ext.args ?: '--cpg --modified-bases C:m'
    """
    modkit pileup ${bam} ${sample_id}_mods.bed --ref ${ref} --threads ${task.cpus ?: 8} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$( modkit --version | sed -e "s/modkit //g" )
    END_VERSIONS
    """
}
