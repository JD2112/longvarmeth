process DORADO_DEMUX {
    tag "$run_id"
    label 'process_medium'

    input:
    tuple val(run_id), path(basecalled_bam), val(kit)

    output:
    tuple val(run_id), path("demux/*.bam"), emit: bams

    script:
    """
    mkdir -p demux
    ${params.dorado_bin} demux \\
        --kit-name ${kit} \\
        ${basecalled_bam} \\
        --output-dir demux \\
        --threads ${task.cpus ?: params.threads}
    """
}
