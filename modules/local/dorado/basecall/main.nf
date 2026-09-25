process DORADO_BASECALL {
    tag "$run_id"
    label 'process_gpu'

    input:
    tuple val(run_id), path(pod5_dir), val(kit), val(model), val(modified_bases_model)

    output:
    tuple val(run_id), path("basecalled_reads.bam"), emit: bam

    script:
    def mod_bases_arg = modified_bases_model ? "--modified-bases-models ${modified_bases_model}" : ""
    """
    mkdir -p models
    ${params.dorado_bin} basecaller \\
        --models-directory models \\
        ${model} \\
        ${pod5_dir} \\
        -x ${params.gpu_devices} \\
        --emit-moves \\
        ${mod_bases_arg} \\
        > basecalled_reads.bam
    """
}
