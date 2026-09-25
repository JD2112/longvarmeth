process FILTLONG {
    tag "$sample_id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/filtlong:0.2.1--he1b5a44_3"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: filtered_reads
    path "versions.yml"                                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    filtlong \\
        --min_length ${params.filt_min_len} \\
        --keep_percent ${params.filt_keep_pct} \\
        ${args} \\
        ${reads} | gzip > ${sample_id}_filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filtlong: \$( filtlong --version | sed -e "s/Filtlong v//g" )
    END_VERSIONS
    """
}
