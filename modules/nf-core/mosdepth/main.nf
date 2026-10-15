process MOSDEPTH {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta

    output:
    tuple val(sample_id), path("${sample_id}.mosdepth.global.dist.txt"), emit: global_txt
    tuple val(sample_id), path("${sample_id}.mosdepth.region.dist.txt"), emit: regions_txt, optional: true
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt")    , emit: summary_txt
    tuple val(sample_id), path("${sample_id}.per-base.bed.gz")        , emit: per_base_bed, optional: true
    tuple val(sample_id), path("${sample_id}.regions.bed.gz")         , emit: regions_bed, optional: true
    tuple val(sample_id), path("${sample_id}.quantized.bed.gz")       , emit: quantized_bed, optional: true
    path "versions.yml"                                               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def reference = fasta ? "--fasta ${fasta}" : ""
    """
    mosdepth \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${args} \\
        ${sample_id} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$( mosdepth --version 2>&1 | sed 's/mosdepth //g' )
    END_VERSIONS
    """
}
