process HIFINCNV {
    tag "$sample_id"
    label 'process_medium'

    conda "bioconda::hificnv=1.0.0"
    container "quay.io/biocontainers/hificnv:1.0.0--h9ee0642_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path excluded_regions

    output:
    tuple val(sample_id), path("${sample_id}.hificnv.vcf.gz"), path("${sample_id}.hificnv.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(sample_id), path("${sample_id}.hificnv.copynumber.bedgraph")                             , optional: true, emit: bedgraph
    path "versions.yml"                                                                                , emit: versions

    script:
    def exclude_arg = excluded_regions ? "--exclude ${excluded_regions}" : ""
    """
    hificnv \\
        --bam ${bam} \\
        --ref ${reference} \\
        ${exclude_arg} \\
        --threads ${task.cpus} \\
        --output-prefix ${sample_id}.hificnv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: \$( hificnv --version | sed 's/hificnv //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.hificnv.vcf.gz
    touch ${sample_id}.hificnv.vcf.gz.tbi
    touch ${sample_id}.hificnv.copynumber.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hificnv: 1.0.0
    END_VERSIONS
    """
}
