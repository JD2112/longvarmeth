process TRGT {
    tag "$sample_id"
    label 'process_medium'

    conda "bioconda::trgt=1.0.0"
    container "quay.io/biocontainers/trgt:1.0.0--h9ee0642_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path repeat_bed

    output:
    tuple val(sample_id), path("${sample_id}.trgt.vcf.gz"), path("${sample_id}.trgt.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(sample_id), path("${sample_id}.spanning.bam")                                     , emit: spanning_bam
    path "versions.yml"                                                                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    trgt \\
        --genome ${reference} \\
        --repeats ${repeat_bed} \\
        --reads ${bam} \\
        --threads ${task.cpus} \\
        --output-prefix ${sample_id}.trgt \\
        ${args}

    mv ${sample_id}.trgt.spanning.bam ${sample_id}.spanning.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trgt: \$( trgt --version | sed 's/trgt //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.trgt.vcf.gz
    touch ${sample_id}.trgt.vcf.gz.tbi
    touch ${sample_id}.spanning.bam
    touch versions.yml
    """
}
