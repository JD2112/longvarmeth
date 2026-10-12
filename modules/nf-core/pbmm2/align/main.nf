process PBMM2_ALIGN {
    tag "$sample_id"
    label 'process_high'

    conda "bioconda::pbmm2=1.13.0 bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/pbmm2:1.13.0--h9ee0642_0"

    input:
    tuple val(sample_id), path(reads)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: bam_bai
    path "versions.yml"                                                                          , emit: versions

    script:
    def args = task.ext.args ?: '--preset HIFI --sort'
    """
    pbmm2 align ${args} -j ${task.cpus} ${reference} ${reads} ${sample_id}.aligned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$( pbmm2 --version 2>&1 | head -n1 | sed 's/pbmm2 //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.aligned.bam
    touch ${sample_id}.aligned.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: 1.13.0
    END_VERSIONS
    """
}
