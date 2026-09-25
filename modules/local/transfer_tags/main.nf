process TRANSFER_TAGS {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pysam:0.21.0--py310hc2720c6_1"

    input:
    tuple val(sample_id), path(mapped_bam), path(mapped_bai), path(source_bam)

    output:
    tuple val(sample_id), path("${sample_id}.mmtags.sorted.bam"), path("${sample_id}.mmtags.sorted.bam.bai"), emit: tagged_bam_bai
    path "versions.yml"                                                                                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    transfer_mod_tags.py --source ${source_bam} --mapped ${mapped_bam} --out ${sample_id}.mmtags.sorted.bam ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$( python -c "import pysam; print(pysam.__version__)" )
    END_VERSIONS
    """
}
