process DEEPVARIANT {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/deepvariant:1.6.1--h9ee0642_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    val model_type // 'PACBIO', 'ONT_R104', 'WGS'

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"), path("${sample_id}.deepvariant.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(sample_id), path("${sample_id}.deepvariant.g.vcf.gz"), path("${sample_id}.deepvariant.g.vcf.gz.tbi"), optional: true, emit: gvcf_tbi
    path "versions.yml"                                                                                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def model = model_type ?: 'ONT_R104'
    """
    run_deepvariant \\
        --model_type=${model} \\
        --ref=${reference} \\
        --reads=${bam} \\
        --output_vcf=${sample_id}.deepvariant.vcf.gz \\
        --output_gvcf=${sample_id}.deepvariant.g.vcf.gz \\
        --num_shards=${task.cpus ?: 4} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$( run_deepvariant --version 2>&1 | head -n1 | sed -e 's/.*version //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.deepvariant.vcf.gz
    touch ${sample_id}.deepvariant.vcf.gz.tbi
    touch ${sample_id}.deepvariant.g.vcf.gz
    touch ${sample_id}.deepvariant.g.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: 1.6.1
    END_VERSIONS
    """
}

