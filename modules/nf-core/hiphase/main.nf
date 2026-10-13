process HIPHASE {
    tag "$sample_id"
    label 'process_high'

    conda "bioconda::hiphase=1.4.0"
    container "quay.io/biocontainers/hiphase:1.4.0--h9ee0642_0"

    input:
    tuple val(sample_id), path(bam), path(bai), path(small_vcf), path(small_tbi), path(sv_vcf), path(sv_tbi)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.phased.bam"), path("${sample_id}.phased.bam.bai")      , emit: phased_bam_bai
    tuple val(sample_id), path("${sample_id}.small_phased.vcf.gz"), path("${sample_id}.small_phased.vcf.gz.tbi"), emit: small_vcf_phased
    tuple val(sample_id), path("${sample_id}.sv_phased.vcf.gz"), path("${sample_id}.sv_phased.vcf.gz.tbi")      , optional: true, emit: sv_vcf_phased
    path "versions.yml"                                                                             , emit: versions

    script:
    def sv_args = sv_vcf ? "--vcf-sv ${sv_vcf}" : ""
    """
    hiphase \\
        --bam ${bam} \\
        --vcf ${small_vcf} \\
        ${sv_args} \\
        --reference ${reference} \\
        --output-bam ${sample_id}.phased.bam \\
        --output-vcf ${sample_id}.small_phased.vcf.gz \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$( hiphase --version | sed 's/hiphase //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.phased.bam
    touch ${sample_id}.phased.bam.bai
    touch ${sample_id}.small_phased.vcf.gz
    touch ${sample_id}.small_phased.vcf.gz.tbi
    touch ${sample_id}.sv_phased.vcf.gz
    touch ${sample_id}.sv_phased.vcf.gz.tbi
    touch versions.yml
    """
}
