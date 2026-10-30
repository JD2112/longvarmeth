process HIPHASE {
    tag "$sample_id"
    label 'process_high'

    conda "bioconda::hiphase=1.4.0"
    container "quay.io/biocontainers/hiphase:1.4.0--h9ee0642_0"

    input:
    tuple val(sample_id), path(bam), path(bai), path(small_vcf), path(small_tbi), path(sv_vcf), path(sv_tbi), path(str_vcf), path(str_tbi)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.phased.bam"), path("${sample_id}.phased.bam.bai")                    , emit: phased_bam_bai
    tuple val(sample_id), path("${sample_id}.small_phased.vcf.gz"), path("${sample_id}.small_phased.vcf.gz.tbi")  , emit: small_vcf_phased
    tuple val(sample_id), path("${sample_id}.sv_phased.vcf.gz"), path("${sample_id}.sv_phased.vcf.gz.tbi")        , optional: true, emit: sv_vcf_phased
    tuple val(sample_id), path("${sample_id}.str_phased.vcf.gz"), path("${sample_id}.str_phased.vcf.gz.tbi")      , optional: true, emit: str_vcf_phased
    tuple val(sample_id), path("${sample_id}.hiphase.stats.csv")                                                   , optional: true, emit: stats_csv
    path "versions.yml"                                                                                           , emit: versions

    script:
    def sv_args = (sv_vcf && sv_vcf.name != 'NO_FILE') ? "--vcf ${sv_vcf}:${sample_id}.sv_phased.vcf.gz" : ""
    def str_args = (str_vcf && str_vcf.name != 'NO_FILE') ? "--vcf ${str_vcf}:${sample_id}.str_phased.vcf.gz" : ""
    """
    hiphase \\
        --bam ${bam} \\
        --output-bam ${sample_id}.phased.bam \\
        --reference ${reference} \\
        --vcf ${small_vcf}:${sample_id}.small_phased.vcf.gz \\
        ${sv_args} \\
        ${str_args} \\
        --summary-file ${sample_id}.hiphase.stats.csv \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: \$( hiphase --version | sed 's/hiphase //' )
    END_VERSIONS
    """

    stub:
    def sv_stub = (sv_vcf && sv_vcf.name != 'NO_FILE') ? "touch ${sample_id}.sv_phased.vcf.gz; touch ${sample_id}.sv_phased.vcf.gz.tbi" : ""
    def str_stub = (str_vcf && str_vcf.name != 'NO_FILE') ? "touch ${sample_id}.str_phased.vcf.gz; touch ${sample_id}.str_phased.vcf.gz.tbi" : ""
    """
    touch ${sample_id}.phased.bam
    touch ${sample_id}.phased.bam.bai
    touch ${sample_id}.small_phased.vcf.gz
    touch ${sample_id}.small_phased.vcf.gz.tbi
    touch ${sample_id}.hiphase.stats.csv
    ${sv_stub}
    ${str_stub}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hiphase: 1.4.0
    END_VERSIONS
    """
}
