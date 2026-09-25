process SNIFFLES {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/sniffles:2.8.0--pyhdfd78af_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_svs.vcf.gz"), emit: sv_vcf
    path "versions.yml"                                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    sniffles --input ${bam} --vcf ${sample_id}_svs.vcf.gz --reference ${reference} --threads ${task.cpus ?: 4} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$( sniffles --version )
    END_VERSIONS
    """
}
