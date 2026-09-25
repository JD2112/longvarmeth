process BCFTOOLS_CONCAT {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0"

    input:
    tuple val(sample_id), path(vcfs), path(tbis)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf_tbi
    path "versions.yml"                                                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bcftools concat --threads ${task.cpus} -a -d none ${args} -Oz -o ${sample_id}.vcf.gz ${vcfs}
    bcftools index -t ${sample_id}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version 2>&1 | head -n 1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """
}
