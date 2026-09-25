process BCFTOOLS_CALL {
    tag "$sample_id ($region)"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0"

    input:
    tuple val(sample_id), path(bam), path(bai), val(region)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_${region}.vcf.gz"), path("${sample_id}_${region}.vcf.gz.tbi"), emit: vcf_tbi
    path "versions.yml"                                                                                    , emit: versions

    script:
    def ploidy_flag = params.ploidy == 1 ? '--ploidy 1' : ''
    def args = task.ext.args ?: ''
    """
    bcftools mpileup --threads ${task.cpus} -Ou -r "${region}" -f ${reference} -a FORMAT/DP ${bam} | \\
    bcftools call --threads ${task.cpus} -mv ${ploidy_flag} -f GQ -Ob -o ${sample_id}_${region}.raw.bcf
    
    bcftools filter -e 'QUAL<40 || INFO/DP<3 || INFO/DP>50 || INFO/MQ<40 || FMT/GQ<30' ${args} -Oz -o ${sample_id}_${region}.vcf.gz ${sample_id}_${region}.raw.bcf
    bcftools index -t ${sample_id}_${region}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version 2>&1 | head -n 1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """
}
