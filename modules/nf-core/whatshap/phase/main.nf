process WHATSHAP_PHASE {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/whatshap:2.8--py39h2de1943_0"

    input:
    tuple val(sample_id), path(vcf), path(bam), path(bai)
    path reference_files

    output:
    tuple val(sample_id), path("${sample_id}.phased.vcf.gz"), path("${sample_id}.phased.vcf.gz.tbi"), emit: phased_vcf
    path "versions.yml"                                                                              , emit: versions

    script:
    def reference = (reference_files instanceof List) ? reference_files.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') || it.name.endsWith('.fna') } : reference_files
    def args = task.ext.args ?: ''
    """
    whatshap phase \\
        --reference ${reference} \\
        --ignore-read-groups \\
        ${args} \\
        -o ${sample_id}.phased.vcf \\
        ${vcf} \\
        ${bam}
    
    bgzip ${sample_id}.phased.vcf
    tabix -p vcf ${sample_id}.phased.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
    END_VERSIONS
    """
}
