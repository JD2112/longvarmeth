process WHATSHAP_HAPLOTAG {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/whatshap:2.8--py39h2de1943_0"

    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi), path(bam), path(bai)
    path reference_files

    output:
    tuple val(sample_id), path("${sample_id}.haplotagged.bam"), path("${sample_id}.haplotagged.bam.bai"), emit: haplotagged_bam_bai
    path "versions.yml"                                                                                  , emit: versions

    script:
    def reference = (reference_files instanceof List) ? reference_files.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') || it.name.endsWith('.fna') } : reference_files
    def args = task.ext.args ?: ''
    """
    whatshap haplotag \\
        --reference ${reference} \\
        --ignore-read-groups \\
        ${args} \\
        -o ${sample_id}.haplotagged.bam \\
        ${vcf} \\
        ${bam}
    
    python -c "import pysam; pysam.index('${sample_id}.haplotagged.bam')"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$( whatshap --version )
    END_VERSIONS
    """
}
