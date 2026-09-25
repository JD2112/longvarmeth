process JBROWSE_CREATE {
    tag "$sample_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/jbrowse2:4.3.0--hdfd78af_0"

    input:
    tuple val(sample_id), path(assembly_files), path(bam), path(bai), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}"), emit: jbrowse_dir
    path "versions.yml"                       , emit: versions

    script:
    def fasta = (assembly_files instanceof List) ? assembly_files.find { it.name.endsWith('.fa') || it.name.endsWith('.fasta') || it.name.endsWith('.fna') } : assembly_files
    """
    jbrowse create ${sample_id}
    if [ -f "${fasta}" ]; then
        jbrowse add-assembly ${fasta} --load copy --name "${sample_id}_assembly" --out ${sample_id}
    fi
    if [ -f "${bam}" ]; then
        jbrowse add-track ${bam} --load copy --name "${sample_id}_bam" --out ${sample_id}
    fi
    if [ -f "${vcf}" ]; then
        jbrowse add-track ${vcf} --load copy --name "${sample_id}_vcf" --out ${sample_id}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jbrowse: \$( jbrowse --version )
    END_VERSIONS
    """
}
