process CNVKIT {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/cnvkit:0.9.13--pyhdfd78af_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.cns"), path("${sample_id}-diagram.pdf"), emit: cnv_files
    path "versions.yml"                                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    cnvkit.py batch ${bam} --normal -m wgs --fasta ${reference} --output-dir . --diagram -p ${task.cpus} ${args}
    
    bam_base=\$(basename ${bam} .bam)
    if [ "\${bam_base}" != "${sample_id}" ]; then
        mv "\${bam_base}.cns" "${sample_id}.cns"
        mv "\${bam_base}-diagram.pdf" "${sample_id}-diagram.pdf"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$( cnvkit.py version | sed -e "s/cnvkit //g" )
    END_VERSIONS
    """
}
