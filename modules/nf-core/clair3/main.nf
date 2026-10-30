process CLAIR3 {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/clair3:1.0.10--py39h6d43e06_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    val platform // 'ont', 'hifi'

    output:
    tuple val(sample_id), path("${sample_id}.clair3.vcf.gz"), path("${sample_id}.clair3.vcf.gz.tbi"), emit: vcf_tbi
    path "versions.yml"                                                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def plat = platform ?: 'ont'
    """
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus ?: 4} \\
        --platform=${plat} \\
        --model_path="\${CONDA_PREFIX}/bin/models/${plat}" \\
        --output=clair3_out \\
        --enable_long_indel \\
        ${args}

    mv clair3_out/merge_output.vcf.gz ${sample_id}.clair3.vcf.gz
    mv clair3_out/merge_output.vcf.gz.tbi ${sample_id}.clair3.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$( run_clair3.sh --version 2>&1 | head -n1 | sed -e 's/Clair3 //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.clair3.vcf.gz
    touch ${sample_id}.clair3.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: 1.0.10
    END_VERSIONS
    """
}

