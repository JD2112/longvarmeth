process FASTQC {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: qc_files
    path "versions.yml"        , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    export _JAVA_OPTIONS="-Xmx64g"
    fastqc --threads ${task.cpus ?: 12} ${args} ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}
