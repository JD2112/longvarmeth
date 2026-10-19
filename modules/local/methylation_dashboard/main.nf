process METHYLATION_DASHBOARD {
    tag "$sample_id"
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(sample_id), path(profile_tsv)
    path template_html

    output:
    tuple val(sample_id), path("*.html"), emit: html
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    """
    generate_methylation_dashboard.py \\
        --profile "${profile_tsv}" \\
        --template "${template_html}" \\
        --sample_id "${sample_id}" \\
        --output "${prefix}_methylation_dashboard.html" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        gms_nallo_methylation_dashboard: "0.3.0"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample_id}"
    """
    touch "${prefix}_methylation_dashboard.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.10"
        gms_nallo_methylation_dashboard: "0.3.0"
    END_VERSIONS
    """
}
