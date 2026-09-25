process MERGE_BAMS {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"

    input:
    tuple val(sample_id), path(bam_input)

    output:
    tuple val(sample_id), path("${sample_id}_merged.bam"), emit: bam
    path "versions.yml"                                  , emit: versions

    script:
    """
    if [ -d "${bam_input}" ]; then
        find -L "${bam_input}" -type f -name "*.bam" > bam_list.txt
        num_bams=\$(wc -l < bam_list.txt)
        if [ "\$num_bams" -eq 0 ]; then
            echo "Error: No BAM files found in directory ${bam_input}" >&2
            exit 1
        elif [ "\$num_bams" -eq 1 ]; then
            cp "\$(cat bam_list.txt)" "${sample_id}_merged.bam"
        else
            samtools merge -@ ${task.cpus ?: 4} "${sample_id}_merged.bam" -b bam_list.txt
        fi
        rm -f bam_list.txt
    else
        cp "${bam_input}" "${sample_id}_merged.bam"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n 1 | sed -e "s/samtools //g" )
    END_VERSIONS
    """
}
