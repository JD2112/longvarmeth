process VEP {
    tag "$sample_id"
    label 'process_high'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/ensembl-vep:116.0--pl5321h2a3209d_0"

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path ready_signal
    
    output:
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz"), path("${sample_id}.vep.vcf.gz.tbi"), emit: vcf_tbi
    path "${sample_id}.vep.summary.html"                                                      , emit: summary_html
    path "versions.yml"                                                                       , emit: versions

    script:
    def cache_dir = params.vep_cache ? file(params.vep_cache).toAbsolutePath().toString() : "${projectDir}/db/vep"
    def args = task.ext.args ?: ''
    """
    vep \\
        -i ${vcf} \\
        -o ${sample_id}.vep.vcf \\
        --vcf \\
        --offline \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --dir_cache "${cache_dir}" \\
        --fork ${task.cpus ?: 4} \\
        --stats_file ${sample_id}.vep.summary.html \\
        ${args}

    bgzip -c ${sample_id}.vep.vcf > ${sample_id}.vep.vcf.gz
    tabix -p vcf ${sample_id}.vep.vcf.gz
    rm ${sample_id}.vep.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( vep --help 2>&1 | grep "ensembl-vep" | head -n 1 | sed -e "s/.*: //g" )
    END_VERSIONS
    """
}
