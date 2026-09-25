process VEP_DOWNLOAD {
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/ensembl-vep:116.0--pl5321h2a3209d_0"
    
    output:
    path "vep_db_ready.txt", emit: ready
    path "versions.yml"    , emit: versions
    
    script:
    def cache_dir = params.vep_cache ? file(params.vep_cache).toAbsolutePath().toString() : "${projectDir}/db/vep"
    """
    mkdir -p "${cache_dir}"
    if [ ! -d "${cache_dir}/${params.vep_species}/${params.vep_version}_${params.vep_assembly}" ]; then
        echo "VEP cache not found in ${cache_dir}. Downloading..."
        vep_install -a cf -s ${params.vep_species} -y ${params.vep_assembly} -c "${cache_dir}" --CONVERT --NO_HTSLIB --NO_UPDATE
    else
        echo "VEP cache already exists in ${cache_dir}."
    fi
    touch vep_db_ready.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( vep --help 2>&1 | grep "ensembl-vep" | head -n 1 | sed -e "s/.*: //g" )
    END_VERSIONS
    """
}
