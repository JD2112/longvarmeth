process SNPEFF_DOWNLOAD {
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/snpeff:5.2--hdfd78af_1"
    
    output:
    path "snpeff_db_ready.txt", emit: ready
    path "versions.yml"        , emit: versions
    
    script:
    def raw_cache_dir = params.snpeff_cache ? file(params.snpeff_cache).toAbsolutePath().toString() : "${projectDir}/db/snpeff"
    def cache_dir = file("${raw_cache_dir}/data/${params.snpeff_db}").exists() ? "${raw_cache_dir}/data" : raw_cache_dir
    """
    mkdir -p "${cache_dir}"
    if [ ! -d "${cache_dir}/${params.snpeff_db}" ]; then
        echo "Database ${params.snpeff_db} not found in ${cache_dir}. Downloading..."
        snpEff -Xmx8g download -dataDir "${cache_dir}" ${params.snpeff_db}
    else
        echo "Database ${params.snpeff_db} already exists in ${cache_dir}."
    fi
    touch snpeff_db_ready.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$( snpEff -version 2>&1 | head -n 1 | sed -e "s/snpEff //g" )
    END_VERSIONS
    """
}
