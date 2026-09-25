process SNPEFF {
    tag "$sample_id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "quay.io/biocontainers/snpeff:5.2--hdfd78af_1"

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path ready_signal
    
    output:
    tuple val(sample_id), path("${sample_id}.snpeff.vcf"), emit: vcf
    path "${sample_id}.snpeff.genes.txt"                 , emit: genes_txt
    path "${sample_id}.snpeff.summary.html"              , emit: summary_html
    path "versions.yml"                                  , emit: versions

    script:
    def raw_cache_dir = params.snpeff_cache ? file(params.snpeff_cache).toAbsolutePath().toString() : "${projectDir}/db/snpeff"
    def cache_dir = file("${raw_cache_dir}/data/${params.snpeff_db}").exists() ? "${raw_cache_dir}/data" : raw_cache_dir
    def args = task.ext.args ?: ''
    """
    export JAVA_TOOL_OPTIONS="-XX:-UsePerfData"

    snpEff -Xmx8g -quiet -dataDir "${cache_dir}" ${args} ${params.snpeff_db} ${vcf} > ${sample_id}.snpeff.vcf

    mv snpEff_genes.txt ${sample_id}.snpeff.genes.txt
    mv snpEff_summary.html ${sample_id}.snpeff.summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$( snpEff -version 2>&1 | head -n 1 | sed -e "s/snpEff //g" )
    END_VERSIONS
    """
}
