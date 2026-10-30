process HAPPY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap.py:0.3.15--py27h5c963e8_0':
        'quay.io/biocontainers/hap.py:0.3.15--py27h5c963e8_0' }"

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    path truth_vcf
    path truth_tbi
    path truth_bed
    path ref_fasta
    path ref_fai

    output:
    tuple val(meta), path("${prefix}_happy.summary.csv") , emit: summary
    tuple val(meta), path("${prefix}_happy.extended.csv"), emit: extended
    tuple val(meta), path("${prefix}_happy.vcf.gz")      , emit: vcf
    tuple val(meta), path("${prefix}_happy_mqc.tsv")     , emit: mqc
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    hap.py \\
        ${truth_vcf} \\
        ${query_vcf} \\
        -f ${truth_bed} \\
        -r ${ref_fasta} \\
        -o ${prefix}_happy \\
        --threads ${task.cpus} \\
        ${args}

    python3 -c "
import csv
try:
    with open('${prefix}_happy.summary.csv') as f:
        reader = csv.DictReader(f)
        with open('${prefix}_happy_mqc.tsv', 'w') as out:
            out.write('# id: \\'happy_${prefix}\\'\\n')
            out.write('# section_name: \\'Hap.py Truth-Set Benchmark (${prefix})\\'\\n')
            out.write('# plot_type: \\'table\\'\\n')
            out.write('Type\\tFilter\\tRecall\\tPrecision\\tF1_Score\\tTruth_TP\\tQuery_FP\\tTruth_FN\\n')
            for row in reader:
                out.write(f\\\"{row.get('Type','')}\\t{row.get('Filter','')}\\t{row.get('METRIC.Recall','')}\\t{row.get('METRIC.Precision','')}\\t{row.get('METRIC.F1_Score','')}\\t{row.get('TRUTH.TP','')}\\t{row.get('QUERY.FP','')}\\t{row.get('TRUTH.FN','')}\\n\\\")
except Exception as e:
    with open('${prefix}_happy_mqc.tsv', 'w') as out:
        out.write('# id: \\'happy_${prefix}\\'\\n# section_name: \\'Hap.py Truth-Set Benchmark (${prefix})\\'\\n# plot_type: \\'table\\'\\nType\\tRecall\\tPrecision\\n')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        happy: 0.3.15
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_happy.summary.csv
    touch ${prefix}_happy.extended.csv
    echo "" | gzip > ${prefix}_happy.vcf.gz
    cat <<EOF > ${prefix}_happy_mqc.tsv
# id: 'happy_${prefix}'
# section_name: 'Hap.py Truth-Set Benchmark (${prefix})'
# plot_type: 'table'
Type\tFilter\tRecall\tPrecision\tF1_Score\tTruth_TP\tQuery_FP\tTruth_FN
SNP\tALL\t0.9985\t0.9991\t0.9988\t3500000\t3150\t5250
INDEL\tALL\t0.8850\t0.8920\t0.8885\t510000\t61700\t66200
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        happy: 0.3.15
    END_VERSIONS
    """
}

