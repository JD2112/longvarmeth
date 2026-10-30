process TRUVARI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/truvari:4.2.2--pyhdfd78af_0':
        'quay.io/biocontainers/truvari:4.2.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    path truth_vcf
    path truth_tbi
    path truth_bed
    path ref_fasta
    path ref_fai

    output:
    tuple val(meta), path("${prefix}_truvari/summary.json"), emit: summary
    tuple val(meta), path("${prefix}_truvari_mqc.tsv")    , emit: mqc
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    truvari bench \\
        -b ${truth_vcf} \\
        -s ${query_vcf} \\
        --includebed ${truth_bed} \\
        -f ${ref_fasta} \\
        -o ${prefix}_truvari/ \\
        ${args}

    python3 -c "
import json
try:
    with open('${prefix}_truvari/summary.json') as f:
        d = json.load(f)
    with open('${prefix}_truvari_mqc.tsv', 'w') as out:
        out.write('# id: \\'truvari_${prefix}\\'\\n')
        out.write('# section_name: \\'Truvari SV Truth-Set Benchmark (${prefix})\\'\\n')
        out.write('# plot_type: \\'table\\'\\n')
        out.write('Metric\\tValue\\n')
        for k in ['TP-base', 'TP-comp', 'FP', 'FN', 'precision', 'recall', 'f1']:
            if k in d:
                val = d[k]
                if isinstance(val, float):
                    out.write(f'{k}\\t{val:.4f}\\n')
                else:
                    out.write(f'{k}\\t{val}\\n')
except Exception as e:
    with open('${prefix}_truvari_mqc.tsv', 'w') as out:
        out.write('# id: \\'truvari_${prefix}\\'\\n# section_name: \\'Truvari SV Truth-Set Benchmark (${prefix})\\'\\n# plot_type: \\'table\\'\\nMetric\\tValue\\n')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: 4.2.2
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_truvari
    cat <<EOF > ${prefix}_truvari/summary.json
{
  "TP-base": 8420,
  "TP-comp": 8415,
  "FP": 412,
  "FN": 560,
  "precision": 0.9533,
  "recall": 0.9376,
  "f1": 0.9454
}
EOF

    cat <<EOF > ${prefix}_truvari_mqc.tsv
# id: 'truvari_${prefix}'
# section_name: 'Truvari SV Truth-Set Benchmark (${prefix})'
# plot_type: 'table'
Metric\tValue
TP-base\t8420
TP-comp\t8415
FP\t412
FN\t560
precision\t0.9533
recall\t0.9376
f1\t0.9454
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: 4.2.2
    END_VERSIONS
    """
}

