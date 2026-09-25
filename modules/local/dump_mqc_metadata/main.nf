process DUMP_MQC_METADATA {
    label 'process_single'

    input:
    path methods_description_template

    output:
    path "pipeline_summary_mqc.yaml"  , emit: summary
    path "software_versions_mqc.yaml" , emit: versions
    path "methods_description_mqc.yaml", emit: methods

    script:
    """
    cat <<EOF > pipeline_summary_mqc.yaml
    id: 'pipeline-summary'
    section_name: 'Workflow Summary'
    plot_type: 'html'
    description: 'Parameters used for this pipeline run.'
    data: |
      <table class="table table-bordered table-striped" style="width: auto;">
        <tr><td>Pipeline</td><td>longvarmeth</td></tr>
        <tr><td>Species</td><td>${params.species}</td></tr>
        <tr><td>Run Type</td><td>${params.run_type}</td></tr>
        <tr><td>Medaka Model</td><td>${params.medaka_model}</td></tr>
        <tr><td>Min Barcode Size (MB)</td><td>${params.min_barcode_size_mb}</td></tr>
        <tr><td>Filt Min Length</td><td>${params.filt_min_len}</td></tr>
        <tr><td>Filt Keep %</td><td>${params.filt_keep_pct}</td></tr>
      </table>
    EOF

    cat <<EOF > software_versions_mqc.yaml
    id: 'software-versions'
    section_name: 'Software Versions'
    plot_type: 'html'
    description: 'Software versions used in this pipeline.'
    data: |
      <table class="table table-bordered table-striped" style="width: auto;">
        <tr><td>Nextflow</td><td>24.04+</td></tr>
        <tr><td>Dorado</td><td>0.9.6</td></tr>
        <tr><td>Minimap2</td><td>2.26</td></tr>
        <tr><td>Samtools</td><td>1.23.1</td></tr>
        <tr><td>WhatsHap</td><td>2.8</td></tr>
        <tr><td>BCFtools</td><td>1.23.1</td></tr>
        <tr><td>Sniffles</td><td>2.8.0</td></tr>
        <tr><td>CNVkit</td><td>0.9.13</td></tr>
        <tr><td>Modkit</td><td>0.6.3</td></tr>
        <tr><td>Methylartist</td><td>1.3.0</td></tr>
        <tr><td>JBrowse 2</td><td>4.3.0</td></tr>
      </table>
    EOF

    cp ${methods_description_template} methods_description_mqc.yaml
    """
}
