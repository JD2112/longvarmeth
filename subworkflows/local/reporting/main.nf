// subworkflows/local/reporting/main.nf
// Step 12: Reporting with MultiQC and Quarto executive summary

include { MULTIQC } from '../../../modules/nf-core/multiqc'
include { DUMP_MQC_METADATA } from '../../../modules/local/dump_mqc_metadata'

workflow REPORTING {
    take:
    qc_files // channel of QC files
    versions // channel of collected versions

    main:
    ch_versions = Channel.empty()

    // 1. Dump Pipeline Metadata for MultiQC
    def methods_template = file("${projectDir}/assets/methods_description.yml")
    DUMP_MQC_METADATA(methods_template)

    // 2. MultiQC Dashboard
    def multiqc_files = Channel.empty()
        .concat(qc_files)
        .concat(DUMP_MQC_METADATA.out.summary)
        .concat(DUMP_MQC_METADATA.out.versions)
        .concat(DUMP_MQC_METADATA.out.methods)
        .concat(versions)
        .collect()

    def multiqc_config = file("${projectDir}/assets/multiqc_config.yml")
    MULTIQC(multiqc_files, multiqc_config)

    emit:
    multiqc_report = MULTIQC.out.report
    multiqc_data   = MULTIQC.out.data
}
