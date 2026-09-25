include { BASECALL_DEMUX } from '../subworkflows/local/basecall_demux'
include { QC_ASSEMBLY } from '../subworkflows/local/qc_assembly'
include { ALIGN_VARIANT } from '../subworkflows/local/align_variant'
include { ANNOTATION_VIS } from '../subworkflows/local/annotation_vis'
include { MULTIQC } from '../modules/nf-core/multiqc'
include { DUMP_MQC_METADATA } from '../modules/local/dump_mqc_metadata'
include { MERGE_BAMS } from '../modules/local/merge_bams'

workflow NANOPORE_PIPELINE {
    // 1. Parse samplesheet
    def samplesheet_file = file(params.input)
    if (!samplesheet_file.exists()) {
        error "Samplesheet file not found at: ${params.input}"
    }

    ch_versions = Channel.empty()

    def samples_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.run_id,
                row.input_path,
                row.input_type.toLowerCase().trim(),
                row.kit,
                row.model,
                row.modified_bases_model
            )
        }

    // Split input channels based on input type
    def raw_runs = samples_ch.filter { it[2] == 'pod5' || it[2] == 'fast5' }
    def direct_bam_runs = samples_ch.filter { it[2] == 'bam' }

    // Determine if there are raw (pod5/fast5) runs to basecall
    def has_raw_runs = false
    samplesheet_file.splitCsv(header: true).each { row ->
        def type = row.input_type?.toLowerCase()?.trim()
        if (type == 'pod5' || type == 'fast5') {
            has_raw_runs = true
        }
    }

    // 2. Basecall & Demux (for raw pod5/fast5) only if raw runs exist
    def basecalled_bams = Channel.empty()
    if (has_raw_runs) {
        BASECALL_DEMUX(raw_runs)
        basecalled_bams = BASECALL_DEMUX.out.demux_bams
        ch_versions = ch_versions.mix(BASECALL_DEMUX.out.versions)
    }

    // Format direct BAMs to match (sample_id, bam)
    def direct_bams_raw = direct_bam_runs.map { run_id, input_path, input_type, kit, model, modified_bases_model ->
        tuple(run_id, file(input_path))
    }

    MERGE_BAMS(direct_bams_raw)
    ch_versions = ch_versions.mix(MERGE_BAMS.out.versions)
    def direct_bams = MERGE_BAMS.out.bam

    // Combine all demultiplexed/input BAMs
    def all_demux_bams = basecalled_bams.concat(direct_bams)

    // 3. QC & Optional Assembly
    QC_ASSEMBLY(all_demux_bams)
    ch_versions = ch_versions.mix(QC_ASSEMBLY.out.versions)
    def fastq_ch = QC_ASSEMBLY.out.fastq

    // 4. Alignment & Variant Calling
    ALIGN_VARIANT(fastq_ch, all_demux_bams)
    ch_versions = ch_versions.mix(ALIGN_VARIANT.out.versions)

    // 5. Visualisation / Annotation
    def assembly_ch = Channel.empty()
    if (params.run_type == 'de_novo_assembly') {
        assembly_ch = QC_ASSEMBLY.out.assembly
    } else {
        if (params.ref_nuc) {
            def ref_nuc_file = file(params.ref_nuc)
            def ref_nuc_fai = file("${params.ref_nuc}.fai")
            assembly_ch = fastq_ch.map { sample_id, fastq ->
                tuple(sample_id, [ref_nuc_file, ref_nuc_fai])
            }
        }
    }

    if (params.run_type == 'reference_based' && params.ref_nuc) {
        ANNOTATION_VIS(
            assembly_ch,
            ALIGN_VARIANT.out.tagged_nuc_bam_bai.map { tuple(it[0], it[1], it[2]) },
            ALIGN_VARIANT.out.vcf_nuc.map { tuple(it[0], it[1], it[2]) }
        )
        ch_versions = ch_versions.mix(ANNOTATION_VIS.out.versions)
    }

    // 6. MQC Metadata Generation
    def methods_template = file("${projectDir}/assets/methods_description.yml")
    DUMP_MQC_METADATA(methods_template)

    // 7. MultiQC Report
    def multiqc_files = Channel.empty()
        .concat(QC_ASSEMBLY.out.fastqc_files)
        .concat(QC_ASSEMBLY.out.nanoplot_files)
        .concat(QC_ASSEMBLY.out.cramino_files)
        .concat(QC_ASSEMBLY.out.quast_files)
        .concat(QC_ASSEMBLY.out.busco_files)
        .concat(ALIGN_VARIANT.out.vcf_nuc.map { it[1] })
        .concat(DUMP_MQC_METADATA.out.summary)
        .concat(DUMP_MQC_METADATA.out.versions)
        .concat(DUMP_MQC_METADATA.out.methods)
        .concat(ch_versions)
        .collect()

    def multiqc_config = file("${projectDir}/assets/multiqc_config.yml")
    MULTIQC(multiqc_files, multiqc_config)
}
