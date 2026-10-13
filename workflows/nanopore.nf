// workflows/nanopore.nf
// Native Oxford Nanopore pipeline entrypoint implementing the 12-step architecture

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { READ_QC } from '../subworkflows/local/qc'
include { BASECALL_DEMUX } from '../subworkflows/local/basecall_demux'
include { DE_NOVO_ASSEMBLY } from '../subworkflows/local/assembly'
include { ALIGNMENT } from '../subworkflows/local/alignment'
include { VARIANT_CALLING } from '../subworkflows/local/variant_calling'
include { STRUCTURAL_VARIANTS } from '../subworkflows/local/structural_variants'
include { METHYLATION } from '../subworkflows/local/methylation'
include { PHASING } from '../subworkflows/local/phasing'
include { ANNOTATION } from '../subworkflows/local/annotation'
include { CLINICAL_GOVERNANCE } from '../subworkflows/local/governance'
include { VISUALIZATION } from '../subworkflows/local/visualization'
include { REPORTING } from '../subworkflows/local/reporting'
include { MERGE_BAMS } from '../modules/local/merge_bams'

workflow NANOPORE_PIPELINE {
    def samplesheet_file = file(params.input)
    if (!samplesheet_file.exists()) {
        error "Samplesheet file not found at: ${params.input}"
    }

    ch_versions = Channel.empty()
    ch_qc_files = Channel.empty()

    def samples_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.run_id,
                row.input_path,
                row.input_type.toLowerCase().trim(),
                row.kit ?: '',
                row.model ?: '',
                row.modified_bases_model ?: ''
            )
        }

    // Step 1: Input Checksum
    INPUT_CHECK(samples_ch)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    def raw_runs = samples_ch.filter { it[2] == 'pod5' || it[2] == 'fast5' }
    def direct_bam_runs = samples_ch.filter { it[2] == 'bam' }

    // Step 3 (Raw ONT): Dorado GPU Basecalling & Demultiplexing
    BASECALL_DEMUX(raw_runs)
    def basecalled_bams = BASECALL_DEMUX.out.demux_bams
    ch_versions = ch_versions.mix(BASECALL_DEMUX.out.versions)

    def direct_bams_raw = direct_bam_runs.map { run_id, input_path, input_type, kit, model, mod_model ->
        tuple(run_id, file(input_path))
    }
    MERGE_BAMS(direct_bams_raw)
    ch_versions = ch_versions.mix(MERGE_BAMS.out.versions)
    def direct_bams = MERGE_BAMS.out.bam

    def all_demux_bams = basecalled_bams.concat(direct_bams)

    // Step 2: Read QC
    READ_QC(all_demux_bams)
    ch_versions = ch_versions.mix(READ_QC.out.versions)
    def fastq_ch = READ_QC.out.fastq
    ch_qc_files = ch_qc_files
        .concat(READ_QC.out.fastqc_files)
        .concat(READ_QC.out.nanoplot_files)
        .concat(READ_QC.out.cramino_files)

    // Step 3: De Novo Assembly (Flye)
    if (params.run_type == 'de_novo_assembly') {
        DE_NOVO_ASSEMBLY(fastq_ch)
        ch_versions = ch_versions.mix(DE_NOVO_ASSEMBLY.out.versions)
        ch_qc_files = ch_qc_files
            .concat(DE_NOVO_ASSEMBLY.out.quast_reports)
            .concat(DE_NOVO_ASSEMBLY.out.busco_reports)
    }

    // Reference-guided Steps 4 - 11
    if (params.run_type == 'reference_based' && params.ref_nuc) {
        def ref_nuc_path = file(params.ref_nuc)
        def ref_fai_path = file("${params.ref_nuc}.fai")

        // Step 4: Minimap2 Alignment & Tag Restoration (skip if input is pre-aligned)
        def aligned_bams = Channel.empty()
        def pre_aligned_runs = samples_ch.filter { it[2] == 'aligned_bam' }
        def need_alignment_runs = samples_ch.filter { it[2] != 'aligned_bam' }

        if (pre_aligned_runs) {
            def pre_aligned_bams = pre_aligned_runs.map { run_id, input_path, input_type, kit, model, mod_model ->
                def bam_file = file(input_path)
                def bai_file = file("${input_path}.bai")
                if (!bai_file.exists()) {
                    bai_file = file(input_path.replaceAll(/\.bam$/, '.bai'))
                }
                tuple(run_id, bam_file, bai_file)
            }
            aligned_bams = aligned_bams.mix(pre_aligned_bams)
        }

        if (need_alignment_runs && !params.skip_alignment) {
            ALIGNMENT(fastq_ch, all_demux_bams, ref_nuc_path)
            ch_versions = ch_versions.mix(ALIGNMENT.out.versions)
            aligned_bams = aligned_bams.mix(ALIGNMENT.out.aligned_bam_bai)
        }

        // Step 5: Small Variants
        def small_vcfs = Channel.empty()
        if (!params.skip_variants) {
            VARIANT_CALLING(aligned_bams, ref_nuc_path, ref_fai_path)
            ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)
            small_vcfs = VARIANT_CALLING.out.vcf_tbi
            ch_qc_files = ch_qc_files.concat(small_vcfs.map { it[1] })
        }

        // Step 6: Structural Variants & CNVs (Sniffles2 / CNVkit)
        def sv_vcfs = Channel.empty()
        if (!params.skip_variants) {
            STRUCTURAL_VARIANTS(aligned_bams, ref_nuc_path)
            ch_versions = ch_versions.mix(STRUCTURAL_VARIANTS.out.versions)
            sv_vcfs = STRUCTURAL_VARIANTS.out.sv_vcf
        }

        // Step 7: Methylation (Modkit pileup)
        def bed_methylation = Channel.empty()
        if (!params.skip_methylation) {
            METHYLATION(aligned_bams, ref_nuc_path)
            ch_versions = ch_versions.mix(METHYLATION.out.versions)
            bed_methylation = METHYLATION.out.bed_methylation
        }

        // Step 8: WhatsHap Phasing & Haplotagging
        def phased_bams = aligned_bams
        def phased_vcfs = small_vcfs
        if (!params.skip_variants && (params.species == 'human' || params.enable_phasing)) {
            PHASING(aligned_bams, small_vcfs, sv_vcfs, ref_nuc_path, ref_fai_path)
            ch_versions = ch_versions.mix(PHASING.out.versions)
            phased_bams = PHASING.out.phased_bam_bai
            phased_vcfs = PHASING.out.phased_vcf
        }

        // Step 9: Annotation
        def final_vcfs = phased_vcfs
        if (!params.skip_variants) {
            ANNOTATION(phased_vcfs)
            ch_versions = ch_versions.mix(ANNOTATION.out.versions)
            final_vcfs = ANNOTATION.out.annotated_vcf
        }

        // Step 10: Clinical Governance
        if (params.enable_audit_trail || params.enable_phi_sanitization) {
            def audit_targets = aligned_bams
                .map { sample_id, bam, bai ->
                    tuple(sample_id, [bam, bai])
                }
            CLINICAL_GOVERNANCE(aligned_bams, audit_targets)
            ch_versions = ch_versions.mix(CLINICAL_GOVERNANCE.out.versions)
        }

        // Step 11: Visualization (JBrowse 2, Methylartist & Methylation Dashboard)
        def ref_tuple_ch = fastq_ch.map { sample_id, fastq ->
            tuple(sample_id, [ref_nuc_path, ref_fai_path])
        }
        VISUALIZATION(ref_tuple_ch, phased_bams, final_vcfs, bed_methylation)
        ch_versions = ch_versions.mix(VISUALIZATION.out.versions)
        ch_qc_files = ch_qc_files.mix(VISUALIZATION.out.methylation_html.map { meta, html -> html })
    }

    // Step 12: Reporting
    REPORTING(ch_qc_files, ch_versions)
}
