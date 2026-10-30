// workflows/longvarmeth.nf
// Unified Long-Read Analysis Pipeline (Oxford Nanopore & PacBio HiFi)
// Structured across 12 distinct, clinical-grade subworkflows

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
include { MOSDEPTH } from '../modules/nf-core/mosdepth'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index'
include { BENCHMARKING } from '../subworkflows/local/benchmarking'

workflow LONGVARMETH {
    // Check samplesheet exists
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

    // =========================================================================
    // STEP 1: Input File Checksum & Ingestion Verification
    // =========================================================================
    INPUT_CHECK(samples_ch)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    def ch_input_checksums = INPUT_CHECK.out.checksums

    // Branch input samples by input type
    def ch_branched = samples_ch.branch {
        raw_signal:  it[2] in ['pod5', 'fast5']
        unaligned:   it[2] in ['bam', 'pacbio_bam', 'pacbio_hifi', 'hifi', 'ubam', 'unaligned_bam']
        aligned:     it[2] == 'aligned_bam'
        fastq:       it[2] in ['fastq', 'fastq.gz', 'fq', 'fq.gz']
    }

    // =========================================================================
    // STEP 3 (Raw ONT): Host GPU Basecalling & Demultiplexing
    // =========================================================================
    BASECALL_DEMUX(ch_branched.raw_signal)
    def basecalled_bams = BASECALL_DEMUX.out.demux_bams
    ch_versions = ch_versions.mix(BASECALL_DEMUX.out.versions)

    // Format and merge direct input unaligned BAMs
    def direct_bams_raw = ch_branched.unaligned.map { run_id, input_path, input_type, kit, model, modified_bases_model ->
        tuple(run_id, file(input_path))
    }
    MERGE_BAMS(direct_bams_raw)
    ch_versions = ch_versions.mix(MERGE_BAMS.out.versions)
    def direct_bams = MERGE_BAMS.out.bam

    def all_demux_bams = basecalled_bams.mix(direct_bams)

    // =========================================================================
    // STEP 2: Pre-Alignment Read QC (FastQC, Cramino, NanoPlot)
    // =========================================================================
    READ_QC(all_demux_bams)
    ch_versions = ch_versions.mix(READ_QC.out.versions)
    def fastq_ch = READ_QC.out.fastq.mix(
        ch_branched.fastq.map { run_id, input_path, input_type, kit, model, mod_model ->
            tuple(run_id, file(input_path))
        }
    )
    ch_qc_files = ch_qc_files
        .concat(READ_QC.out.fastqc_files)
        .concat(READ_QC.out.nanoplot_files)
        .concat(READ_QC.out.cramino_files)

    // =========================================================================
    // STEP 3 (De Novo): De Novo Assembly & Consensus Polishing
    // =========================================================================
    def assembly_fasta = Channel.empty()
    if (params.run_type == 'de_novo_assembly') {
        DE_NOVO_ASSEMBLY(fastq_ch)
        ch_versions = ch_versions.mix(DE_NOVO_ASSEMBLY.out.versions)
        assembly_fasta = DE_NOVO_ASSEMBLY.out.assembly
        ch_qc_files = ch_qc_files
            .concat(DE_NOVO_ASSEMBLY.out.quast_reports)
            .concat(DE_NOVO_ASSEMBLY.out.busco_reports)
    }

    // =========================================================================
    // REFERENCE-GUIDED BRANCH (Steps 4 - 11)
    // =========================================================================
    if (params.run_type == 'reference_based') {
        def ref_nuc_path = file(params.ref_nuc)
        def ref_fai_path = file("${params.ref_nuc}.fai")

        // Pre-aligned BAMs: verify index existence or auto-generate with SAMTOOLS_INDEX
        def ch_prealigned_check = ch_branched.aligned
            .map { run_id, input_path, input_type, kit, model, mod_model ->
                def bamFile = file(input_path)
                def baiFile = file("${input_path}.bai")
                def altBai  = file(input_path.replaceAll(/\.bam$/, '.bai'))
                def csiFile = file("${input_path}.csi")
                def hasIdx = baiFile.exists() || altBai.exists() || csiFile.exists()
                def resolvedBai = baiFile.exists() ? baiFile : (altBai.exists() ? altBai : (csiFile.exists() ? csiFile : null))
                tuple(run_id, bamFile, resolvedBai, hasIdx)
            }
            .branch {
                has_index:   it[3] == true
                needs_index: it[3] == false
            }

        def ch_prealigned_ready = ch_prealigned_check.has_index
            .map { run_id, bam, bai, hasIdx ->
                tuple(run_id, bam, bai, 'pre_aligned')
            }

        SAMTOOLS_INDEX(ch_prealigned_check.needs_index.map { run_id, bam, bai, hasIdx -> tuple(run_id, bam) })
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
        def ch_prealigned_indexed = SAMTOOLS_INDEX.out.bam_bai
            .map { run_id, bam, bai ->
                tuple(run_id, bam, bai, 'pre_aligned_indexed')
            }

        def ch_all_prealigned = ch_prealigned_ready.mix(ch_prealigned_indexed)

        // ---------------------------------------------------------------------
        // STEP 4: Reference Alignment & Kinetic Tag Preservation (pbmm2 / minimap2)
        // ---------------------------------------------------------------------
        def ch_newly_aligned = Channel.empty()
        if (!params.skip_alignment) {
            ALIGNMENT(fastq_ch, all_demux_bams, ref_nuc_path)
            ch_versions = ch_versions.mix(ALIGNMENT.out.versions)
            ch_newly_aligned = ALIGNMENT.out.aligned_bam_bai
                .map { run_id, bam, bai ->
                    tuple(run_id, bam, bai, 'pipeline_aligned')
                }
        }

        def ch_all_aligned_with_prov = ch_all_prealigned.mix(ch_newly_aligned)
        def aligned_bams = ch_all_aligned_with_prov.map { run_id, bam, bai, prov -> tuple(run_id, bam, bai) }

        // Post-alignment Coverage QC (Mosdepth)
        MOSDEPTH(aligned_bams, ref_nuc_path)
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
        ch_qc_files = ch_qc_files.concat(MOSDEPTH.out.global_txt.map { it[1] })
                                 .concat(MOSDEPTH.out.summary_txt.map { it[1] })

        // ---------------------------------------------------------------------
        // STEP 5: Small Variant Calling (SNVs & InDels via BCFtools / DeepVariant)
        // ---------------------------------------------------------------------
        def small_vcfs = Channel.empty()
        if (!params.skip_variants) {
            VARIANT_CALLING(aligned_bams, ref_nuc_path, ref_fai_path)
            ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)
            small_vcfs = VARIANT_CALLING.out.vcf_tbi
            ch_qc_files = ch_qc_files.concat(small_vcfs.map { it[1] })
        }

        // ---------------------------------------------------------------------
        // STEP 6: Structural Variants, CNVs & Repeat Expansions (Sniffles2 / TRGT / HiFiCNV)
        // ---------------------------------------------------------------------
        def sv_vcfs = Channel.empty()
        def sv_vcfs_tbi = Channel.empty()
        def trgt_vcfs_tbi = Channel.empty()
        if (!params.skip_variants) {
            STRUCTURAL_VARIANTS(aligned_bams, ref_nuc_path)
            ch_versions   = ch_versions.mix(STRUCTURAL_VARIANTS.out.versions)
            sv_vcfs       = STRUCTURAL_VARIANTS.out.sv_vcf
            sv_vcfs_tbi   = STRUCTURAL_VARIANTS.out.sv_vcf_tbi
            trgt_vcfs_tbi = STRUCTURAL_VARIANTS.out.trgt_vcf
        }

        // ---------------------------------------------------------------------
        // STEP 7: Epigenomics & Base Modifications (modkit pileup / pbcpgtools)
        // ---------------------------------------------------------------------
        def bed_methylation = Channel.empty()
        if (!params.skip_methylation) {
            METHYLATION(aligned_bams, ref_nuc_path)
            ch_versions     = ch_versions.mix(METHYLATION.out.versions)
            bed_methylation = METHYLATION.out.bed_methylation
            ch_qc_files     = ch_qc_files.mix(METHYLATION.out.qc_report.map { id, txt -> txt })
        }

        // ---------------------------------------------------------------------
        // STEP 8: Multi-Variant Joint Co-Phasing (HiPhase / WhatsHap)
        // ---------------------------------------------------------------------
        def phased_bams = aligned_bams
        def phased_vcfs = small_vcfs
        if (!params.skip_variants && (params.species == 'human' || params.enable_phasing)) {
            PHASING(aligned_bams, small_vcfs, sv_vcfs_tbi, trgt_vcfs_tbi, ref_nuc_path, ref_fai_path)
            ch_versions     = ch_versions.mix(PHASING.out.versions)
            phased_bams     = PHASING.out.phased_bam_bai
            phased_vcfs     = PHASING.out.phased_vcf
            ch_qc_files     = ch_qc_files.mix(PHASING.out.phasing_stats.map { id, stats -> stats })
        }

        // ---------------------------------------------------------------------
        // STEP 9: Functional & Clinical Variant Annotation (SnpEff / VEP)
        // ---------------------------------------------------------------------
        def final_vcfs = phased_vcfs
        if (!params.skip_variants) {
            ANNOTATION(phased_vcfs)
            ch_versions = ch_versions.mix(ANNOTATION.out.versions)
            final_vcfs = ANNOTATION.out.annotated_vcf
        }

        // ---------------------------------------------------------------------
        // STEP 10: Empirical Truth-Set Benchmarking (GIAB HG002)
        // ---------------------------------------------------------------------
        if (params.run_benchmarking) {
            def truth_snv_vcf = params.benchmark_truth_vcf ? file(params.benchmark_truth_vcf) : null
            def truth_snv_bed = params.benchmark_truth_bed ? file(params.benchmark_truth_bed) : null
            def truth_sv_vcf  = params.benchmark_sv_truth_vcf ? file(params.benchmark_sv_truth_vcf) : null
            def truth_sv_bed  = params.benchmark_sv_truth_bed ? file(params.benchmark_sv_truth_bed) : null

            BENCHMARKING(
                phased_vcfs,
                sv_vcfs_tbi,
                ref_nuc_path,
                ref_fai_path,
                truth_snv_vcf,
                truth_snv_bed,
                truth_sv_vcf,
                truth_sv_bed
            )
            ch_versions = ch_versions.mix(BENCHMARKING.out.versions)
            ch_qc_files = ch_qc_files.mix(BENCHMARKING.out.qc_files)
        }

        // ---------------------------------------------------------------------
        // STEP 11: Clinical Governance & Forensic Audit Trail (milou pattern)
        // ---------------------------------------------------------------------
        def ch_sanitized_bams = phased_bams
        def ch_sanitized_vcfs = final_vcfs
        if (params.enable_audit_trail || params.enable_phi_sanitization) {
            def extra_artifacts = bed_methylation.map { id, bed -> tuple(id, [bed]) }
            CLINICAL_GOVERNANCE(phased_bams, final_vcfs, ch_input_checksums, ref_nuc_path, extra_artifacts)
            ch_versions       = ch_versions.mix(CLINICAL_GOVERNANCE.out.versions)
            ch_sanitized_bams = CLINICAL_GOVERNANCE.out.sanitized_bams
            ch_sanitized_vcfs = CLINICAL_GOVERNANCE.out.sanitized_vcfs
        }

        // ---------------------------------------------------------------------
        // STEP 11: Interactive Web Visualization (JBrowse 2, Methylartist & Methylation Dashboard)
        // ---------------------------------------------------------------------
        def ref_tuple_ch = fastq_ch.map { sample_id, fastq ->
            tuple(sample_id, [ref_nuc_path, ref_fai_path])
        }
        VISUALIZATION(ref_tuple_ch, ch_sanitized_bams, ch_sanitized_vcfs, bed_methylation)
        ch_versions = ch_versions.mix(VISUALIZATION.out.versions)
        ch_qc_files = ch_qc_files.mix(VISUALIZATION.out.methylation_html.map { meta, html -> html })
    }

    // =========================================================================
    // STEP 12: Reporting Dashboard (MultiQC & Quarto)
    // =========================================================================
    REPORTING(ch_qc_files, ch_versions)
}
