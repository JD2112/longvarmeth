// subworkflows/local/structural_variants/main.nf
// Step 6: Structural variants (SVs), copy number alterations (CNVs), and tandem repeats (STRs)

include { SNIFFLES } from '../../../modules/nf-core/sniffles'
include { CNVKIT } from '../../../modules/nf-core/cnvkit'
include { HIFINCNV } from '../../../modules/nf-core/hificnv'
include { TRGT } from '../../../modules/nf-core/trgt'

workflow STRUCTURAL_VARIANTS {
    take:
    aligned_bam_bai // tuple(sample_id, bam, bai)
    ref_nuc_path    // path(reference_fasta)

    main:
    ch_versions = Channel.empty()
    def ch_sv_vcf = Channel.empty()
    def ch_cnv_calls = Channel.empty()
    def ch_trgt_vcf = Channel.empty()

    // 1. Structural Variants (Sniffles2)
    SNIFFLES(aligned_bam_bai, ref_nuc_path)
    ch_versions = ch_versions.mix(SNIFFLES.out.versions)
    ch_sv_vcf = SNIFFLES.out.sv_vcf


    // 2. Copy Number Alterations (CNVkit / HiFiCNV)
    if (params.platform == 'pacbio' && params.cnv_caller == 'hificnv') {
        HIFINCNV(aligned_bam_bai, ref_nuc_path, [])
        ch_versions = ch_versions.mix(HIFINCNV.out.versions)
        ch_cnv_calls = HIFINCNV.out.vcf_tbi
    } else {
        CNVKIT(aligned_bam_bai, ref_nuc_path)
        ch_versions = ch_versions.mix(CNVKIT.out.versions)
        ch_cnv_calls = CNVKIT.out.cnv_files
    }


    // 3. Short Tandem Repeats (TRGT)
    if (params.enable_trgt && params.repeat_bed) {
        def repeat_bed_file = file(params.repeat_bed)
        TRGT(aligned_bam_bai, ref_nuc_path, repeat_bed_file)
        ch_versions = ch_versions.mix(TRGT.out.versions)
        ch_trgt_vcf = TRGT.out.vcf_tbi
    }

    emit:
    sv_vcf   = ch_sv_vcf
    cnv_calls= ch_cnv_calls
    trgt_vcf = ch_trgt_vcf
    versions = ch_versions
}
