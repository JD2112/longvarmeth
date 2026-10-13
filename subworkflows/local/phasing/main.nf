// subworkflows/local/phasing/main.nf
// Step 8: Haplotype phasing and read haplotagging (WhatsHap / HiPhase)

include { WHATSHAP_PHASE } from '../../../modules/nf-core/whatshap/phase'
include { WHATSHAP_HAPLOTAG } from '../../../modules/nf-core/whatshap/haplotag'
include { HIPHASE } from '../../../modules/nf-core/hiphase'

workflow PHASING {
    take:
    aligned_bam_bai // tuple(sample_id, bam, bai)
    small_vcf_tbi   // tuple(sample_id, vcf, tbi)
    sv_vcf_tbi      // tuple(sample_id, vcf, tbi)
    ref_nuc_path    // path(reference_fasta)
    ref_fai_path    // path(reference_fai)

    main:
    ch_versions = Channel.empty()
    def ch_phased_bam_bai = Channel.empty()
    def ch_phased_vcf = Channel.empty()

    if (params.platform == 'pacbio' || params.phasing_tool == 'hiphase') {
        // HiPhase Co-Phasing (SNVs + SVs + STRs)
        def hiphase_input = aligned_bam_bai
            .join(small_vcf_tbi)
            .map { sample_id, bam, bai, vcf, tbi ->
                tuple(sample_id, bam, bai, vcf, tbi, [], [])
            }
        HIPHASE(hiphase_input, ref_nuc_path)
        ch_versions = ch_versions.mix(HIPHASE.out.versions)
        ch_phased_bam_bai = HIPHASE.out.phased_bam_bai
        ch_phased_vcf     = HIPHASE.out.small_vcf_phased
    } else {
        // WhatsHap Phasing & Read Haplotagging
        def ref_with_fai = [ref_nuc_path, ref_fai_path]
        def phase_input = small_vcf_tbi
            .join(aligned_bam_bai)
            .map { sample_id, vcf, tbi, bam, bai ->
                tuple(sample_id, vcf, bam, bai)
            }
        WHATSHAP_PHASE(phase_input, ref_with_fai)
        ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions)
        ch_phased_vcf = WHATSHAP_PHASE.out.phased_vcf

        def haplotag_input = WHATSHAP_PHASE.out.phased_vcf
            .join(aligned_bam_bai)
            .map { sample_id, vcf, tbi, bam, bai ->
                tuple(sample_id, vcf, tbi, bam, bai)
            }
        WHATSHAP_HAPLOTAG(haplotag_input, ref_with_fai)
        ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)
        ch_phased_bam_bai = WHATSHAP_HAPLOTAG.out.haplotagged_bam_bai
    }

    emit:
    phased_bam_bai = ch_phased_bam_bai
    phased_vcf     = ch_phased_vcf
    versions       = ch_versions
}
