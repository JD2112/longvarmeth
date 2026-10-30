// subworkflows/local/phasing/main.nf
// Step 8: Multi-variant joint co-phasing and read haplotagging (WhatsHap / HiPhase)

include { WHATSHAP_PHASE } from '../../../modules/nf-core/whatshap/phase'
include { WHATSHAP_HAPLOTAG } from '../../../modules/nf-core/whatshap/haplotag'
include { WHATSHAP_STATS } from '../../../modules/nf-core/whatshap/stats'
include { HIPHASE } from '../../../modules/nf-core/hiphase'

workflow PHASING {
    take:
    aligned_bam_bai // tuple(sample_id, bam, bai)
    small_vcf_tbi   // tuple(sample_id, vcf, tbi)
    sv_vcf_tbi      // tuple(sample_id, vcf, tbi)
    trgt_vcf_tbi    // tuple(sample_id, vcf, tbi)
    ref_nuc_path    // path(reference_fasta)
    ref_fai_path    // path(reference_fai)

    main:
    ch_versions = Channel.empty()
    def ch_phased_bam_bai = Channel.empty()
    def ch_phased_vcf = Channel.empty()
    def ch_phased_sv_vcf = Channel.empty()
    def ch_phasing_stats = Channel.empty()

    if (params.platform == 'pacbio' || params.phasing_tool == 'hiphase') {
        // HiPhase Joint Multi-Variant Co-Phasing (SNVs + SVs + STRs)
        def hiphase_input = aligned_bam_bai
            .join(small_vcf_tbi)
            .join(sv_vcf_tbi, remainder: true)
            .join(trgt_vcf_tbi, remainder: true)
            .map { row ->
                def sample_id = row[0]
                def bam       = row[1]
                def bai       = row[2]
                def vcf       = row[3]
                def tbi       = row[4]
                def sv_vcf    = (row.size() > 5 && row[5] != null) ? row[5] : []
                def sv_tbi    = (row.size() > 6 && row[6] != null) ? row[6] : []
                def str_vcf   = (row.size() > 7 && row[7] != null) ? row[7] : []
                def str_tbi   = (row.size() > 8 && row[8] != null) ? row[8] : []
                tuple(sample_id, bam, bai, vcf, tbi, sv_vcf, sv_tbi, str_vcf, str_tbi)
            }

        HIPHASE(hiphase_input, ref_nuc_path)
        ch_versions       = ch_versions.mix(HIPHASE.out.versions)
        ch_phased_bam_bai = HIPHASE.out.phased_bam_bai
        ch_phased_vcf     = HIPHASE.out.small_vcf_phased
        ch_phased_sv_vcf  = HIPHASE.out.sv_vcf_phased
        ch_phasing_stats  = HIPHASE.out.stats_csv
    } else {
        // WhatsHap Phasing & Read Haplotagging
        def ref_with_fai = [ref_nuc_path, ref_fai_path]
        def phase_input = small_vcf_tbi
            .join(aligned_bam_bai)
            .map { sample_id, vcf, tbi, bam, bai ->
                tuple(sample_id, vcf, bam, bai)
            }
        WHATSHAP_PHASE(phase_input, ref_with_fai)
        ch_versions   = ch_versions.mix(WHATSHAP_PHASE.out.versions)
        ch_phased_vcf = WHATSHAP_PHASE.out.phased_vcf

        def haplotag_input = WHATSHAP_PHASE.out.phased_vcf
            .join(aligned_bam_bai)
            .map { sample_id, vcf, tbi, bam, bai ->
                tuple(sample_id, vcf, tbi, bam, bai)
            }
        WHATSHAP_HAPLOTAG(haplotag_input, ref_with_fai)
        ch_versions       = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)
        ch_phased_bam_bai = WHATSHAP_HAPLOTAG.out.haplotagged_bam_bai

        // Phasing metrics and block statistics
        WHATSHAP_STATS(WHATSHAP_PHASE.out.phased_vcf)
        ch_versions      = ch_versions.mix(WHATSHAP_STATS.out.versions)
        ch_phasing_stats = WHATSHAP_STATS.out.stats_tsv
    }

    emit:
    phased_bam_bai = ch_phased_bam_bai
    phased_vcf     = ch_phased_vcf
    phased_sv_vcf  = ch_phased_sv_vcf
    phasing_stats  = ch_phasing_stats
    versions       = ch_versions
}
