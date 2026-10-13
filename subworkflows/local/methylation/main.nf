// subworkflows/local/methylation/main.nf
// Step 7: Epigenomics and base modifications (modkit for ONT; pb-cpg-tools + methbat for PacBio)

include { MODKIT_PILEUP as MOD_NUC } from '../../../modules/nf-core/modkit/pileup'

workflow METHYLATION {
    take:
    aligned_bam_bai // tuple(sample_id, bam, bai)
    ref_nuc_path    // path(reference_fasta)

    main:
    ch_versions = Channel.empty()
    def ch_bed_methylation = Channel.empty()

    // 1. Modkit Whole-Genome 5mC / 5hmC / 6mA Pileup
    MOD_NUC(aligned_bam_bai, ref_nuc_path)
    ch_versions = ch_versions.mix(MOD_NUC.out.versions)
    ch_bed_methylation = MOD_NUC.out.bed

    emit:
    bed_methylation = ch_bed_methylation
    versions        = ch_versions
}
