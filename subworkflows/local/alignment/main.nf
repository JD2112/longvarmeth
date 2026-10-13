// subworkflows/local/alignment/main.nf
// Step 4: Reference alignment, kinetic tag preservation, and post-alignment QC

include { MINIMAP2 as ALIGN_NUC } from '../../../modules/nf-core/minimap2'
include { PBMM2_ALIGN } from '../../../modules/nf-core/pbmm2/align'
include { SAMTOOLS_SORT as SORT_NUC } from '../../../modules/nf-core/samtools/sort'
include { TRANSFER_TAGS as TAG_NUC } from '../../../modules/local/transfer_tags'

workflow ALIGNMENT {
    take:
    fastq_ch      // tuple(sample_id, fastq)
    demux_bams    // tuple(sample_id, bam)
    ref_nuc_path  // path(reference_fasta)

    main:
    ch_versions = Channel.empty()
    def ch_aligned_bam_bai = Channel.empty()

    if (params.platform == 'pacbio') {
        // Native PacBio HiFi Alignment (preserves kinetic tags MM/ML)
        PBMM2_ALIGN(demux_bams, ref_nuc_path)
        ch_versions = ch_versions.mix(PBMM2_ALIGN.out.versions)
        ch_aligned_bam_bai = PBMM2_ALIGN.out.bam_bai
    } else {
        // Oxford Nanopore Alignment + MM/ML Tag Restoration
        ALIGN_NUC(fastq_ch, ref_nuc_path)
        ch_versions = ch_versions.mix(ALIGN_NUC.out.versions)

        SORT_NUC(ALIGN_NUC.out.sam)
        ch_versions = ch_versions.mix(SORT_NUC.out.versions)
        
        def tag_nuc_input = SORT_NUC.out.bam_bai
            .join(demux_bams)
            .map { sample_id, sorted_bam, sorted_bai, raw_bam ->
                tuple(sample_id, sorted_bam, sorted_bai, raw_bam)
            }
        TAG_NUC(tag_nuc_input)
        ch_versions = ch_versions.mix(TAG_NUC.out.versions)
        ch_aligned_bam_bai = TAG_NUC.out.tagged_bam_bai
    }

    emit:
    aligned_bam_bai = ch_aligned_bam_bai
    versions        = ch_versions
}
