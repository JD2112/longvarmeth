// subworkflows/local/alignment/main.nf
// Step 4: Reference alignment and kinetic tag streaming (pbmm2 / minimap2 -y)

include { MINIMAP2 as ALIGN_NUC } from '../../../modules/nf-core/minimap2'
include { PBMM2_ALIGN } from '../../../modules/nf-core/pbmm2/align'
include { SAMTOOLS_SORT as SORT_NUC } from '../../../modules/nf-core/samtools/sort'

workflow ALIGNMENT {
    take:
    fastq_ch      // tuple(sample_id, fastq)
    demux_bams    // tuple(sample_id, bam)
    ref_nuc_path  // path(reference_fasta)

    main:
    ch_versions = Channel.empty()
    def ch_aligned_bam_bai = Channel.empty()

    if (params.platform == 'pacbio') {
        // Native PacBio HiFi Alignment (pbmm2 preserves kinetic tags MM/ML)
        PBMM2_ALIGN(demux_bams, ref_nuc_path)
        ch_versions = ch_versions.mix(PBMM2_ALIGN.out.versions)
        ch_aligned_bam_bai = PBMM2_ALIGN.out.bam_bai
    } else {
        // Oxford Nanopore Alignment with Native Tag Streaming (-y preserves MM/ML)
        // Prefer unaligned BAM to stream all base-modification tags directly into alignment
        def align_input = demux_bams.mix(fastq_ch)
            .groupTuple(by: 0)
            .map { sample_id, files ->
                def bamFile = files.find { it.name.endsWith('.bam') }
                tuple(sample_id, bamFile ?: files[0])
            }

        ALIGN_NUC(align_input, ref_nuc_path)
        ch_versions = ch_versions.mix(ALIGN_NUC.out.versions)

        SORT_NUC(ALIGN_NUC.out.sam)
        ch_versions = ch_versions.mix(SORT_NUC.out.versions)
        ch_aligned_bam_bai = SORT_NUC.out.bam_bai
    }

    emit:
    aligned_bam_bai = ch_aligned_bam_bai
    versions        = ch_versions
}
