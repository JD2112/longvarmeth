// subworkflows/local/annotation_vis/main.nf
include { METHYLARTIST_LOCUS } from '../../../modules/local/methylartist'
include { JBROWSE_CREATE } from '../../../modules/local/jbrowse'

workflow ANNOTATION_VIS {
    take:
    assembly_ch // tuple(sample_id, assembly_fasta)
    bam_bai_ch // tuple(sample_id, bam, bai)
    vcf_csi_ch // tuple(sample_id, vcf, csi)

    main:
    ch_versions = Channel.empty()

    if (params.methylartist_locus && params.ref_nuc) {
        def ref_nuc_path = file(params.ref_nuc)
        METHYLARTIST_LOCUS(bam_bai_ch, ref_nuc_path, params.methylartist_locus)
        ch_versions = ch_versions.mix(METHYLARTIST_LOCUS.out.versions)
    }

    def jbrowse_input = assembly_ch
        .join(bam_bai_ch)
        .join(vcf_csi_ch)
        .map { sample_id, assembly, bam, bai, vcf, csi ->
            tuple(sample_id, assembly, bam, bai, vcf, csi)
        }

    JBROWSE_CREATE(jbrowse_input)
    ch_versions = ch_versions.mix(JBROWSE_CREATE.out.versions)

    emit:
    jbrowse_dir = JBROWSE_CREATE.out.jbrowse_dir
    versions    = ch_versions
}
