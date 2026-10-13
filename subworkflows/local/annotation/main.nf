// subworkflows/local/annotation/main.nf
// Step 9: Functional and clinical variant effect annotation (SnpEff / VEP)

include { SNPEFF_DOWNLOAD } from '../../../modules/nf-core/snpeff/download'
include { SNPEFF } from '../../../modules/nf-core/snpeff/annotate'
include { BCFTOOLS_CONCAT as COMPRESS_INDEX_SNPEFF } from '../../../modules/nf-core/bcftools/concat'
include { VEP_DOWNLOAD } from '../../../modules/nf-core/vep/download'
include { VEP } from '../../../modules/nf-core/vep/annotate'

workflow ANNOTATION {
    take:
    vcf_tbi // tuple(sample_id, vcf, tbi)

    main:
    ch_versions = Channel.empty()
    def ch_annotated_vcf = vcf_tbi

    if (params.enable_snpeff) {
        SNPEFF_DOWNLOAD()
        ch_versions = ch_versions.mix(SNPEFF_DOWNLOAD.out.versions)

        SNPEFF(ch_annotated_vcf, SNPEFF_DOWNLOAD.out.ready)
        ch_versions = ch_versions.mix(SNPEFF.out.versions)

        COMPRESS_INDEX_SNPEFF(SNPEFF.out.vcf)
        ch_annotated_vcf = COMPRESS_INDEX_SNPEFF.out.vcf_tbi
    }

    if (params.enable_vep) {
        VEP_DOWNLOAD()
        ch_versions = ch_versions.mix(VEP_DOWNLOAD.out.versions)

        VEP(ch_annotated_vcf, VEP_DOWNLOAD.out.ready)
        ch_versions = ch_versions.mix(VEP.out.versions)
        ch_annotated_vcf = VEP.out.vcf_tbi
    }

    emit:
    annotated_vcf = ch_annotated_vcf
    versions      = ch_versions
}
