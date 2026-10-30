// subworkflows/local/variant_calling/main.nf
// Step 5: Small variant calling (SNVs and InDels) with modular caller dispatch

include { BCFTOOLS_CALL as CALL_NUC } from '../../../modules/nf-core/bcftools/call'
include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat'
include { DEEPVARIANT } from '../../../modules/nf-core/deepvariant'
include { CLAIR3 } from '../../../modules/nf-core/clair3'

workflow VARIANT_CALLING {
    take:
    aligned_bam_bai // tuple(sample_id, bam, bai)
    ref_nuc_path    // path(reference_fasta)
    ref_fai_path    // path(reference_fai)

    main:
    ch_versions = Channel.empty()
    def ch_vcf_tbi = Channel.empty()

    if (params.snv_caller == 'deepvariant') {
        // High-precision DeepVariant calling with chemistry-specific neural network model
        def model_type = params.platform == 'pacbio' ? 'PACBIO' : 'ONT_R104'
        DEEPVARIANT(aligned_bam_bai, ref_nuc_path, ref_fai_path, model_type)
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
        ch_vcf_tbi  = DEEPVARIANT.out.vcf_tbi
    } else if (params.snv_caller == 'clair3') {
        // Clair3 neural network long-read variant caller
        def clair3_plat = params.platform == 'pacbio' ? 'hifi' : 'ont'
        CLAIR3(aligned_bam_bai, ref_nuc_path, ref_fai_path, clair3_plat)
        ch_versions = ch_versions.mix(CLAIR3.out.versions)
        ch_vcf_tbi  = CLAIR3.out.vcf_tbi
    } else {
        // Fallback: Parallel chromosome calling with BCFtools mpileup/call
        def ch_chroms = Channel.fromPath(ref_fai_path)
            .splitCsv(sep: '\t')
            .map { row -> row[0] }

        def call_input = aligned_bam_bai.combine(ch_chroms)

        CALL_NUC(call_input, ref_nuc_path)
        ch_versions = ch_versions.mix(CALL_NUC.out.versions)

        def concat_input = CALL_NUC.out.vcf_tbi
            .groupTuple(by: 0)
            .map { sample_id, vcfs, tbis ->
                def sorted = [vcfs, tbis].transpose().sort { a, b -> a[0].name <=> b[0].name }
                tuple(sample_id, sorted.collect { it[0] }, sorted.collect { it[1] })
            }

        BCFTOOLS_CONCAT(concat_input)
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        ch_vcf_tbi  = BCFTOOLS_CONCAT.out.vcf_tbi
    }

    emit:
    vcf_tbi  = ch_vcf_tbi
    versions = ch_versions
}
