// subworkflows/local/variant_calling/main.nf
// Step 5: Small variant calling (SNVs and InDels)

include { BCFTOOLS_CALL as CALL_NUC } from '../../../modules/nf-core/bcftools/call'
include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat'

workflow VARIANT_CALLING {
    take:
    aligned_bam_bai // tuple(sample_id, bam, bai)
    ref_nuc_path    // path(reference_fasta)
    ref_fai_path    // path(reference_fai)

    main:
    ch_versions = Channel.empty()

    // 1. Generate chromosome calling partitions
    def ch_chroms = Channel.fromPath(ref_fai_path)
        .splitCsv(sep: '\t')
        .map { row -> row[0] }
        .filter { chrom -> 
            if (params.species == 'human') {
                return chrom =~ /^chr[0-9XYM]+$/ || chrom =~ /^[0-9XYM]+$/
            }
            return true
        }

    def call_input = aligned_bam_bai.combine(ch_chroms)

    // 2. Parallel small variant calling
    CALL_NUC(call_input, ref_nuc_path)
    ch_versions = ch_versions.mix(CALL_NUC.out.versions)

    // 3. Genome-wide concatenation & indexing
    def concat_input = CALL_NUC.out.vcf_tbi
        .groupTuple(by: 0)
        .map { sample_id, vcfs, tbis ->
            def sorted = [vcfs, tbis].transpose().sort { a, b -> a[0].name <=> b[0].name }
            tuple(sample_id, sorted.collect { it[0] }, sorted.collect { it[1] })
        }

    BCFTOOLS_CONCAT(concat_input)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    emit:
    vcf_tbi  = BCFTOOLS_CONCAT.out.vcf_tbi
    versions = ch_versions
}
