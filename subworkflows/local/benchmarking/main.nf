// subworkflows/local/benchmarking/main.nf
// Benchmarking against truth sets (GIAB HG002) using hap.py (SNVs/indels) and Truvari (SVs)

include { HAPPY } from '../../../modules/nf-core/happy'
include { TRUVARI } from '../../../modules/nf-core/truvari'

workflow BENCHMARKING {
    take:
    small_vcf_tbi       // tuple(sample_id, vcf, tbi)
    sv_vcf_tbi          // tuple(sample_id, vcf, tbi)
    ref_nuc_path        // path(ref_fasta)
    ref_fai_path        // path(ref_fai)
    truth_snv_vcf       // path(truth_snv_vcf) or null
    truth_snv_bed       // path(truth_snv_bed) or null
    truth_sv_vcf        // path(truth_sv_vcf) or null
    truth_sv_bed        // path(truth_sv_bed) or null

    main:
    ch_versions = Channel.empty()
    ch_qc_files = Channel.empty()
    ch_happy_summary = Channel.empty()
    ch_truvari_summary = Channel.empty()

    // 1. SNV / Indel Benchmarking with hap.py
    if (truth_snv_vcf && truth_snv_bed) {
        def happy_input = small_vcf_tbi.map { id, vcf, tbi ->
            def meta = (id instanceof Map) ? id : [id: id]
            tuple(meta, vcf, tbi)
        }
        def truth_snv_tbi = file("${truth_snv_vcf}.tbi").exists() ? file("${truth_snv_vcf}.tbi") : file("${truth_snv_vcf}.csi")
        HAPPY(
            happy_input,
            truth_snv_vcf,
            truth_snv_tbi,
            truth_snv_bed,
            ref_nuc_path,
            ref_fai_path
        )
        ch_versions      = ch_versions.mix(HAPPY.out.versions)
        ch_happy_summary = HAPPY.out.summary
        ch_qc_files      = ch_qc_files.mix(HAPPY.out.mqc.map { meta, tsv -> tsv })
    }

    // 2. Structural Variant Benchmarking with Truvari
    if (truth_sv_vcf && truth_sv_bed) {
        def truvari_input = sv_vcf_tbi.map { id, vcf, tbi ->
            def meta = (id instanceof Map) ? id : [id: id]
            tuple(meta, vcf, tbi)
        }
        def truth_sv_tbi = file("${truth_sv_vcf}.tbi").exists() ? file("${truth_sv_vcf}.tbi") : file("${truth_sv_vcf}.csi")
        TRUVARI(
            truvari_input,
            truth_sv_vcf,
            truth_sv_tbi,
            truth_sv_bed,
            ref_nuc_path,
            ref_fai_path
        )
        ch_versions        = ch_versions.mix(TRUVARI.out.versions)
        ch_truvari_summary = TRUVARI.out.summary
        ch_qc_files        = ch_qc_files.mix(TRUVARI.out.mqc.map { meta, tsv -> tsv })
    }

    emit:
    happy_summary   = ch_happy_summary
    truvari_summary = ch_truvari_summary
    qc_files        = ch_qc_files
    versions        = ch_versions
}

