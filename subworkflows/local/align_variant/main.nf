// subworkflows/local/align_variant/main.nf
include { MINIMAP2 as ALIGN_NUC; MINIMAP2 as ALIGN_MT_AP; MINIMAP2 as ALIGN_MT_AU } from '../../../modules/nf-core/minimap2'
include { SAMTOOLS_SORT as SORT_NUC; SAMTOOLS_SORT as SORT_MT_AP; SAMTOOLS_SORT as SORT_MT_AU } from '../../../modules/nf-core/samtools/sort'
include { TRANSFER_TAGS as TAG_NUC; TRANSFER_TAGS as TAG_MT_AP; TRANSFER_TAGS as TAG_MT_AU } from '../../../modules/local/transfer_tags'
include { BCFTOOLS_CALL as CALL_NUC; BCFTOOLS_CALL as CALL_MT_AP; BCFTOOLS_CALL as CALL_MT_AU } from '../../../modules/nf-core/bcftools/call'
include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat'
include { SNIFFLES } from '../../../modules/nf-core/sniffles'
include { CNVKIT } from '../../../modules/nf-core/cnvkit'
include { WHATSHAP_PHASE } from '../../../modules/nf-core/whatshap/phase'
include { WHATSHAP_HAPLOTAG } from '../../../modules/nf-core/whatshap/haplotag'
include { MODKIT_PILEUP as MOD_NUC; MODKIT_PILEUP as MOD_MT_AP; MODKIT_PILEUP as MOD_MT_AU } from '../../../modules/nf-core/modkit/pileup'
include { SNPEFF_DOWNLOAD } from '../../../modules/nf-core/snpeff/download'
include { SNPEFF } from '../../../modules/nf-core/snpeff/annotate'
include { BCFTOOLS_CONCAT as COMPRESS_INDEX_SNPEFF } from '../../../modules/nf-core/bcftools/concat'
include { VEP_DOWNLOAD } from '../../../modules/nf-core/vep/download'
include { VEP } from '../../../modules/nf-core/vep/annotate'

workflow ALIGN_VARIANT {
    take:
    fastq_ch // tuple(sample_id, fastq)
    demux_bams // tuple(sample_id, bam)

    main:
    ch_versions = Channel.empty()

    def ch_tagged_nuc_bam_bai = Channel.empty()
    def ch_tagged_mt_ap_bam_bai = Channel.empty()
    def ch_tagged_mt_au_bam_bai = Channel.empty()

    def ch_vcf_nuc = Channel.empty()
    def ch_vcf_mt_ap = Channel.empty()
    def ch_vcf_mt_au = Channel.empty()

    def ch_bed_nuc = Channel.empty()
    def ch_bed_mt_ap = Channel.empty()
    def ch_bed_mt_au = Channel.empty()

    def ch_phased_bam_bai = Channel.empty()

    // 1. Nuclear Alignment & Tag Transfer
    if (params.ref_nuc) {
        def ref_nuc_path = file(params.ref_nuc)
        def ref_nuc_with_fai = [file(params.ref_nuc), file("${params.ref_nuc}.fai")]
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
        ch_tagged_nuc_bam_bai = TAG_NUC.out.tagged_bam_bai

        // Variant Calling - Nuclear
        def ch_chroms = Channel.fromPath("${params.ref_nuc}.fai")
            .splitCsv(sep: '\t')
            .map { row -> row[0] }
            .filter { chrom -> 
                if (params.species == 'human') {
                    return chrom =~ /^chr[0-9XYM]+$/ || chrom =~ /^[0-9XYM]+$/
                }
                return true
            }

        def call_nuc_input = ch_tagged_nuc_bam_bai
            .combine(ch_chroms)

        CALL_NUC(call_nuc_input, ref_nuc_path)
        ch_versions = ch_versions.mix(CALL_NUC.out.versions)

        def concat_input = CALL_NUC.out.vcf_tbi
            .groupTuple(by: 0)
            .map { sample_id, vcfs, tbis ->
                def sorted = [vcfs, tbis].transpose().sort { a, b -> a[0].name <=> b[0].name }
                tuple(sample_id, sorted.collect { it[0] }, sorted.collect { it[1] })
            }

        BCFTOOLS_CONCAT(concat_input)
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        ch_vcf_nuc = BCFTOOLS_CONCAT.out.vcf_tbi

        // Annotation
        if (params.enable_snpeff) {
            SNPEFF_DOWNLOAD()
            ch_versions = ch_versions.mix(SNPEFF_DOWNLOAD.out.versions)

            SNPEFF(ch_vcf_nuc, SNPEFF_DOWNLOAD.out.ready)
            ch_versions = ch_versions.mix(SNPEFF.out.versions)

            COMPRESS_INDEX_SNPEFF(SNPEFF.out.vcf)
            ch_vcf_nuc = COMPRESS_INDEX_SNPEFF.out.vcf_tbi
        }
        if (params.enable_vep) {
            VEP_DOWNLOAD()
            ch_versions = ch_versions.mix(VEP_DOWNLOAD.out.versions)

            VEP(ch_vcf_nuc, VEP_DOWNLOAD.out.ready)
            ch_versions = ch_versions.mix(VEP.out.versions)
            ch_vcf_nuc = VEP.out.vcf_tbi
        }

        // Methylation - Nuclear
        MOD_NUC(ch_tagged_nuc_bam_bai, ref_nuc_path)
        ch_versions = ch_versions.mix(MOD_NUC.out.versions)
        ch_bed_nuc = MOD_NUC.out.bed

        // Phasing & Haplotagging
        if (params.species == 'human' || params.enable_phasing) {
            def phase_input = ch_vcf_nuc
                .join(ch_tagged_nuc_bam_bai)
                .map { sample_id, vcf, tbi, bam, bai ->
                    tuple(sample_id, vcf, bam, bai)
                }
            WHATSHAP_PHASE(phase_input, ref_nuc_with_fai)
            ch_versions = ch_versions.mix(WHATSHAP_PHASE.out.versions)

            def haplotag_input = WHATSHAP_PHASE.out.phased_vcf
                .join(ch_tagged_nuc_bam_bai)
                .map { sample_id, vcf, tbi, bam, bai ->
                    tuple(sample_id, vcf, tbi, bam, bai)
                }
            WHATSHAP_HAPLOTAG(haplotag_input, ref_nuc_with_fai)
            ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)
            ch_phased_bam_bai = WHATSHAP_HAPLOTAG.out.haplotagged_bam_bai
        }

        // Structural Variants & CNV
        if (params.species == 'eukaryote' || params.species == 'human') {
            SNIFFLES(ch_tagged_nuc_bam_bai, ref_nuc_path)
            ch_versions = ch_versions.mix(SNIFFLES.out.versions)

            CNVKIT(ch_tagged_nuc_bam_bai, ref_nuc_path)
            ch_versions = ch_versions.mix(CNVKIT.out.versions)
        }
    }

    // 2. Mitochondrial AP Alignment & Tag Transfer
    if (params.ref_mt_ap) {
        def ref_mt_ap_path = file(params.ref_mt_ap)
        ALIGN_MT_AP(fastq_ch, ref_mt_ap_path)
        SORT_MT_AP(ALIGN_MT_AP.out.sam)

        def tag_mt_ap_input = SORT_MT_AP.out.bam_bai
            .join(demux_bams)
            .map { sample_id, sorted_bam, sorted_bai, raw_bam ->
                tuple(sample_id, sorted_bam, sorted_bai, raw_bam)
            }
        TAG_MT_AP(tag_mt_ap_input)
        ch_tagged_mt_ap_bam_bai = TAG_MT_AP.out.tagged_bam_bai

        CALL_MT_AP(ch_tagged_mt_ap_bam_bai, ref_mt_ap_path)
        ch_vcf_mt_ap = CALL_MT_AP.out.vcf_tbi

        MOD_MT_AP(ch_tagged_mt_ap_bam_bai, ref_mt_ap_path)
        ch_bed_mt_ap = MOD_MT_AP.out.bed
    }

    // 3. Mitochondrial AU Alignment & Tag Transfer
    if (params.ref_mt_au) {
        def ref_mt_au_path = file(params.ref_mt_au)
        ALIGN_MT_AU(fastq_ch, ref_mt_au_path)
        SORT_MT_AU(ALIGN_MT_AU.out.sam)

        def tag_mt_au_input = SORT_MT_AU.out.bam_bai
            .join(demux_bams)
            .map { sample_id, sorted_bam, sorted_bai, raw_bam ->
                tuple(sample_id, sorted_bam, sorted_bai, raw_bam)
            }
        TAG_MT_AU(tag_mt_au_input)
        ch_tagged_mt_au_bam_bai = TAG_MT_AU.out.tagged_bam_bai

        CALL_MT_AU(ch_tagged_mt_au_bam_bai, ref_mt_au_path)
        ch_vcf_mt_au = CALL_MT_AU.out.vcf_tbi

        MOD_MT_AU(ch_tagged_mt_au_bam_bai, ref_mt_au_path)
        ch_bed_mt_au = MOD_MT_AU.out.bed
    }

    emit:
    tagged_nuc_bam_bai   = ch_tagged_nuc_bam_bai
    tagged_mt_ap_bam_bai = ch_tagged_mt_ap_bam_bai
    tagged_mt_au_bam_bai = ch_tagged_mt_au_bam_bai
    vcf_nuc              = ch_vcf_nuc
    vcf_mt_ap            = ch_vcf_mt_ap
    vcf_mt_au            = ch_vcf_mt_au
    bed_nuc              = ch_bed_nuc
    bed_mt_ap            = ch_bed_mt_ap
    bed_mt_au            = ch_bed_mt_au
    phased_bam_bai       = ch_phased_bam_bai
    versions             = ch_versions
}
