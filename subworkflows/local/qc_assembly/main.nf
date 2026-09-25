// subworkflows/local/qc_assembly/main.nf
include { CRAMINO } from '../../../modules/nf-core/cramino'
include { SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq'
include { FASTQC } from '../../../modules/nf-core/fastqc'
include { NANOPLOT } from '../../../modules/nf-core/nanoplot'
include { FILTLONG } from '../../../modules/nf-core/filtlong'
include { FLYE } from '../../../modules/nf-core/flye'
include { MEDAKA } from '../../../modules/nf-core/medaka'
include { QUAST } from '../../../modules/nf-core/quast'
include { BUSCO } from '../../../modules/nf-core/busco'

workflow QC_ASSEMBLY {
    take:
    demux_bams // tuple(sample_id, bam)

    main:
    ch_versions = Channel.empty()

    CRAMINO(demux_bams)
    ch_versions = ch_versions.mix(CRAMINO.out.versions)

    SAMTOOLS_FASTQ(demux_bams)
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)
    def fastq_ch = SAMTOOLS_FASTQ.out.fastq

    FASTQC(fastq_ch)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    NANOPLOT(fastq_ch)
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)

    def polished_assembly = Channel.empty()
    def quast_ch = Channel.empty()
    def busco_ch = Channel.empty()

    if (params.run_type == 'de_novo_assembly') {
        FILTLONG(fastq_ch)
        ch_versions = ch_versions.mix(FILTLONG.out.versions)

        FLYE(FILTLONG.out.filtered_reads)
        ch_versions = ch_versions.mix(FLYE.out.versions)
        
        def medaka_input = FILTLONG.out.filtered_reads.join(FLYE.out.assembly)
        MEDAKA(medaka_input)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
        polished_assembly = MEDAKA.out.polished_assembly

        QUAST(polished_assembly)
        ch_versions = ch_versions.mix(QUAST.out.versions)

        BUSCO(polished_assembly)
        ch_versions = ch_versions.mix(BUSCO.out.versions)
        
        quast_ch = QUAST.out.report
        busco_ch = BUSCO.out.report
    }

    emit:
    fastq          = fastq_ch
    assembly       = polished_assembly
    fastqc_files   = FASTQC.out.qc_files
    nanoplot_files = NANOPLOT.out.report
    cramino_files  = CRAMINO.out.report
    quast_files    = quast_ch
    busco_files    = busco_ch
    versions       = ch_versions
}
