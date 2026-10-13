// subworkflows/local/qc/main.nf
// Step 2: Pre-alignment read quality assessment

include { CRAMINO } from '../../../modules/nf-core/cramino'
include { SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq'
include { FASTQC } from '../../../modules/nf-core/fastqc'
include { NANOPLOT } from '../../../modules/nf-core/nanoplot'

workflow READ_QC {
    take:
    demux_bams // tuple(sample_id, bam)

    main:
    ch_versions = Channel.empty()

    // 1. Yield & BAM Quality metrics
    CRAMINO(demux_bams)
    ch_versions = ch_versions.mix(CRAMINO.out.versions)

    // 2. Extract FASTQ for QC tools that require sequence streams
    SAMTOOLS_FASTQ(demux_bams)
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)
    def fastq_ch = SAMTOOLS_FASTQ.out.fastq

    // 3. FastQC
    FASTQC(fastq_ch)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // 4. NanoPlot Long-Read Length & Q-Score Profiling
    NANOPLOT(fastq_ch)
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)

    emit:
    fastq          = fastq_ch
    cramino_files  = CRAMINO.out.tsv
    fastqc_files   = FASTQC.out.zip
    nanoplot_files = NANOPLOT.out.txt
    versions       = ch_versions
}
