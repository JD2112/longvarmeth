// subworkflows/local/assembly/main.nf
// Step 3 (de novo): De novo assembly, consensus polishing, and completeness QC

include { FILTLONG } from '../../../modules/nf-core/filtlong'
include { FLYE } from '../../../modules/nf-core/flye'
include { MEDAKA } from '../../../modules/nf-core/medaka'
include { QUAST } from '../../../modules/nf-core/quast'
include { BUSCO } from '../../../modules/nf-core/busco'

workflow DE_NOVO_ASSEMBLY {
    take:
    fastq_ch // tuple(sample_id, fastq)

    main:
    ch_versions = Channel.empty()

    // 1. Quality & Length filtering
    FILTLONG(fastq_ch)
    ch_versions = ch_versions.mix(FILTLONG.out.versions)

    // 2. Long-read de novo assembly
    FLYE(FILTLONG.out.filtered_reads)
    ch_versions = ch_versions.mix(FLYE.out.versions)
    
    // 3. Consensus Polishing (Medaka)
    def medaka_input = FILTLONG.out.filtered_reads.join(FLYE.out.assembly)
    MEDAKA(medaka_input)
    ch_versions = ch_versions.mix(MEDAKA.out.versions)
    def polished_assembly = MEDAKA.out.polished_assembly

    // 4. Assembly Metrics (QUAST)
    QUAST(polished_assembly)
    ch_versions = ch_versions.mix(QUAST.out.versions)

    // 5. Gene Completeness Benchmark (BUSCO)
    BUSCO(polished_assembly)
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    emit:
    assembly      = polished_assembly
    quast_reports = QUAST.out.report
    busco_reports = BUSCO.out.report
    versions      = ch_versions
}
