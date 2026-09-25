// subworkflows/local/basecall_demux/main.nf
include { POD5_CONVERT } from '../../../modules/nf-core/pod5/convert'
include { DORADO_BASECALL } from '../../../modules/local/dorado/basecall'
include { DORADO_DEMUX } from '../../../modules/local/dorado/demux'

workflow BASECALL_DEMUX {
    take:
    samples_ch // tuple(run_id, input_path, input_type, kit, model, modified_bases_model)

    main:
    ch_versions = Channel.empty()

    def pod5_runs = samples_ch.filter { it[2] == 'pod5' }
    def fast5_runs = samples_ch.filter { it[2] == 'fast5' }

    POD5_CONVERT(fast5_runs.map { tuple(it[0], it[1]) })
    ch_versions = ch_versions.mix(POD5_CONVERT.out.versions)

    def converted_runs = POD5_CONVERT.out.pod5_dir
        .join(fast5_runs.map { tuple(it[0], it[2], it[3], it[4], it[5]) })
        .map { tuple(it[0], it[1], it[3], it[4], it[5]) }

    def all_pod5_runs = pod5_runs.map { tuple(it[0], it[1], it[3], it[4], it[5]) }
        .concat(converted_runs)

    DORADO_BASECALL(all_pod5_runs)
    
    def demux_inputs = DORADO_BASECALL.out.bam
        .join(samples_ch.map { tuple(it[0], it[3]) })
    
    DORADO_DEMUX(demux_inputs)

    def demux_bams = DORADO_DEMUX.out.bams
        .transpose()
        .map { run_id, bam ->
            def barcode = bam.name.replaceAll(/.*_(barcode\d+)\.bam/, '$1')
            if (barcode == bam.name) {
                barcode = "unclassified"
            }
            return tuple("${run_id}_${barcode}", bam)
        }
        .filter { it[0] !=~ /.*unclassified.*/ }

    emit:
    demux_bams // tuple(sample_id, bam)
    versions   = ch_versions
}
