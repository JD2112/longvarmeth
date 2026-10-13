// subworkflows/local/visualization/main.nf
// Step 11: Interactive visualization (JBrowse 2 and Methylartist locus plots)

include { METHYLARTIST_LOCUS } from '../../../modules/local/methylartist'
include { JBROWSE_CREATE } from '../../../modules/local/jbrowse'
include { METHYLATION_DASHBOARD } from '../../../modules/local/methylation_dashboard'

workflow VISUALIZATION {
    take:
    assembly_ref     // tuple(sample_id, [ref_fasta, ref_fai])
    tagged_bams      // tuple(sample_id, bam, bai)
    vcf_tbi          // tuple(sample_id, vcf, tbi)
    methylation_data // tuple(sample_id, bed_or_tsv) (optional, can be empty channel)

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // 1. Methylartist Locus Plotting
    if (params.methylartist_locus) {
        def methylartist_input = tagged_bams
            .join(assembly_ref)
            .map { sample_id, bam, bai, ref_files ->
                tuple(
                    tuple(sample_id, bam, bai),
                    ref_files[0],
                    params.methylartist_locus
                )
            }
        METHYLARTIST_LOCUS(
            methylartist_input.map { bam_tuple, ref, locus -> bam_tuple },
            methylartist_input.map { bam_tuple, ref, locus -> ref },
            methylartist_input.map { bam_tuple, ref, locus -> locus }
        )
        ch_versions = ch_versions.mix(METHYLARTIST_LOCUS.out.versions)
    }

    // 2. Automated JBrowse 2 Directory Generation
    def jbrowse_input = tagged_bams
        .join(vcf_tbi)
        .join(assembly_ref)
        .map { sample_id, bam, bai, vcf, tbi, ref_files ->
            tuple(sample_id, ref_files, bam, bai, vcf, tbi)
        }

    JBROWSE_CREATE(jbrowse_input)
    ch_versions = ch_versions.mix(JBROWSE_CREATE.out.versions)


    // 3. Clinical Methylation Dashboard (adapted from gms-nallo)
    if (params.enable_methylation_dashboard) {
        def template_file = file("${projectDir}/assets/methylation_dashboard/dashboard-template.html")
        METHYLATION_DASHBOARD(methylation_data, template_file)
        ch_versions = ch_versions.mix(METHYLATION_DASHBOARD.out.versions)
        ch_reports  = ch_reports.mix(METHYLATION_DASHBOARD.out.html)
    }

    emit:
    jbrowse_dir       = JBROWSE_CREATE.out.jbrowse_dir
    methylation_html  = ch_reports
    versions          = ch_versions
}
