// subworkflows/local/methylation/main.nf
// Step 7: Epigenomics and base modifications (modkit for ONT; pb-cpg-tools for PacBio HiFi)

include { MODKIT_PILEUP as MOD_NUC } from '../../../modules/nf-core/modkit/pileup'
include { PBCPGTOOLS } from '../../../modules/nf-core/pbcpgtools'

process CHECK_MOD_TAGS {
    tag "$sample_id"
    label 'process_single'

    conda "bioconda::samtools=1.23.1"
    container "quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.mod_tags_qc.txt"), emit: qc_report
    path "versions.yml"                                       , emit: versions

    script:
    """
    echo "=== Kinetic Modification Tag Check for ${sample_id} ===" > ${sample_id}.mod_tags_qc.txt
    MM_COUNT=\$(samtools view ${bam} | head -n 1000 | grep -c "MM:Z:" || true)
    ML_COUNT=\$(samtools view ${bam} | head -n 1000 | grep -c "ML:B:" || true)
    echo "Sampled first 1000 alignments:" >> ${sample_id}.mod_tags_qc.txt
    echo "  Alignments with MM tag: \${MM_COUNT} / 1000" >> ${sample_id}.mod_tags_qc.txt
    echo "  Alignments with ML tag: \${ML_COUNT} / 1000" >> ${sample_id}.mod_tags_qc.txt
    if [ "\${MM_COUNT}" -eq 0 ]; then
        echo "[WARN] No MM base-modification tags detected in first 1000 reads of ${bam}!" >> ${sample_id}.mod_tags_qc.txt
    else
        echo "[INFO] Modification tags validated successfully." >> ${sample_id}.mod_tags_qc.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$( samtools --version 2>&1 | head -n1 | sed -e 's/samtools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.mod_tags_qc.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """
}

workflow METHYLATION {
    take:
    aligned_bam_bai // tuple(sample_id, bam, bai)
    ref_nuc_path    // path(reference_fasta)

    main:
    ch_versions = Channel.empty()
    def ch_bed_methylation = Channel.empty()
    def ch_methylation_qc  = Channel.empty()

    // 1. Tag validation QC
    CHECK_MOD_TAGS(aligned_bam_bai)
    ch_versions       = ch_versions.mix(CHECK_MOD_TAGS.out.versions)
    ch_methylation_qc = CHECK_MOD_TAGS.out.qc_report

    // 2. Platform-specialized Base Modification Analysis
    if (params.platform == 'pacbio') {
        // Native PacBio HiFi CpG aggregation with pb-cpg-tools
        PBCPGTOOLS(aligned_bam_bai, ref_nuc_path)
        ch_versions        = ch_versions.mix(PBCPGTOOLS.out.versions)
        ch_bed_methylation = PBCPGTOOLS.out.bed
    } else {
        // Oxford Nanopore modkit whole-genome 5mC / 5hmC pileup
        MOD_NUC(aligned_bam_bai, ref_nuc_path)
        ch_versions        = ch_versions.mix(MOD_NUC.out.versions)
        ch_bed_methylation = MOD_NUC.out.bed
    }

    emit:
    bed_methylation = ch_bed_methylation
    qc_report       = ch_methylation_qc
    versions        = ch_versions
}
