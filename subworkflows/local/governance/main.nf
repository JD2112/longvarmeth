// subworkflows/local/governance/main.nf
// Strict clinical governance, PHI sanitization, and deterministic audit trail
// Inspired by JD2112/milou architecture

process SANITIZE_METADATA {
    tag "$sample_id"
    label 'process_single'

    conda "bioconda::samtools=1.19.2 bioconda::bcftools=1.19"
    container "quay.io/biocontainers/mulled-v2-42f3ea21636d68b3a3822d52504b022089b21618:0f8ab8c73856d10f2142e05f63d662c5fb3ea8c7-0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    val enable_sanitization

    output:
    tuple val(sample_id), path("${sample_id}.sanitized.bam"), path("${sample_id}.sanitized.bam.bai"), emit: sanitized_bam_bai
    path "versions.yml"                                                                             , emit: versions

    script:
    if (enable_sanitization) {
        """
        # Strip PHI/patient-identifying comments and unneeded headers from BAM
        samtools view -H ${bam} | grep -v -E "(@CO|SM:|DS:.*patient|DS:.*hospital)" > clean_header.sam || true
        samtools reheader clean_header.sam ${bam} > ${sample_id}.sanitized.bam
        samtools index ${sample_id}.sanitized.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$( samtools --version | head -n1 | sed 's/samtools //' )
        END_VERSIONS
        """
    } else {
        """
        ln -s ${bam} ${sample_id}.sanitized.bam
        ln -s ${bai} ${sample_id}.sanitized.bam.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$( samtools --version | head -n1 | sed 's/samtools //' )
        END_VERSIONS
        """
    }

    stub:
    """
    touch ${sample_id}.sanitized.bam
    touch ${sample_id}.sanitized.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.19.2
    END_VERSIONS
    """
}

process GENERATE_AUDIT_MANIFEST {
    tag "$sample_id"
    label 'process_single'

    input:
    tuple val(sample_id), path(artifacts)
    val workflow_params

    output:
    tuple val(sample_id), path("${sample_id}_manifest.json"), emit: manifest
    path "${sample_id}_checksums.sha256"                    , emit: checksums

    script:
    def artifacts_list = artifacts instanceof List ? artifacts.join(' ') : artifacts.toString()
    """
    # 1. Compute SHA-256 checksums for each output artifact
    sha256sum ${artifacts_list} > ${sample_id}_checksums.sha256

    # 2. Build structured JSON audit record
    cat <<-EOF > ${sample_id}_manifest.json
    {
        "pipeline": "longvarmeth",
        "sample_id": "${sample_id}",
        "timestamp_utc": "\$(date -u +'%Y-%m-%dT%H:%M:%SZ')",
        "workflow_commit": "${workflow.commit ?: 'dev'}",
        "workflow_version": "${workflow.manifest.version ?: '1.0.0'}",
        "parameters": {
            "platform": "${params.platform}",
            "species": "${params.species}",
            "run_type": "${params.run_type}",
            "snv_caller": "${params.snv_caller}",
            "sv_caller": "${params.sv_caller}",
            "phasing_tool": "${params.phasing_tool}"
        },
        "artifacts": [
            \$(awk '{ printf "{\\"file\\": \\"%s\\", \\"sha256\\": \\"%s\\"},\\n", \$2, \$1 }' ${sample_id}_checksums.sha256 | sed '\$ s/,\$//')
        ]
    }
    EOF
    """

    stub:
    """
    touch ${sample_id}_checksums.sha256
    echo '{"pipeline":"longvarmeth","sample_id":"${sample_id}"}' > ${sample_id}_manifest.json
    """
}

workflow CLINICAL_GOVERNANCE {
    take:
    ch_bams           // tuple(sample_id, bam, bai)
    ch_final_outputs  // tuple(sample_id, [files])

    main:
    ch_versions = Channel.empty()

    SANITIZE_METADATA(ch_bams, params.enable_phi_sanitization ?: false)
    ch_versions = ch_versions.mix(SANITIZE_METADATA.out.versions)

    GENERATE_AUDIT_MANIFEST(ch_final_outputs, params)

    emit:
    sanitized_bams = SANITIZE_METADATA.out.sanitized_bam_bai
    manifests      = GENERATE_AUDIT_MANIFEST.out.manifest
    checksums      = GENERATE_AUDIT_MANIFEST.out.checksums
    versions       = ch_versions
}
