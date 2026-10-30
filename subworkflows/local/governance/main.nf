// subworkflows/local/governance/main.nf
// Strict clinical governance, end-to-end PHI sanitization, and deterministic forensic audit trail
// Inspired by JD2112/milou architecture

process SANITIZE_BAM {
    tag "$anon_id"
    label 'process_single'

    conda "bioconda::samtools=1.23.1"
    container "quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"

    input:
    tuple val(sample_id), val(anon_id), path(bam), path(bai)
    val enable_sanitization

    output:
    tuple val(anon_id), path("${anon_id}.sanitized.bam"), path("${anon_id}.sanitized.bam.bai"), emit: sanitized_bam_bai
    path "versions.yml"                                                                      , emit: versions

    script:
    if (enable_sanitization) {
        """
        # 1. Extract header and strip PHI, local paths, and hospital identifiers
        samtools view -H ${bam} \\
            | grep -v -E "(@CO|DS:.*patient|DS:.*hospital|DS:.*barcode)" \\
            | sed "s/SM:${sample_id}/SM:${anon_id}/g" > clean_header.sam

        # 2. Reheader and index BAM
        samtools reheader clean_header.sam ${bam} > ${anon_id}.sanitized.bam
        samtools index -@ ${task.cpus ?: 2} ${anon_id}.sanitized.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$( samtools --version 2>&1 | head -n1 | sed -e 's/samtools //' )
        END_VERSIONS
        """
    } else {
        """
        ln -s ${bam} ${anon_id}.sanitized.bam
        ln -s ${bai} ${anon_id}.sanitized.bam.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$( samtools --version 2>&1 | head -n1 | sed -e 's/samtools //' )
        END_VERSIONS
        """
    }

    stub:
    """
    touch ${anon_id}.sanitized.bam
    touch ${anon_id}.sanitized.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """
}

process SANITIZE_VCF {
    tag "$anon_id"
    label 'process_single'

    conda "bioconda::bcftools=1.23.1"
    container "quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0"

    input:
    tuple val(sample_id), val(anon_id), path(vcf), path(tbi)
    val enable_sanitization

    output:
    tuple val(anon_id), path("${anon_id}.sanitized.vcf.gz"), path("${anon_id}.sanitized.vcf.gz.tbi"), emit: sanitized_vcf_tbi
    path "versions.yml"                                                                             , emit: versions

    script:
    if (enable_sanitization) {
        """
        # 1. Create sample replacement map
        echo "${sample_id} ${anon_id}" > sample_name_map.txt

        # 2. Extract and sanitize header comments (strip command lines revealing host paths)
        bcftools view -h ${vcf} | grep -v -E "(##CL=|##commandline=|##bcftools_callCommand=|##contig.*URL=)" > clean_vcf_header.txt

        # 3. Reheader with clean headers and anonymized sample identifier
        bcftools reheader -h clean_vcf_header.txt -s sample_name_map.txt ${vcf} -o ${anon_id}.sanitized.vcf.gz
        bcftools index -t ${anon_id}.sanitized.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$( bcftools --version 2>&1 | head -n1 | sed -e 's/bcftools //' )
        END_VERSIONS
        """
    } else {
        """
        ln -s ${vcf} ${anon_id}.sanitized.vcf.gz
        ln -s ${tbi} ${anon_id}.sanitized.vcf.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$( bcftools --version 2>&1 | head -n1 | sed -e 's/bcftools //' )
        END_VERSIONS
        """
    }

    stub:
    """
    touch ${anon_id}.sanitized.vcf.gz
    touch ${anon_id}.sanitized.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.23.1
    END_VERSIONS
    """
}

process GENERATE_AUDIT_MANIFEST {
    tag "$anon_id"
    label 'process_single'

    input:
    tuple val(sample_id), val(anon_id), path(artifacts)
    path ref_nuc
    path input_checksums
    val workflow_params

    output:
    tuple val(anon_id), path("${anon_id}_manifest.json"), emit: manifest
    path "${anon_id}_checksums.sha256"                  , emit: checksums
    path "sample_mapping.audit.tsv"                     , optional: true, emit: mapping_audit

    script:
    def artifacts_list = artifacts instanceof List ? artifacts.join(' ') : artifacts.toString()
    def input_chk_list = input_checksums instanceof List ? input_checksums.join(' ') : input_checksums.toString()
    """
    # 1. Compute reference FASTA checksum
    REF_SHA=\$(sha256sum ${ref_nuc} | awk '{print \$1}')

    # 2. Extract input checksum recorded in Step 1
    INPUT_SHA=\$(cat ${input_chk_list} | awk '{print \$1}' | head -n1)

    # 3. Compute SHA-256 checksums for each final output artifact
    sha256sum ${artifacts_list} > ${anon_id}_checksums.sha256

    # 4. Record restricted de-identification audit mapping
    echo -e "${sample_id}\\t${anon_id}\\t\$(date -u +'%Y-%m-%dT%H:%M:%SZ')" > sample_mapping.audit.tsv

    # 5. Build comprehensive, publication-ready forensic JSON manifest
    cat <<-EOF > ${anon_id}_manifest.json
    {
        "pipeline": "longvarmeth",
        "sample_id": "${sample_id}",
        "anonymized_id": "${anon_id}",
        "timestamp_utc": "\$(date -u +'%Y-%m-%dT%H:%M:%SZ')",
        "workflow": {
            "version": "${workflow.manifest.version ?: '1.0.0'}",
            "commit": "${workflow.commitId ?: 'dev'}",
            "container_engine": "${workflow.containerEngine ?: 'none'}",
            "profile": "${workflow.profile ?: 'standard'}"
        },
        "provenance": {
            "input_checksum_sha256": "\${INPUT_SHA}",
            "reference_fasta": "${ref_nuc.name}",
            "reference_sha256": "\${REF_SHA}"
        },
        "parameters": {
            "platform": "${params.platform}",
            "species": "${params.species}",
            "run_type": "${params.run_type}",
            "snv_caller": "${params.snv_caller}",
            "sv_caller": "${params.sv_caller}",
            "phasing_tool": "${params.phasing_tool}",
            "phi_sanitization_enabled": ${params.enable_phi_sanitization ? true : false}
        },
        "artifacts": [
            \$(awk '{ printf "{\\"file\\": \\"%s\\", \\"sha256\\": \\"%s\\"},\\n", \$2, \$1 }' ${anon_id}_checksums.sha256 | sed '\$ s/,\$//')
        ]
    }
    EOF
    """

    stub:
    """
    touch ${anon_id}_checksums.sha256
    touch sample_mapping.audit.tsv
    echo '{"pipeline":"longvarmeth","anonymized_id":"${anon_id}"}' > ${anon_id}_manifest.json
    """
}

workflow CLINICAL_GOVERNANCE {
    take:
    ch_bams            // tuple(sample_id, bam, bai)
    ch_vcfs            // tuple(sample_id, vcf, tbi)
    ch_input_checksums // tuple(sample_id, input_sha256_file)
    ref_nuc_path       // path(reference_fasta)
    ch_extra_artifacts // tuple(sample_id, [files])

    main:
    ch_versions = Channel.empty()

    // 1. Generate deterministic anonymized identifier (anon_id)
    def ch_sample_tuples = ch_bams
        .map { sample_id, bam, bai ->
            def anon_id = params.enable_phi_sanitization ?
                "ANON_" + java.security.MessageDigest.getInstance("SHA-256").digest(sample_id.bytes).encodeHex().toString().substring(0, 10) :
                sample_id
            tuple(sample_id, anon_id, bam, bai)
        }

    // 2. Sanitize BAM headers
    SANITIZE_BAM(ch_sample_tuples, params.enable_phi_sanitization ?: false)
    ch_versions = ch_versions.mix(SANITIZE_BAM.out.versions)

    // 3. Sanitize VCF headers and replace sample names
    def ch_vcf_tuples = ch_vcfs
        .map { sample_id, vcf, tbi ->
            def anon_id = params.enable_phi_sanitization ?
                "ANON_" + java.security.MessageDigest.getInstance("SHA-256").digest(sample_id.bytes).encodeHex().toString().substring(0, 10) :
                sample_id
            tuple(sample_id, anon_id, vcf, tbi)
        }
    SANITIZE_VCF(ch_vcf_tuples, params.enable_phi_sanitization ?: false)
    ch_versions = ch_versions.mix(SANITIZE_VCF.out.versions)

    // 4. Build comprehensive forensic audit manifest
    def ch_all_input_checksums = ch_input_checksums.map { id, file -> file }.collect()
    
    // Join sanitized BAM and VCF outputs by sample index
    def ch_sanitized_pairs = ch_sample_tuples
        .map { sample_id, anon_id, bam, bai -> tuple(anon_id, sample_id) }
        .join(SANITIZE_BAM.out.sanitized_bam_bai)
        .join(SANITIZE_VCF.out.sanitized_vcf_tbi)
        .map { anon_id, sample_id, s_bam, s_bai, s_vcf, s_tbi ->
            tuple(sample_id, anon_id, [s_bam, s_bai, s_vcf, s_tbi])
        }

    def ch_manifest_input = ch_sanitized_pairs
        .join(ch_extra_artifacts, remainder: true)
        .map { row ->
            def sample_id     = row[0]
            def anon_id       = row[1]
            def artifacts     = row[2] instanceof List ? row[2] : [row[2]]
            def extra         = (row.size() > 3 && row[3] != null) ? (row[3] instanceof List ? row[3] : [row[3]]) : []
            def all_artifacts = (artifacts + extra).flatten().findAll { it != null }
            tuple(sample_id, anon_id, all_artifacts)
        }

    GENERATE_AUDIT_MANIFEST(ch_manifest_input, ref_nuc_path, ch_all_input_checksums, params)

    emit:
    sanitized_bams    = SANITIZE_BAM.out.sanitized_bam_bai
    sanitized_vcfs    = SANITIZE_VCF.out.sanitized_vcf_tbi
    manifests         = GENERATE_AUDIT_MANIFEST.out.manifest
    checksums         = GENERATE_AUDIT_MANIFEST.out.checksums
    mapping_audit     = GENERATE_AUDIT_MANIFEST.out.mapping_audit
    versions          = ch_versions
}
