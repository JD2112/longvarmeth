// subworkflows/local/input_check/main.nf
// Step 1: Input file verification and cryptographic SHA-256 ingestion checksums

process CHECK_INPUT_CHECKSUM {
    tag "$sample_id"
    label 'process_single'

    input:
    tuple val(sample_id), path(input_file), val(input_type)

    output:
    tuple val(sample_id), path("${sample_id}_input.sha256"), emit: checksum
    path "versions.yml"                                    , emit: versions

    script:
    """
    if [ -d "${input_file}" ]; then
        find -L ${input_file} -type f -exec sha256sum {} + > ${sample_id}_input.sha256 || true
    else
        sha256sum ${input_file} > ${sample_id}_input.sha256
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$( sha256sum --version | head -n1 | sed 's/sha256sum (GNU coreutils) //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_input.sha256
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: 9.1
    END_VERSIONS
    """
}

workflow INPUT_CHECK {
    take:
    samples_ch // tuple(sample_id, input_path, input_type, kit, model, modified_bases_model)

    main:
    ch_versions = Channel.empty()

    def check_inputs = samples_ch.map { sample_id, input_path, input_type, kit, model, mod_model ->
        tuple(sample_id, file(input_path), input_type)
    }

    CHECK_INPUT_CHECKSUM(check_inputs)
    ch_versions = ch_versions.mix(CHECK_INPUT_CHECKSUM.out.versions)

    emit:
    checksums = CHECK_INPUT_CHECKSUM.out.checksum
    versions  = ch_versions
}
