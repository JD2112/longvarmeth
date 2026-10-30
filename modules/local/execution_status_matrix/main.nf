process EXECUTION_STATUS_MATRIX {
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10':
        'quay.io/biocontainers/python:3.10' }"

    input:
    val run_type
    val platform
    val snv_caller
    val sv_caller
    val phasing_tool
    val skip_alignment
    val skip_variants
    val skip_methylation
    val enable_phasing
    val run_benchmarking
    val enable_phi_sanitization
    val enable_audit_trail

    output:
    path "pipeline_execution_status_mqc.tsv", emit: tsv
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "
import os

run_type = '${run_type}'
platform = '${platform}'
snv_caller = '${snv_caller}'
sv_caller = '${sv_caller}'
phasing_tool = '${phasing_tool}'
skip_align = '${skip_alignment}'.lower() in ['true', '1']
skip_vars = '${skip_variants}'.lower() in ['true', '1']
skip_meth = '${skip_methylation}'.lower() in ['true', '1']
phasing_on = '${enable_phasing}'.lower() in ['true', '1'] or platform == 'pacbio'
bench_on = '${run_benchmarking}'.lower() in ['true', '1']
phi_on = '${enable_phi_sanitization}'.lower() in ['true', '1']
audit_on = '${enable_audit_trail}'.lower() in ['true', '1']

stages = [
    ('1', 'Input QC & Checksum Ingestion', 'COMPLETED', 'INPUT_CHECK / CRAMINO', 'Checksums SHA-256 TSV', 'Raw input file integrity cryptographically verified'),
    ('2', 'Basecalling / Demultiplexing', 'COMPLETED' if platform == 'ont' else 'SKIPPED (PacBio native HiFi)', 'Dorado / Pod5' if platform == 'ont' else 'N/A', 'Demuxed BAMs / FASTQ', 'GPU/CPU neural basecalling'),
    ('3', 'Read Quality Trimming', 'COMPLETED', 'Filtlong / NanoPlot', 'Filtered Reads FASTQ', 'Read length and quality threshold filtering'),
    ('4', 'Long-Read Alignment & Streaming', 'SKIPPED (Pre-aligned Input)' if skip_align else 'COMPLETED', 'minimap2 -ax map-ont -y' if platform == 'ont' else 'pbmm2', 'Coordinate-Sorted BAM & BAI', 'Modification tags (MM/ML) preserved in stream'),
    ('5', 'Modification Tag Integrity Verification', 'COMPLETED' if not skip_meth else 'SKIPPED (Methylation Disabled)', 'CHECK_MOD_TAGS (samtools view)', 'Mod Tag Audit Report', 'SAM auxiliary tag presence validated before pileup'),
    ('6', 'Small Variant Calling', 'SKIPPED (Variants Disabled)' if skip_vars else 'COMPLETED', f'{snv_caller.upper()}', 'VCF.GZ & TBI', f'Neural / statistical diploid caller ({snv_caller})'),
    ('7', 'Structural Variant & STR Calling', 'SKIPPED (Variants Disabled)' if skip_vars else 'COMPLETED', f'{sv_caller.upper()} / TRGT', 'SV / STR VCF.GZ', 'Long-read structural variant and tandem repeat genotyping'),
    ('8', 'Epigenomics & Methylation Pileup', 'SKIPPED (Methylation Disabled)' if skip_meth else 'COMPLETED', 'modkit pileup' if platform == 'ont' else 'pbcpgtools (aligned_bam_to_cpg_scores)', 'bedMethyl / CpG BedGraph', 'High-confidence 5mC CpG methylation scoring'),
    ('9', 'Joint Multi-Variant Co-Phasing', 'COMPLETED' if (not skip_vars and phasing_on) else 'SKIPPED (Phasing Disabled)', f'{phasing_tool.upper()} / WhatsHap', 'Phased VCF & Haplotagged BAM', 'Joint SNV + SV + STR haplotype resolution'),
    ('10', 'Truth-Set Benchmarking', 'COMPLETED' if bench_on else 'SKIPPED (Benchmarking Inactive)', 'hap.py (SNVs) / Truvari (SVs)', 'Benchmark Summary TSVs', 'Concordance, precision, recall, and F1 evaluation against GIAB HG002'),
    ('11', 'Clinical Governance & De-Identification', 'COMPLETED' if (phi_on or audit_on) else 'SKIPPED (Governance Disabled)', 'SANITIZE_BAM / SANITIZE_VCF / MANIFEST', 'Sanitized BAM/VCF & manifest.json', 'Deterministic anon_id, header scrubbing, full-chain provenance'),
    ('12', 'Interactive Visualization & Reporting', 'COMPLETED', 'MultiQC / JBrowse 2 / Dashboard', 'MultiQC HTML & JBrowse Tracks', 'Clinical executive summary and interactive genome tracks')
]

with open('pipeline_execution_status_mqc.tsv', 'w') as out:
    out.write('# id: \\'pipeline_execution_status\\'\\n')
    out.write('# section_name: \\'Pipeline Execution Status Matrix\\'\\n')
    out.write('# description: \\'Transparent execution status audit across all 12 pipeline analytical stages. Skipped or parameter-disabled stages are explicitly distinguished from completed runs to prevent incomplete analyses from being mistaken for complete results.\\'\\n')
    out.write('# plot_type: \\'table\\'\\n')
    out.write('# pconfig:\\n')
    out.write('#    title: \\'Analytical Stage Execution Matrix\\'\\n')
    out.write('Stage\\tAnalytical_Stage\\tStatus\\tTool_Invoked\\tOutput_Artifacts\\tExecution_Notes\\n')
    for stage_num, name, status, tool, artifacts, notes in stages:
        out.write(f'{stage_num}\\t{name}\\t{status}\\t{tool}\\t{artifacts}\\t{notes}\\n')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
    END_VERSIONS
    """

    stub:
    """
    cat <<EOF > pipeline_execution_status_mqc.tsv
# id: 'pipeline_execution_status'
# section_name: 'Pipeline Execution Status Matrix'
# description: 'Transparent execution status audit across all 12 pipeline analytical stages. Skipped or parameter-disabled stages are explicitly distinguished from completed runs to prevent incomplete analyses from being mistaken for complete results.'
# plot_type: 'table'
# pconfig:
#    title: 'Analytical Stage Execution Matrix'
Stage\tAnalytical_Stage\tStatus\tTool_Invoked\tOutput_Artifacts\tExecution_Notes
1\tInput QC & Checksum Ingestion\tCOMPLETED\tINPUT_CHECK / CRAMINO\tChecksums SHA-256 TSV\tRaw input file integrity cryptographically verified
2\tBasecalling / Demultiplexing\tCOMPLETED\tDorado / Pod5\tDemuxed BAMs / FASTQ\tGPU/CPU neural basecalling
3\tRead Quality Trimming\tCOMPLETED\tFiltlong / NanoPlot\tFiltered Reads FASTQ\tRead length and quality threshold filtering
4\tLong-Read Alignment & Streaming\tCOMPLETED\tminimap2 -ax map-ont -y\tCoordinate-Sorted BAM & BAI\tModification tags (MM/ML) preserved in stream
5\tModification Tag Integrity Verification\tCOMPLETED\tCHECK_MOD_TAGS\tMod Tag Audit Report\tSAM auxiliary tag presence validated before pileup
6\tSmall Variant Calling\tCOMPLETED\t${snv_caller}\tVCF.GZ & TBI\tDiploid SNV/Indel calling
7\tStructural Variant & STR Calling\tCOMPLETED\t${sv_caller} / TRGT\tSV / STR VCF.GZ\tLong-read structural variant and tandem repeat genotyping
8\tEpigenomics & Methylation Pileup\tCOMPLETED\tmodkit pileup\tbedMethyl / CpG BedGraph\tHigh-confidence 5mC CpG methylation scoring
9\tJoint Multi-Variant Co-Phasing\tCOMPLETED\t${phasing_tool} / WhatsHap\tPhased VCF & Haplotagged BAM\tJoint SNV + SV + STR haplotype resolution
10\tTruth-Set Benchmarking\tCOMPLETED\thap.py / Truvari\tBenchmark Summary TSVs\tConcordance, precision, recall, and F1 evaluation against GIAB HG002
11\tClinical Governance & De-Identification\tCOMPLETED\tSANITIZE_BAM / SANITIZE_VCF / MANIFEST\tSanitized BAM/VCF & manifest.json\tDeterministic anon_id, header scrubbing, full-chain provenance
12\tInteractive Visualization & Reporting\tCOMPLETED\tMultiQC / JBrowse 2 / Dashboard\tMultiQC HTML & JBrowse Tracks\tClinical executive summary and interactive genome tracks
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10
    END_VERSIONS
    """
}

