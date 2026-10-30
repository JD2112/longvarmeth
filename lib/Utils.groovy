// lib/Utils.groovy
// Helper utilities for parameter validation, samplesheet auditing, and run completion

class Utils {

    private static final List<String> VALID_RUN_TYPES = ['reference_based', 'de_novo_assembly']
    private static final List<String> VALID_PLATFORMS = ['ont', 'pacbio', 'auto']
    private static final List<String> VALID_SPECIES   = ['human', 'eukaryote', 'prokaryote', 'virus']

    private static final List<String> ONT_VALID_INPUT_TYPES = [
        'pod5', 'fast5', 'bam', 'ubam', 'unaligned_bam', 'aligned_bam', 'fastq', 'fastq.gz', 'fq', 'fq.gz'
    ]
    private static final List<String> PACBIO_VALID_INPUT_TYPES = [
        'pacbio_bam', 'bam', 'ubam', 'hifi', 'pacbio_hifi', 'unaligned_bam', 'aligned_bam', 'fastq', 'fastq.gz', 'fq', 'fq.gz'
    ]

    public static void validateParameters(params, log) {
        def errors = []

        // 1. Validate platform
        if (params.platform && !VALID_PLATFORMS.contains(params.platform.toLowerCase())) {
            errors << "Invalid --platform '${params.platform}'. Supported: ${VALID_PLATFORMS.join(', ')}"
        }

        // 2. Validate run_type
        if (params.run_type && !VALID_RUN_TYPES.contains(params.run_type.toLowerCase())) {
            errors << "Invalid --run_type '${params.run_type}'. Supported: ${VALID_RUN_TYPES.join(', ')}"
        }

        // 3. Validate species profile
        if (params.species && !VALID_SPECIES.contains(params.species.toLowerCase())) {
            errors << "Invalid --species '${params.species}'. Supported: ${VALID_SPECIES.join(', ')}"
        }

        // 4. Strict reference genome guard for reference_based runs
        if (params.run_type == 'reference_based') {
            if (!params.ref_nuc || params.ref_nuc == true || params.ref_nuc.toString().trim().isEmpty()) {
                errors << "--run_type 'reference_based' requires a valid reference genome FASTA path via --ref_nuc. A missing reference genome cannot produce a valid variant calling run."
            } else {
                def refStr = params.ref_nuc.toString().trim()
                def isRemote = refStr.startsWith('http://') || refStr.startsWith('https://') || refStr.startsWith('s3://')
                if (!isRemote) {
                    def refFile = new File(refStr)
                    if (!refFile.exists()) {
                        errors << "Reference FASTA file specified by --ref_nuc does not exist: '${refStr}'"
                    } else if (refFile.length() == 0) {
                        errors << "Reference FASTA file specified by --ref_nuc is empty (0 bytes): '${refStr}'"
                    }
                }
            }
        }

        // 5. Validate samplesheet path
        if (!params.input) {
            errors << "Parameter --input is missing. Please provide a path to a valid samplesheet CSV."
        } else {
            def samplesheetFile = new File(params.input.toString())
            if (!samplesheetFile.exists()) {
                errors << "Samplesheet file '${params.input}' does not exist!"
            } else {
                // 6. Deep samplesheet content validation
                validateSamplesheet(samplesheetFile, params.platform ?: 'auto', errors)
            }
        }

        if (errors) {
            def errorBanner = """
================================================================================
CRITICAL CONFIGURATION ERROR(S):
--------------------------------------------------------------------------------
${errors.collect { "  - " + it }.join('\n')}
================================================================================
Pipeline execution aborted to avoid producing incomplete or corrupted results.
"""
            log.error errorBanner
            System.exit(1)
        }
    }

    private static void validateSamplesheet(File file, String platform, List<String> errors) {
        def lines = file.readLines().findAll { it && !it.trim().startsWith('#') }
        if (lines.isEmpty()) {
            errors << "Samplesheet '${file.path}' is empty!"
            return
        }

        def header = lines[0].split(',').collect { it.trim().toLowerCase() }
        def requiredCols = ['run_id', 'input_path', 'input_type']
        for (col in requiredCols) {
            if (!header.contains(col)) {
                errors << "Samplesheet is missing required column '${col}'. Header contains: [${header.join(', ')}]"
                return
            }
        }

        def runIdIdx    = header.indexOf('run_id')
        def inputPathIdx = header.indexOf('input_path')
        def inputTypeIdx = header.indexOf('input_type')

        def seenSampleIds = new HashSet<String>()

        for (int i = 1; i < lines.size(); i++) {
            def rowNum = i + 1
            def cols = lines[i].split(',').collect { it.trim() }
            if (cols.size() <= Math.max(runIdIdx, Math.max(inputPathIdx, inputTypeIdx))) {
                errors << "Row ${rowNum} in samplesheet does not have enough columns: '${lines[i]}'"
                continue
            }

            def runId = cols[runIdIdx]
            def inputPath = cols[inputPathIdx]
            def inputType = cols[inputTypeIdx].toLowerCase()

            if (!runId) {
                errors << "Row ${rowNum}: 'run_id' is blank."
            } else {
                if (seenSampleIds.contains(runId)) {
                    errors << "Row ${rowNum}: Duplicate 'run_id' '${runId}' detected. Each sample/run must have a unique identifier."
                }
                seenSampleIds.add(runId)
            }

            if (!inputPath) {
                errors << "Row ${rowNum} (${runId}): 'input_path' is blank."
            } else {
                def isRemote = inputPath.startsWith('http://') || inputPath.startsWith('https://') || inputPath.startsWith('s3://')
                if (!isRemote) {
                    def targetFile = new File(inputPath)
                    if (!targetFile.exists()) {
                        errors << "Row ${rowNum} (${runId}): Input file does not exist at '${inputPath}'"
                    }
                }
            }

            // Platform and input type compatibility checks
            if (platform == 'ont') {
                if (!ONT_VALID_INPUT_TYPES.contains(inputType)) {
                    if (inputType in ['pacbio_bam', 'hifi', 'pacbio_hifi']) {
                        errors << "Row ${rowNum} (${runId}): PacBio format '${inputType}' supplied to ONT workflow (--platform ont). Use --platform pacbio or supply ONT data."
                    } else {
                        errors << "Row ${rowNum} (${runId}): Unsupported input_type '${inputType}' for ONT. Supported: ${ONT_VALID_INPUT_TYPES.join(', ')}"
                    }
                }
            } else if (platform == 'pacbio') {
                if (!PACBIO_VALID_INPUT_TYPES.contains(inputType)) {
                    if (inputType in ['pod5', 'fast5']) {
                        errors << "Row ${rowNum} (${runId}): ONT raw signal '${inputType}' supplied to PacBio workflow (--platform pacbio). ONT signal requires the ONT basecalling pipeline."
                    } else {
                        errors << "Row ${rowNum} (${runId}): Unsupported input_type '${inputType}' for PacBio. Supported: ${PACBIO_VALID_INPUT_TYPES.join(', ')}"
                    }
                }
            } else {
                // Auto / unified
                def allSupported = (ONT_VALID_INPUT_TYPES + PACBIO_VALID_INPUT_TYPES).unique()
                if (!allSupported.contains(inputType)) {
                    errors << "Row ${rowNum} (${runId}): Unsupported input_type '${inputType}'. Supported: ${allSupported.join(', ')}"
                }
            }
        }
    }

    public static void printSummary(workflow, params, log) {
        def refDisplay = params.ref_nuc ?: 'None (de novo assembly)'
        if (params.run_type == 'reference_based' && !params.ref_nuc) {
            refDisplay = 'MISSING (Error)'
        }
        def summary = """
--------------------------------------------------------------------------------
Execution Run Summary:
--------------------------------------------------------------------------------
Run Name     : ${workflow.runName}
Profile      : ${workflow.profile}
Species      : ${params.species}
Platform     : ${params.platform}
Run Type     : ${params.run_type}
Samplesheet  : ${params.input}
Output Dir   : ${params.outdir}
Reference    : ${refDisplay}
Phasing      : ${params.species == 'human' || params.enable_phasing ? 'Enabled' : 'Disabled'}
Container    : ${workflow.containerEngine ?: 'None'}
Launch Dir   : ${workflow.launchDir}
Work Dir     : ${workflow.workDir}
--------------------------------------------------------------------------------
"""
        log.info summary
    }

    public static void pipelineCompleted(wf, params, log) {
        if (!wf) return
        if (wf.success) {
            def cpuHours = wf.stats?.computeTimeFmt ?: (wf.stats?.cpuTime ? String.format("%.2f h", wf.stats.cpuTime / 3600000.0) : (wf.cpuHours ?: 'N/A'))
            log.info """
================================================================================
 Pipeline completed successfully!
 Duration : ${wf.duration ?: 'N/A'}
 CPU hours: ${cpuHours}
 Output   : ${wf.launchDir}/${params.outdir}
================================================================================
"""
        } else {
            log.error """
================================================================================
 Pipeline completed with ERRORS!
 Error message: ${wf.errorMessage ?: 'Unknown error'}
 Exit status  : ${wf.exitStatus ?: 'N/A'}
 Work dir     : ${wf.workDir ?: 'N/A'}
================================================================================
"""
        }
    }
}

