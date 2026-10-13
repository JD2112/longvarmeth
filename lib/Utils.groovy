// lib/Utils.groovy
// Helper utilities for parameter validation, summary printing, and run completion

class Utils {

    public static void validateParameters(params, log) {
        if (!params.input) {
            log.error "ERROR: Parameter --input is missing. Please provide a path to a valid samplesheet CSV."
            System.exit(1)
        }
        def samplesheet = new File(params.input)
        if (!samplesheet.exists()) {
            log.error "ERROR: Samplesheet file '${params.input}' does not exist!"
            System.exit(1)
        }

        if (params.run_type == 'reference_based' && !params.ref_nuc) {
            log.warn "WARNING: --run_type reference_based is selected, but --ref_nuc is empty! Alignment and variant calling will be skipped."
        }
    }

    public static void printSummary(workflow, params, log) {
        def summary = """
--------------------------------------------------------------------------------
Execution Run Summary:
--------------------------------------------------------------------------------
Run Name     : ${workflow.runName}
Profile      : ${workflow.profile}
Species      : ${params.species}
Run Type     : ${params.run_type}
Samplesheet  : ${params.input}
Output Dir   : ${params.outdir}
Reference    : ${params.ref_nuc ?: 'None (de novo assembly)'}
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
            log.info """
================================================================================
 Pipeline completed successfully!
 Duration : ${wf.duration ?: 'N/A'}
 CPU hours: ${wf.cpuHours ?: 'N/A'}
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

