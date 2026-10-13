#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { LONGVARMETH } from './workflows/longvarmeth'
include { NANOPORE_PIPELINE } from './workflows/nanopore'
include { PACBIO_PIPELINE } from './workflows/pacbio'

def logo = """
================================================================================
  _                                                    _   _     
 | | ___  _ __   __ _  __   ____ _ _ __ _ __ ___   ___| |_| |__  
 | |/ _ \\| '_ \\ / _` | \\ \\ / / _` | '__| '_ ` _ \\ / _ \\ __| '_ \\ 
 | | (_) | | | | (_| |  \\ V / (_| | |  | | | | | |  __/ |_| | | |
 |_|\\___/|_| |_|\\__, |   \\_/ \\__,_|_|  |_| |_| |_|\\___|\\__|_| |_|
                |___/                                            
================================================================================
  longvarmeth: long-read variant calling & methylation pipeline
  version : ${workflow.manifest.version ?: '1.0.0'}
  species : ${params.species} | mode : ${params.run_type}
================================================================================
"""

log.info logo


workflow {
    Utils.validateParameters(params, log)
    Utils.printSummary(workflow, params, log)

    if (params.platform == 'pacbio') {
        PACBIO_PIPELINE()
    } else if (params.platform == 'ont') {
        NANOPORE_PIPELINE()
    } else {
        LONGVARMETH()
    }
}

workflow.onComplete {
    Utils.pipelineCompleted(workflow, params, log)
}


