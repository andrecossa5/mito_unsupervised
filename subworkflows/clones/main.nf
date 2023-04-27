// measter subworkflow

// Include here
nextflow.enable.dsl = 2
include { LEIDEN } from "./modules/leiden.nf"

// 

//----------------------------------------------------------------------------//
// clustering_clones subworkflow
//----------------------------------------------------------------------------//

workflow clustering_clones {
    
    take:
        ch_samples 

    main:

        // Create options
        options = ch_samples
        .combine(params.filtering)
        .combine(params.dimred)
        .combine(params.min_cell_number)
        .combine(params.k)
        .filter{ !(it[1] != "pegasus" && it[3] != "no_dimred") }

        // Execute jobs
        LEIDEN(options)

    emit:
        job_output = LEIDEN.out.logs

}