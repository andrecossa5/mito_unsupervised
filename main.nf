// mito_unsupervised pipeline
nextflow.enable.dsl = 2
include { infer_clones } from "./subworkflows/clones/main"

// Samples channel
ch_samples = Channel
    .fromPath("${params.path_data}/*", type:'dir') 
    .map{ it.getName() }

//

//----------------------------------------------------------------------------//
// mito_unsupervised entry points
//----------------------------------------------------------------------------//

//

workflow clones {

    infer_clones(ch_samples)
    infer_clones.out.job_output.view()

}

//

// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}