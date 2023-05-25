// mito_unsupervised

// Include here
nextflow.enable.dsl = 2
include { LEIDEN } from "./modules/leiden.nf"
include { PREP_LINEAGE } from "./modules/prep_data_lineage.nf"
include { LINEAGE } from "./modules/lineage.nf"
include { SCORE_LINEAGE } from "./modules/score_lineage.nf"
include { VIREO } from "./modules/vireo.nf"

// 
 
//----------------------------------------------------------------------------//
// clustering_clones subworkflow
//----------------------------------------------------------------------------//

workflow infer_clones {
    
    take:
        ch_samples 

    main:

        if (params.leiden) {

            options = ch_samples
                .combine(params.filtering)
                .combine(params.dimred)
                .combine(params.min_cell_number)
                .combine(params.k)
                .filter{ !(it[1] != "pegasus" && it[2] != "no_dimred") }
            LEIDEN(options)
            results = LEIDEN.out.output

        } 
        if (params.lineage) {
            
            options = ch_samples
                .combine(params.min_cell_number)
            PREP_LINEAGE(options)
            LINEAGE(PREP_LINEAGE.out.complete_matrix)
            SCORE_LINEAGE(LINEAGE.out.labels)
            results = SCORE_LINEAGE.out.output

        }
        if (params.vireo) {
            
            options = ch_samples
                .combine(params.filtering)
                .combine(params.min_cell_number)
                .filter{ it[1] == 'MQuad' }
            VIREO(options)
            results = VIREO.out.output

        }

    emit:
        job_output = results

}