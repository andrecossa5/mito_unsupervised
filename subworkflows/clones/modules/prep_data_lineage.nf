// SCORE_LINEAGE module

nextflow.enable.dsl = 2

//

process PREP_LINEAGE {

    tag "${sample}"

    input:
    tuple val(sample), val(min_cell_number)

    output:
    tuple val(sample), val(min_cell_number), path("formatted.csv"), path("meta.csv"), emit: complete_matrix
    
    script:
    """
    python ${baseDir}/bin/prep_data_for_lineage.py 
    ${params.path_data} \
    ${sample} \
    ${min_cell_number}
    """

    stub:
    """
    touch formatted.csv
    touch meta.csv
    """

}