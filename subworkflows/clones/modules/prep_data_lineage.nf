// SCORE_LINEAGE module

nextflow.enable.dsl = 2

//

process PREP_LINEAGE {

    tag "${sample}"

    input:
    val(sample)

    output:
    tuple val(sample), path("formatted.csv"), path("meta.csv"), emit: complete_matrix
    
    script:
    """
    python ${baseDir}/bin/prep_data_for_lineage.py ${params.path_data} ${sample}
    """

    stub:
    """
    touch formatted.csv
    touch meta.csv
    """

}