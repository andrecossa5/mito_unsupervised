// SCORE_LINEAGE module

nextflow.enable.dsl = 2

//

process SCORE_LINEAGE {

    tag "${sample}_${min_cell_number}"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    tuple val(sample), val(min_cell_number), path(ground_truth), path(lineage)

    output:
    path "out_LINEAGE_${sample}_${min_cell_number}.pickle", emit: output
    path "log_LINEAGE_${sample}_${min_cell_number}.txt", emit: logs
    
    script:
    """
    python ${baseDir}/bin/score_lineage.py \
    ${sample} \
    ${ground_truth} \
    ${lineage} \
    ${min_cell_number} > log_LINEAGE_${sample}_${min_cell_number}.txt
    """

    stub:
    """
    touch "out_LINEAGE_${sample}_${min_cell_number}.pickle"
    touch "log_LINEAGE_${sample}_${min_cell_number}.txt"
    """

}