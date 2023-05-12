// SCORE_LINEAGE module

nextflow.enable.dsl = 2

//

process SCORE_LINEAGE {

    tag "${sample}"
    publishDir "${params.outdir_lineage}/", mode: 'copy'

    input:
    tuple val(sample), path(ground_truth), path(lineage)

    output:
    path "out_${sample}.csv", emit: output
    path "log_${sample}.txt", emit: logs
    
    script:
    """
    python ${baseDir}/bin/score_lineage.py \
    ${sample} \
    ${ground_truth} \
    ${lineage} > log_${sample}.txt
    """

    stub:
    """
    touch "out_${sample}.csv"
    touch "log_${sample}.txt"
    """

}