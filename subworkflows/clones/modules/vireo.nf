// Vire module

nextflow.enable.dsl = 2

//

process VIREO {

    tag "${sample}_${filtering}_${min_cell_number}"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering),
        val(min_cell_number)

    output:
    path "out_vireo_${sample}_${filtering}_${min_cell_number}.pickle", emit: output
    path "log_vireo_${sample}_${filtering}_${min_cell_number}.txt", emit: logs
    
    script:
    """
    python ${baseDir}/bin/vireo.py \
    --min_cov_treshold ${params.min_cov_treshold} \
    --ncores ${task.cpus} \
    -p ${params.path_data} \
    --sample ${sample} \
    --filtering ${filtering} \
    --min_cell_number ${min_cell_number} \
    --range ${params.vireo_clone_range} \
    --p_treshold ${params.vireo_p_assignment} > scam.txt
    """

    stub:
    """
    touch "out_vireo_${sample}_${filtering}_${min_cell_number}.pickle"
    touch "log_vireo_${sample}_${filtering}_${min_cell_number}.txt"
    """

}


