// JOB module

nextflow.enable.dsl = 2

//

process LEIDEN {

    tag "${sample}_${filtering}_${dimred}_${min_cell_number}_${k}"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering),
        val(dimred),
        val(min_cell_number),
        val(k)

    output:
    path "out_${sample}_${filtering}_${dimred}_${min_cell_number}_${k}.csv", emit: output
    path "log_${sample}_${filtering}_${dimred}_${min_cell_number}_${k}.txt", emit: logs
    
    script:
    """
    python ${baseDir}/bin/leiden_clustering.py \
        --min_cov_treshold ${params.min_cov_treshold} \
        --ncores ${task.cpus} \
        --p ${params.path_data} \
        --sample ${sample} \
        --filtering ${filtering} \
        --dimred ${dimred} \
        --min_cell_number ${min_cell_number} \
        --range ${params.resolution_range} \
        --k ${k} \
        --n ${params.n}
    """

    stub:
    """
    touch "out_${sample}_${filtering}_${dimred}_${min_cell_number}_${k}.csv"
    touch "log_${sample}_${filtering}_${dimred}_${min_cell_number}_${k}.txt"
    """

}