// LINEAGE module

nextflow.enable.dsl = 2

//

process LINEAGE {

    tag "${sample}"

    input:
    tuple val(sample), val(min_cell_number), path(afm), path(meta)

    output:
    tuple val(sample), val(min_cell_number), path("ground_truth.csv"), path("lineage_labels.csv"), emit: labels
    
    script:
    """
    Rscript ${baseDir}/bin/lineage_clustering.r \
        ${sample} \
        ${afm} \
        ${meta} \
        ${params.lineage_iters} \
        ${task.cpus}
    """

    stub:
    """
    touch lineage_labels.csv
    touch ground_truth.csv
    """

}