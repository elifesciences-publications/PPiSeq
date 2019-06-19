#!/usr/bin/env nextflow

/////
///// Pipeline for testing how informative different insert sizes will be for
///// the barcode-ORF pairing reaction.
/////

file("./tmp").mkdirs()
file("./reports").mkdirs()

//
//
//
//
//

ppis_mean_signal = Channel.fromPath("data/*mean_fitness_positive.csv")
convert_csv_to_nxp = Channel.fromPath("scripts/convert_csv_to_other_network_files.py")

process csv_to_networkx_pickles {
    publishDir 'tmp'
    container "shub://darachm/singularity_python3_for_darach:latest"
    cpus 1
    memory 2G
    input: 
        each file(script) from convert_csv_to_nxp
        file(csv) from ppis_mean_signal
    output: 
        file("mean.nxp") into mean_nxp
        file("mean.graphml") into mean_graphml
    shell:
    '''
    python3 !{script} --type mean !{csv} mean
    '''
}

//
//
//
//
//


tecan_plate_map = Channel.fromPath("data/tecan_validation_assays/plate_strain_map/*")
    .collect()
tecan_data_folders = Channel.fromPath("data/tecan_validation_assays/tecan_data/*",type:"dir")
    .collect()
calculate_auc = Channel.fromPath("scripts/calculate_AUC.R")


process calculate_tecan_validation_statistic {
    publishDir 'tmp'
    container "shub://darachm/singularity_r_for_darach:v0.0.3"
    cpus 1
    memory 2G
    input:
        each file(script) from calculate_auc
        file('*') from tecan_data_folders 
        file('*') from tecan_plate_map
    output:
        file("tecan_validation_statistics.csv") into tecan_validation_csv
    shell:
    '''
    Rscript !{script}
    '''
}


// container "shub://darachm/singularity_runningJobs:v0.2.0"


// make networks that are 2 or more, 3 or more environments
// make networks that are subsets of these, but delta to DMSO
// calc entropy/specificity
// validation data and predictions

// vim: syntax=javascript
