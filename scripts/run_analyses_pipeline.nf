#!/usr/bin/env nextflow

/////
///// Pipeline for testing how informative different insert sizes will be for
///// the barcode-ORF pairing reaction.
/////

file("./tmp").mkdirs()
file("./reports").mkdirs()

ppis_mean_signal = Channel.fromPath("data/*mean_fitness_positive.csv")
convert_csv_to_nxp = Channel.fromPath("scripts/convert_csv_to_other_network_files.py")

process csv_to_networkx_pickles {
    publishDir 'tmp'
    container "shub://darachm/singularity_python3_for_darach:latest"
    cpus 1
    memory 2G
    input: 
        file(csv) from ppis_mean_signal
        file(script) from convert_csv_to_nxp
    output: 
        file("mean.nxp") into mean_nxp
        file("mean.graphml") into mean_graphml
    shell:
    '''
    python3 !{script} --type mean !{csv} mean
    '''
}

// container "shub://darachm/singularity_runningJobs:v0.2.0"


