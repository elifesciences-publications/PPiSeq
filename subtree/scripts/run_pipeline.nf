#!/usr/bin/env nextflow

/////
/////
///// Pipeline for analyzing this PPISeq data, doing a few analyses
///// 
/////

// Making output directories
file("./tmp").mkdirs()
file("./reports").mkdirs()

//
//
// Arranging the PPISeq data into tables for later use
//
//

// compiling sqlite big long tables for counts
compile_data_counts = Channel.fromPath("scripts/compile_data_counts.ipynb")
counts_data = Channel
    .fromPath("data/Lineage_barcode_fitness_files/*.csv")
    .collect()
    .into{ counts_data_1; counts_data_2; }
process compile_counts {
    memory '8 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_python3_for_darach:v0.0.4"
    input:
        file(script) from compile_data_counts
        file counts_data_1
    output:
        file("counts.sqlite") into counts_database
        file("*.html") into counts_report
    shell:
    '''
    jupyter nbconvert --to=html --ExecutePreprocessor.timeout=-1 --execute !{script}
    '''
}
counts_database
    .into{
        counts_database_1;
        counts_database_2;
        counts_database_3;
        counts_database_4;
        }

// compiling sqlite per barcoded lineage fitness
compile_data_fitness_per_lineage = Channel.fromPath("scripts/compile_data_fitness_per_lineage.ipynb")
process compile_fitness_per_lineage {
    memory '8 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_python3_for_darach:v0.0.4"
    input:
        file(script) from compile_data_fitness_per_lineage
        file counts_data_2
    output:
        file("lineage_fitness.sqlite") into lineage_fitness_database
        file("*.html") into compile_fitness_per_lineage_report
    shell:
    '''
    jupyter nbconvert --to=html --ExecutePreprocessor.timeout=-1 --execute !{script}
    '''
}

// summarizing counts
counts_summarize = Channel.fromPath("scripts/counts_summarize.ipynb")
process summarizing_counts {
    memory '2 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_python3_for_darach:v0.0.4"
    input:
        file script from counts_summarize
        file counts_database_3
    output:
        file("*.html") into reports_summarizing_counts
        file("counts_summary.sqlite") into counts_summary_database
    shell:
    '''
    jupyter nbconvert --to=html --ExecutePreprocessor.timeout=-1 --execute !{script}
    '''
}
counts_summary_database
    .into{
        counts_summary_database_1;
        }

// compiling sqlite big long tables for normalized mean fitness
compile_data_fitness = Channel.fromPath("scripts/compile_data_fitness.ipynb")
fitness_data_pre = Channel
    .fromPath("data/PPI_mean_fitness_calling_files/*.csv")
    .collect()
    .into{ fitness_data; }
process compile_fitness {
    memory '2 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_python3_for_darach:v0.0.4"
    input:
        file(script) from compile_data_fitness
        file fitness_data
    output:
        file("fitness.sqlite") into fitness_database
        file("*.html") into fitness_report
    shell:
    '''
    jupyter nbconvert --to=html --ExecutePreprocessor.timeout=-1 --execute !{script}
    '''
}

fitness_database
    .into{
        fitness_database_1;
        fitness_database_2;
        fitness_database_3;
        fitness_database_4;
        fitness_database_5;
        fitness_database_6;
        fitness_database_7;
        fitness_database_8;
        }





//
//
// Compiling together features of ORFs
//
//

// naming and GO resources, list
sgd_go = Channel.fromPath("data/other_sources/gene_association_191111.sgd")
go_terms = Channel.fromPath("data/other_sources/go_terms_191111.tab")
go_slim = Channel.fromPath("data/other_sources/go_slim_mapping_191111.tab")
go_complex_slim = Channel.fromPath("data/other_sources/go_protein_complex_slim_191111.tab")

compile_data_sgd_go = Channel.fromPath("scripts/compile_data_sgd_go.R")
process compile_data_sgd_go {
    publishDir 'tmp'
    memory '2 GB'
    container "shub://darachm/singularity_r_for_darach:v0.0.15"
    input:
        file(script) from compile_data_sgd_go
        file sgd_go
        file go_terms
        file go_slim
        file go_complex_slim
    output:
        file("sgd_go.sqlite") into sgd_go_data
        file("*.html") into sgd_go_report
    shell:
    '''
    Rscript -e "knitr::stitch_rhtml('!{script}')"
    '''
}

sgd_go_data
    .into{
        sgd_go_data_1;
        sgd_go_data_2;
        sgd_go_data_3;
        sgd_go_data_4;
        }

////
////
//// chong2015 localization data preprocess
////
////

Channel.fromPath("data/other_sources/chongEtAl2015/wt*.csv")
    .collect()
    .into{ chong2015_localization; }
homodimer_localization_pre = Channel.fromPath("scripts/localization_preprocess.R")
process localization_pre_homodimers {
    memory '2 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_r_for_darach:v0.0.15"
    input:
        file(script) from homodimer_localization_pre
        file chong2015_localization 
    output:
        file("localization_data.sqlite") into localization_data
        file("*.html") into localization_report
    shell:
    '''
    Rscript -e "knitr::stitch_rhtml('!{script}')"
    '''
}

localization_data
    .into {
        localization_data_1 ;
        localization_data_2 ;
        }

////
////
//// compiling some tables into sqlite
////
////

Channel.fromPath("data/other_sources/geneFeatures_022415_EK.txt")
    .into{ gene_features }
Channel.fromPath("data/other_sources/protein_properties_190922.tab")
    .into{ protein_properties }
Channel.fromPath("data/other_sources/ho2018unified_ProteinAbundances.csv")
    .into{ ho_abundances;  }
Channel.fromPath("data/other_sources/BIOGRID-PTM-3.5.177_ptm_yeast.txt")
    .into{ biogrid_ptm; }
// Generated by:
// wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.177/BIOGRID-PTMS-3.5.177.ptm.zip 
// unzip BIOGRID-PTMS-3.5.177.ptm.zip 
// head -n 1 BIOGRID-PTM-3.5.177.ptmtab.txt > BIOGRID-PTM-3.5.177_ptm_yeast.txt
// grep " 559292  " BIOGRID-PTM-3.5.177.ptmtab.txt >> BIOGRID-PTM-3.5.177_ptm_yeast.txt
Channel.fromPath("data/other_sources/sgd_homologues_191027.csv")
    .into{ homologues }
Channel.fromPath("data/other_sources/ygob_ohnologs_191028.csv")
    .into{ ohnologs }
Channel.fromPath("data/other_sources/marchant2019_tableS3.csv")
    .into{ marchant_ohnologs }
compile_data_features = Channel.fromPath("scripts/compile_data_features.R")
process datasets_compiling {
    memory '2 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_r_for_darach:v0.0.15"
    input:
        file script from compile_data_features
        file protein_properties
        file ho_abundances
        file biogrid_ptm
        file gene_features
        file homologues
        file ohnologs
        file marchant_ohnologs
    output:
        file("features.sqlite") into datasets_feature_tables
        file("*.html") into data_features_report
    shell:
    '''
    Rscript -e "knitr::stitch_rhtml('!{script}')"
    '''
}
datasets_feature_tables
    .into {
        datasets_feature_tables_1;
        datasets_feature_tables_2;
        datasets_feature_tables_3;
        datasets_feature_tables_4;
        datasets_feature_tables_5;
        }

// per orf database
compile_orf_features = Channel.fromPath("scripts/compile_orf_features.ipynb")
process orf_features_compile {
    memory '8 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_python3_for_darach:v0.0.4"
    input:
        file script from compile_orf_features
        file fitness_database_1
        file datasets_feature_tables_1
        file sgd_go_data_1
    output:
//        file("*.html") into reports_orf_features
        file("orf_features.sqlite") into orf_features_database
    shell:
    '''
    jupyter nbconvert --to=html --ExecutePreprocessor.timeout=-1 --execute !{script}
    '''
}

orf_features_database
    .into {
        orf_features_database_1;
        orf_features_database_2;
        orf_features_database_3;
        orf_features_database_4;
        orf_features_database_5;
        }








//
//
// Homodimers analyses
//
//

strict_thresholds_2 = Channel.fromPath("data/precomp/PPI_environment_count_summary_SD_merge_filter_strict_threshold.csv")

// homodimers modeling
homodimer_modeling = Channel.fromPath("scripts/homodimer_modeling_fitting.R")
process model_homodimers {
    memory '4 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_r_for_darach:v0.0.15"
    input:
        file script from homodimer_modeling
        file fitness_database_3
        file strict_thresholds_2
    output:
        file("*.html") into reports_modeling
        file("homodimer_all_models.csv") into homodimer_dataframe
    shell:
    '''
    Rscript -e "knitr::stitch_rhtml('!{script}')"
    '''
}

// homodimers interpreting explained
homodimer_interpreting = Channel.fromPath("scripts/homodimer_modeling_interpreting.R")
process interpret_homodimers {
    memory '2 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_r_for_darach:v0.0.15"
    input:
        file script from homodimer_interpreting
        file homodimer_dataframe
    output:
        file("*.html") into reports_interpreting
        set file("explained_ppis_modeling.csv"), 
            file("homodimer_all_models.csv") into homodimer_results
    shell:
    '''
    Rscript -e "knitr::stitch_rhtml('!{script}')"
    '''
}

homodimer_results
    .into{
        homodimer_results_1;
        homodimer_results_2;
        homodimer_results_3;
        homodimer_results_4;
        homodimer_results_5;
        homodimer_results_6;
    }

// making homodimer database
compile_homodimer_ppi_features = Channel.fromPath("scripts/compile_homodimer_ppi_features.ipynb")
process homodimer_ppi_features_compile {
    memory '4 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_python3_for_darach:v0.0.4"
    input:
        file script from compile_homodimer_ppi_features
        file homodimer_results_1
        file datasets_feature_tables_2
        file localization_data_1
        file orf_features_database_1
    output:
        file("*.html") into reports_homodimer_ppi_features
        file("homodimer_ppi_features.sqlite") into homodimer_ppi_features_database
    shell:
    '''
    jupyter nbconvert --to=html --ExecutePreprocessor.timeout=-1 --execute !{script}
    '''
}

// Species accumulation analysis
accumulation_curves_crunch = Channel.fromPath("scripts/accumulation_curves.R")
validation_rates = Channel.fromPath("data/precomp/Validated_PPI_environment_count_summary_SD_merge_filter.csv")
process crunch_accumulation_curves {
    cpus '7'
    memory '8 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_r_for_darach:v0.0.15"
    input:
        file script from accumulation_curves_crunch
        file fitness_database_7
        file validation_rates
    output:
        set file("*.RData"),
            file("accumulation_curves.html") into crunch_accumulation_curves
    shell:
    '''
    curdir=$(pwd)
    command="rmarkdown::render(input='!{script}',output_dir='$curdir',knit_root_dir='$curdir')"
    Rscript -e $command
    '''
}

accumulation_curves_report = Channel.fromPath("scripts/accumulation_curves_report.R")
process report_accumulation_curves {
    cpus '1'
    memory '2 GB'
    publishDir 'tmp'
    container "shub://darachm/singularity_r_for_darach:v0.0.15"
    input:
        file script from accumulation_curves_report
        set file(rdata), file(html) from crunch_accumulation_curves
    output:
        set file("accumulation_curves_report.html"), 
            file("*png"), file("*pdf") into report_accumulation_curves
    shell:
    '''
    curdir=$(pwd)
    command="rmarkdown::render(input='!{script}',output_dir='$curdir',knit_root_dir='$curdir')"
    Rscript -e $command
    '''
}



//
//
// analyze and conclude homodimer
//
//

strict_thresholds = Channel.fromPath("data/precomp/PPI_environment_count_summary_SD_merge_filter_strict_threshold.csv")

homodimer_report = Channel.fromPath("scripts/homodimer_report.Rmd")
process report_homodimer {
    memory 4G
    publishDir 'tmp'
    container "shub://darachm/singularity_r_for_darach:v0.0.15"
    input:
        file script from homodimer_report
        file fitness_database_4
        file homodimer_results_4
        file homodimer_ppi_features_database
        file datasets_feature_tables_3
        file sgd_go_data_2
        file orf_features_database_3
        file strict_thresholds
        set file(report), file(plotz), file(plotz2) from report_accumulation_curves
    output:
        set file("homodimer_report.html"), file("*png"), file("*pdf") into report_homodimer
    shell:
    '''
    curdir=$(pwd)
    command="rmarkdown::render(input='!{script}',output_dir='$curdir',knit_root_dir='$curdir')"
    Rscript -e $command
    '''
}



