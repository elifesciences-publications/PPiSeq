.PHONY: all setup \
	tecan_analysis \
	prcomp_analysis 

all: setup tecan_analysis prcomp_analysis

tecan_analysis: output/validations_statistic_against_properties.png
prcomp_analysis: 

setup:
	mkdir -p tmp
	mkdir -p output

rscript=singularity exec shub://darachm/singularity_r_for_darach:v.0.0.2

# RData object for moving between names easily
tmp/NameIDList.RData: scripts/makingLists.R data/SGD_features_170601.tab
	$(rscript) $<

# Setting up RData for hastening subsequent analysis
tmp/all_environments_counts.RData: \
		scripts/reading_counts_to_rdata.R \
		data/*PPI_barcodes_fitness*.csv
	$(rscript) $<
tmp/all_environments_scores_and_stats.RData: \
		scripts/reading_stats_to_rdata.R \
		data/*mean_fitness*.csv
	$(rscript) $<
tmp/all_environments_all_biological_ppi_mean_fitness.RData: \
		scripts/averaging_all_fitnesses.R \
		tmp/all_environments_counts.RData
	$(rscript) $<

# Analyzing tecan validations
tmp/tecan_validation_statistics.csv: \
		scripts/calculate_AUC.R \
		data/tecan_validation_assays/
	$(rscript) $<
output/validations_statistic_against_properties.png: \
		scripts/analyze_tecan.R \
		tmp/all_environments_scores_and_stats.RData tmp/tecan_validation_statistics.csv
	$(rscript) $<

# Principal components analyses
output/prcomp_of_samples_pc12345678.svg: \
		scripts/prcomp_analysis.R \
		tmp/all_environments_counts.RData \
		tmp/NameIDList.RData
	$(rscript) $< 
