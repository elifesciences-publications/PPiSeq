.PHONY: all setup \
	tecan_analysis \
	prcomp_analysis \
	networkx_pickles

input_data_csv = $(shell find data -maxdepth 1 -name "*PPI*.csv")
networkx_pickles = $(subst data,tmp,$(subst _PPI_barcodes_fitness_counts.csv,.nxp,$(input_data_csv)))

tecan_analysis: output/validations_statistic_against_properties.png
prcomp_analysis: 
networkx_pickles: $(networkx_pickles)

all: setup tecan_analysis prcomp_analysis

setup:
	mkdir -p tmp
	mkdir -p output

rscript = singularity exec shub://darachm/singularity_r_for_darach:v0.0.2 Rscript
python3 = singularity exec shub://darachm/singularity_python3_for_darach:latest python3

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
		$(networkx_pickles) 
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

tmp/%.nxp: data/%_PPI_barcodes_fitness_counts.csv \
		scripts/convert_csv_to_networkx_graph.py
	$(python3) $(word 2,$^) $(word 1,$^) $@
