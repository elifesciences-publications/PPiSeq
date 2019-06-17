.PHONY: all setup \
	tecan_analysis \
	prcomp_analysis \
	networks

input_data_csv = $(shell find data -maxdepth 1 -name "*PPI*.csv")
mean_csvs = $(shell find data -maxdepth 1 -name "*_mean_fitness_positive.csv")

tecan_analysis: output/validations_statistic_against_properties.png
prcomp_analysis: 

networks_mean_pickles = $(subst data,tmp,$(subst _PPI_barcodes_fitness_counts.csv,_mean.nxp,$(input_data_csv)))
networks_mean_graphml = $(subst data,tmp,$(subst _PPI_barcodes_fitness_counts.csv,_mean.graphml,$(input_data_csv)))
networks_mean_pajek = $(subst data,tmp,$(subst _PPI_barcodes_fitness_counts.csv,_mean.net,$(input_data_csv)))

networks: $(networks_mean_pickles) $(networks_raw_graphml) $(networks_raw_pajek)

all: setup tecan_analysis prcomp_analysis

setup:
	mkdir -p tmp
	mkdir -p output

rscript = singularity exec shub://darachm/singularity_r_for_darach:v0.0.3 Rscript
python3 = singularity exec shub://darachm/singularity_python3_for_darach:latest python3

# RData object for moving between names easily
tmp/NameIDList.RData: scripts/makingLists.R data/SGD_features_170601.tab
	$(rscript) $<

# Setting up RData for hastening subsequent analysis
tmp/all_environments_counts.RData: \
		scripts/reading_counts_to_rdata.R \
		$(input_data_csvs) 
	$(rscript) $<
tmp/all_environments_scores_and_stats.RData: \
		scripts/reading_stats_to_rdata.R \
		$(mean_csvs) 
	$(rscript) $<

#tmp/all_environments_all_biological_ppi_mean_fitness.RData: \
#		scripts/averaging_all_fitnesses.R \
#		tmp/all_environments_counts.RData
#	$(rscript) $<

# Analyzing tecan validations
tmp/tecan_validation_statistics.csv: \
		scripts/calculate_AUC.R \
		data/tecan_validation_assays/
	$(rscript) $<
output/validations_statistic_against_properties.png: \
		scripts/analyze_tecan.R \
		tmp/all_environments_counts.RData \
		tmp/tecan_validation_statistics.csv
	$(rscript) $<

# Principal components analyses
output/prcomp_of_samples_pc12345678.svg: \
		scripts/prcomp_analysis.R \
		tmp/all_environments_counts.RData \
		tmp/NameIDList.RData
	$(rscript) $< 

#
#
#
#
#

tmp/%_mean.nxp: scripts/convert_csv_to_other_network_files.py \
		data/%_mean_fitness_positive.csv
	$(python3) $(word 1,$^) $(word 2,$^) $(subst .nxp,,$@) --type mean
	# type is mean, as opposed to making networks of the raw multigraph 
	# 	(each barcode data independently)
tmp/%_mean.graphml: scripts/convert_csv_to_other_network_files.py \
		data/%_mean_fitness_positive.csv
	$(python3) $(word 1,$^) $(word 2,$^) $(subst .graphml,,$@) --type mean
tmp/%_mean.net: scripts/convert_csv_to_other_network_files.py \
		data/%_mean_fitness_positive.csv
	$(python3) $(word 1,$^) $(word 2,$^) $(subst .net,,$@) --type mean

#
#
#
#
#

tmp/condition_multigraph.nxp: scripts/make_conditions_multigraph.py \
		$(networks_mean_pickles)
	$(python3) $^ --output_base $(subst .nxp,,$@)
tmp/condition_multigraph.net: scripts/make_conditions_multigraph.py \
		$(networks_mean_pickles)
	$(python3) $^ --output_base $(subst .net,,$@)


tmp/ppi_entropy_per_protein.txt: scripts/make_conditions_multigraph.py tmp/condition_multigraph.nxp
	$(python3) $^ 
