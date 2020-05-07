“1_distance_filter_script.R” is a R script to filter out barcodes with high fitness measurement (Fit-seq) errors. 

“2_sort_PPI_barcode_count_4_timepoints.py” or “2_sort_PPI_barcode_count.py” is a python script to sort barcodes based on PPIs they marked. Basically, barcodes that label the same PPI will be in consecutive rows. “2_sort_PPI.ipynb” is the notebook to run the above two scripts in each environment.

“3_p_value_calculation_script.R” is a R script to calculate a p-value for each PPI based on its replicate fitness measurements in each environment.

“4_reference_generation_script.R” is a R script to generate one positive and 50 random reference sets of PPIs in each environment.

“5_threshold_script.R” is a R script to characterize various dynamic thresholds to call positive PPIs in each environment. 

“6_PPI_calling_script.R” is a R script to call positive PPIs with the optimal dynamic threshold in each environment.

“7_merge_SD_call_PPIs.R” is a R script to merge two SD datasets and re-call positive PPIs with the new merged SD datasets (re-run the above procedures).

“8_transform_data_summarize_filter_promiscuous_script.R” is a R script to  summarize promiscuous proteins detected in each environment, and based on that further filter out some false positive PPIs. Meanwhile, this script also generate summary tables that contain comprehensive information for each environment. 

“9_count_PPI_environment_script.R” is a R script to generate summary for PPI calling in each environment. The output file (“PPI_environment_count_summary_SD_merge_filter.csv”) is a commonly used table to generate several figures.

“10_normalize_combine_fitness_variation_score_script.R” is a R script to normalize fitness values for each environment, average fitness values for the same PPI, and calculate the mutability score for each PPI.  The output file (“Variation_score_PPI_environment_pos_SD_merge_filter.csv”) is commonly used to generate several figures. 

“11_strict_threshold.R” is a R script to choose a stricter threshold in each environment, and re-run the above procedures (6-10) to generate two commonly used tables in 9 and 10. 
