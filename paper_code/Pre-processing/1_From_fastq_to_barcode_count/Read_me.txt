For most sequencing, we pooled different samples into one HiSeq lane. 

"PPiSeq_sequencing_data_annotation.xlsx" contains the information that which samples are contained in each Hiseq data.

"PPiSea_condition_primers.xlsx" contains the unique tags for each sample, which can be used to demultiplex each sample.

"Pipeline.pdf" is a figure to describe the pipeline (details described in the methods part of the paper).

"BC.py" is a home-made script to extract DNA barcode and count them by clustering, but here we use it to demultiplex the pooled samples and then extract the barcode (lintag1 or lintag2) and UMI (seqtag1 or seqtag2).

"Parse one HiSeq data of Hydroxyurea environment.ipynb" is a notebook to run BC.py to parse one Hiseq data of Hydroxyurea environment

After obtaining barcodes (lintag) at each time point in each environment, we merge barcodes at different time points for each environment. Then we run Bartender (Bartender Clustering, https://github.com/LaoZZZZZ/bartender-1.1) to cluster barcodes, obtaining cluster centroids and reads belonging to each cluster. Note that the input of Bartender Clustering (bartender_single_com) needs two columns: barcodes and unique molecular identifier (UMI). Here we use line number to replace UMI. Basically we consider there is no PCR duplicates at this step.

"Match_cluster_known_barcodes_levenshtein_update.py" is a home-made script to map the cluster centroid (bartender output) to the known barcodes, and then count the number of unique UMIs for each double barcode. 

"Map bartender cluster result to known barcodes in Hydroxyurea environment.ipynb" is a notebook to run Match_cluster_known_barcodes_levenshtein_update.py to obtain the final counts for known barcodes at each time point in Hydroxyyrea environment.

We go through the above procedures to obtain barcode counts at different time points for each environment.






