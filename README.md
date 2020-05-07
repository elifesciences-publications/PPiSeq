
# PPiSeq
## Codes for the analyses of PPiSeq data
### Analyses for Figure 1, 2, 3, 4A-C, 5, and Supplementary Figures S1-9
There is a script for each figure panel (e.g. Figure1A). To run the codes please download the addtional datasets (Mendeley Data, V1, doi: 10.17632/9ygwhk5cs3.1), and change the directory in the script accordingly. You can find the explanations for some important datasets through Read_me.txt in the corresponding directory. 
### Analyses for Figure 4D-F, 6, and Supplemenatry Figure S10
(Darach edit)


## network properties

### pickle'd networkx files

These are networkx files, but pickled to store on disk.
To make the pickle'd networkx files, as well as graphml files to open in
Cytoscape and pajek (.net) files for possibly running in that Luvansec program,
run in a shell:

    mkdir ~/Dropbox/PPiSeq_02/Working_data/networks
    find ~/Dropbox/PPiSeq_02/Paper_data/*mean_fitness_positive.csv | xargs -I '{}' sh -c 'python3 working_code/convert_csv_to_other_network_files.py --type mean {} ~/Dropbox/PPiSeq_02/Working_data/networks/$(basename {} .csv)'

(with networkx installed via `pip3`)

### calculate network properties, write json of topological properties

Run in shell:

    find ~/Dropbox/PPiSeq_02/Working_data/networks/*.nxp | xargs -I '{}' sh -c 'python3 working_code/analyzing_network_properties.py {} --figure_base ~/Dropbox/PPiSeq_02/Working_data/networks/ --stats_base ~/Dropbox/PPiSeq_02/Working_data/networks/'


### make multigraph across all conditions

Run in shell:

    python3 working_code/make_conditions_multigraph.py ~/Dropbox/PPiSeq_02/Working_data/networks/*mean_fitness_positive.nxp --output_base ~/Dropbox/PPiSeq_02/Working_data/networks/conditions_multigraph

### calculate variability for nodes across all conditions, all partners

Run in shell:

    python3 working_code/calculating_variability.py ~/Dropbox/PPiSeq_02/Working_data/networks/conditions_multigraph.nxp --output_base ~/Dropbox/PPiSeq_02/Working_data/networks/multigraph

### model how consistency of linkage patterns vary with respect to gene features:

Run in shell:

    Rscript 'rmarkdown::render("working_code/variability_modeling.Rmd")'
