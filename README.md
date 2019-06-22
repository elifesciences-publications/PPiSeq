
# PPiSeq
code for second PPiSeq paper

## network properties

### pickle'd networkx files

These are networkx files, but pickled to store on disk.
To make the pickle'd networkx files, as well as graphml files to open in
Cytoscape and pajek (.net) files for possibly running in that Luvansec program,
run in a shell:

    mkdir ~/Dropbox/PPiSeq_02/Working_data/networks
    find ~/Dropbox/PPiSeq_02/Paper_data/*mean_fitness_positive.csv | xargs -I '{}' sh -c 'python3 working_code/convert_csv_to_other_network_files.py --type mean {} ~/Dropbox/PPiSeq_02/Working_data/networks/$(basename {} .csv)'

(with networkx installed via `pip3`)

### calculate network properties, make hairball network displays

Run in shell:

    find ~/Dropbox/PPiSeq_02/Working_data/networks/*.nxp | xargs -I '{}' sh -c 'python3 working_code/analyzing_network_properties.py {} --figure_base ~/Dropbox/PPiSeq_02/Working_data/networks/$(basename {} .nxp) --stats_base ~/Dropbox/PPiSeq_02/Working_data/networks/$(basename {} .nxp)'

### make multigraph across all conditions

Run in shell:

    python3 working_code/make_conditions_multigraph.py ~/Dropbox/PPiSeq_02/Working_data/networks/*mean_fitness_positive.nxp --output_base ~/Dropbox/PPiSeq_02/Working_data/networks/conditions_multigraph
