
# PPiSeq
code for second PPiSeq paper

## network properties

### pickle'd networkx files

These are networkx files, but pickled to store on disk.
To make the pickle'd networkx files, as well as graphml files to open in
Cytoscape and pajek (.net) files for possibly running in that Luvansec program,
run in a shell:

    find ~/Dropbox/PPiSeq_02/Paper_data/*mean_fitness_positive.csv | xargs -I '{}' sh -c 'python3 working_code/convert_csv_to_other_network_files.py --type mean {} ~/Dropbox/PPiSeq_02/Working_data/networks/$(basename {} .csv)'

(with networkx installed via `pip3`)


