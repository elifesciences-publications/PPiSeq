#!/usr/bin/env python3

import argparse
import csv
import networkx as nx
import pickle
import re
import json
import os

parser = argparse.ArgumentParser()
parser.add_argument("network_pickles",nargs="+",
    help="give me .nxp files"
    )
parser.add_argument("--output_base",type=str,
    help="output basename"
    )
args = parser.parse_args()

g = nx.MultiGraph()
for nxp_path in args.network_pickles:
    # For each path of nxp files given
    condition_name = os.path.splitext(os.path.basename(nxp_path))[0]
    # figure out the name of the condition
    print("Processing ",condition_name)
    # report that
    ingraph = nx.read_gpickle(nxp_path)
    # read in the gpickle
    for s_node in ingraph.adj:
        # iterate through all source nodes
        for t_node,edge_datar in ingraph.adj[s_node].items():
            # grab all target nodes and edges
            try:
                # try incase you aren't yet defined
                for current_edge in g.adj[s_node][t_node]:
                    # for each edge
                    same = (
                        g.adj[s_node][t_node][current_edge]["condition"] == condition_name and
                        g.adj[s_node][t_node][current_edge]['fitness'] == edge_datar['fitness'] and 
                        g.adj[s_node][t_node][current_edge]['error'] == edge_datar['error'] and 
                        g.adj[s_node][t_node][current_edge]['fdr'] == edge_datar['fdr']
                        )
                    # are they the same edge?
                    if same:
                        break
                    # if so, break with that value
            except:
                same = False
                # otherwise, you're false
            if not same:
                # if you're false, then it's new
                g.add_edge(s_node,t_node,
                    fitness=edge_datar['fitness'],
                    error=edge_datar['error'],
                    fdr=edge_datar['fdr'],
                    condition=condition_name
                    )
                # add edge
            else:
                pass
                # otherwise, pass

nx.write_gpickle(g,args.output_base+".nxp")
nx.write_graphml(g,args.output_base+".graphml")
