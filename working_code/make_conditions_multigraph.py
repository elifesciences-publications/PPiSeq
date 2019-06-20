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
    condition_name = os.path.splitext(os.path.basename(nxp_path))[0]
    ingraph = nx.read_gpickle(nxp_path)
    for s_node in ingraph.adj:
        for t_node,edge_datar in ingraph.adj[s_node].items():
            g.add_edge(s_node,t_node,
                fitness=edge_datar['fitness'],
                error=edge_datar['error'],
                fdr=edge_datar['fdr'],
                condition=condition_name
                )

nx.write_gpickle(g,args.output_base+".nxp")
nx.write_graphml(g,args.output_base+".graphml")
