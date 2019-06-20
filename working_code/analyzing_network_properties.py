#!/usr/bin/env python3

import argparse
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as axs
import matplotlib
import os
import re
import csv

def csvprinter(obj,header,name_suffix):
    """
    This prints obj to a file using the stats_base argument as a prefix then 
    the name_suffix to end the file. header is a tuple of column names.
    Prints a CSV to that file.
    """
    if args.stats_base:
        with open(args.stats_base+name_suffix,"w") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(obj)

parser = argparse.ArgumentParser()
parser.add_argument("input_nxp",type=str,
    help="I expect a networkx network as a pickle, so nxp"
    )
parser.add_argument("--figure_base",default=False,
    help="basename for output figures"
    )
parser.add_argument("--stats_base",default=False,
    help="basename for output stats files"
    )
args = parser.parse_args()


this_net = nx.read_gpickle(args.input_nxp)
net_name = re.sub(".nxp","",os.path.basename(args.input_nxp))

print("Proc'ing",net_name)

de = [ i[1] for i in this_net.degree()]

csvprinter(this_net.degree(),("YORF","degree"),"_node_degrees.csv")

degree_hist = plt.figure()
plt.hist(de)
plt.title(net_name+": barplot of node degrees")
plt.ylabel("Count of nodes")
plt.xlabel("Degree of node")
if args.figure_base:
    degree_hist.savefig(args.figure_base+net_name+"_degree_histogram.svg")

hairball = plt.figure(figsize=(32, 32))
pos = nx.layout.spring_layout(this_net)
node_sizes = np.sqrt(de)
nodes = nx.draw_networkx_nodes(this_net, pos, node_size=node_sizes, node_color='blue')
try:
    line_widths = [ float(this_net.get_edge_data(u,v)['fitness']) for u,v in this_net.edges() ]
    edges = nx.draw_networkx_edges(this_net, pos, node_size=node_sizes, width=line_widths, alpha=0.5)
except:
    edges = nx.draw_networkx_edges(this_net, pos, node_size=node_sizes, width=0.2, alpha=0.5)
plt.title(net_name+": net of all connections")
if args.figure_base:
    hairball.savefig(args.figure_base+net_name+"_hairball.svg")

if type(this_net) == type(nx.Graph()):
    ecorig = nx.eigenvector_centrality(this_net)
    ec = list(ecorig.values())
    
    csvprinter(ecorig,("YORF","eigenvector_centrality"),"_node_eigenvectorcentrality.csv")
    
    pos = nx.layout.spring_layout(this_net)
    eigenvector_centrality_hist = plt.figure()
    plt.hist(ec)
    #plt.semilogx()
    plt.title("Barplot of node eigenvector centrality")
    plt.ylabel("Count of nodes")
    plt.xlabel("Eigenvector centrality")
    if args.figure_base:
        eigenvector_centrality_hist.savefig(args.figure_base+net_name+"_eigencenter_histogram.svg")
    
    
    pos = nx.layout.spring_layout(this_net)
    eigencenter_hairball = plt.figure(figsize=(16,8))
    nodes = nx.draw_networkx_nodes(this_net, pos, node_size=5, node_color=ec)
    edges = nx.draw_networkx_edges(this_net, pos, width=0.1, alpha=0.2)
    nodes.set_norm(matplotlib.colors.Normalize())
    plt.colorbar(nodes)
    plt.legend()
    axes = plt.gca()
    axes.set_xlim([-.2,.2])
    axes.set_ylim([-.2,.2])
    if args.figure_base:
        eigencenter_hairball.savefig(args.figure_base+net_name+"_eigencenter_hairball.svg")
