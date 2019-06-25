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
import pandas as pd

def process_graph(this_net):

    degrees = pd.Series(dict(this_net.degree()))
    clustering_triangles = pd.Series(nx.algorithms.cluster.clustering(this_net))
    clustering_squares = pd.Series(nx.algorithms.cluster.square_clustering(this_net))
    eigen_centrality = pd.Series(nx.algorithms.centrality.eigenvector_centrality(this_net))
    betweenness_centrality = pd.Series(nx.algorithms.centrality.betweenness_centrality(this_net))
    together = pd.DataFrame(
        dict(
            degree=degrees,
            clustering_coef=clustering_triangles,
            clustering_coef_squares=clustering_squares,
            eigen_centrality=eigen_centrality,
            betweenness_centrality=betweenness_centrality
            )
        )
    
    if args.stats_base:
        together.to_json(args.stats_base+net_name+"_topological_properties.json",orient='index')
    
    this_fig = plt.figure()
    together.plot(kind="hist",y="degree",bins=40)
    plt.yscale("log")
    plt.title(net_name+": barplot of node degrees")
    plt.ylabel("Count of nodes")
    plt.xlabel("Degree of node")
    if args.figure_base:
        plt.savefig(args.figure_base+net_name+"_degree_histogram.svg")
    
    this_fig = plt.figure()
    together.plot(kind="hist",y="eigen_centrality",bins=40)
    plt.yscale("log")
    plt.title(net_name+": barplot of ")
    plt.ylabel("Count of nodes")
    plt.xlabel(" of node")
    if args.figure_base:
        plt.savefig(args.figure_base+net_name+"_eigenvector.svg")

    this_fig = plt.figure()
    together.plot(kind="hist",y="betweenness_centrality",bins=40)
    plt.yscale("log")
    plt.title(net_name+": barplot of ")
    plt.ylabel("Count of nodes")
    plt.xlabel(" of node")
    if args.figure_base:
        plt.savefig(args.figure_base+net_name+"_betweenness.svg")

    this_fig = plt.figure()
    together.plot(kind="scatter",x="degree",y="clustering_coef")
    plt.title(net_name+": scatter")
    plt.xlabel("Degree of node")
    plt.ylabel("Clustering coefficient, triangles")
    if args.figure_base:
        plt.savefig(args.figure_base+net_name+"_degree_vs_clustering_triangles.svg")

    this_fig = plt.figure()
    together.plot(kind="scatter",x="degree",y="clustering_coef_squares")
    plt.title(net_name+": scatter")
    plt.xlabel("Degree of node")
    plt.ylabel("Clustering coefficient, squares")
    if args.figure_base:
        plt.savefig(args.figure_base+net_name+"_degree_vs_clustering_squares.svg")

    #
    #
    #
    
    label_communities = nx.algorithms.community.label_propagation.label_propagation_communities(this_net)
    modularity_communities = nx.algorithms.community.modularity_max.greedy_modularity_communities(this_net)


#
#
#

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

# this_net = nx.read_gpickle("/home/zed/Dropbox/PPiSeq_02/Working_data/networks/NaCl_mean_fitness_positive.nxp")

print("Proc'ing",net_name)

if type(this_net) is type(nx.Graph()):
    process_graph(this_net)




#
#
#
#
#hairball = plt.figure(figsize=(32, 32))
#pos = nx.layout.spring_layout(this_net)
#node_sizes = np.sqrt(de)
#nodes = nx.draw_networkx_nodes(this_net, pos, node_size=node_sizes, node_color='blue')
#try:
#    line_widths = [ float(this_net.get_edge_data(u,v)['fitness']) for u,v in this_net.edges() ]
#    edges = nx.draw_networkx_edges(this_net, pos, node_size=node_sizes, width=line_widths, alpha=0.5)
#except:
#    edges = nx.draw_networkx_edges(this_net, pos, node_size=node_sizes, width=0.2, alpha=0.5)
#plt.title(net_name+": net of all connections")
#if args.figure_base:
#    hairball.savefig(args.figure_base+net_name+"_hairball.svg")
#
##
##
##
#
#if type(this_net) == type(nx.Graph()):
#    ecorig = nx.eigenvector_centrality(this_net)
#    ec = list(ecorig.values())
#    
#    csvprinter(ecorig,("YORF","eigenvector_centrality"),"_node_eigenvectorcentrality.csv")
#    
#    pos = nx.layout.spring_layout(this_net)
#    eigenvector_centrality_hist = plt.figure()
#    plt.hist(ec)
#    #plt.semilogx()
#    plt.title("Barplot of node eigenvector centrality")
#    plt.ylabel("Count of nodes")
#    plt.xlabel("Eigenvector centrality")
#    if args.figure_base:
#        eigenvector_centrality_hist.savefig(args.figure_base+net_name+"_eigencenter_histogram.svg")
#    
#    
#    pos = nx.layout.spring_layout(this_net)
#    eigencenter_hairball = plt.figure(figsize=(16,8))
#    nodes = nx.draw_networkx_nodes(this_net, pos, node_size=5, node_color=ec)
#    edges = nx.draw_networkx_edges(this_net, pos, width=0.1, alpha=0.2)
#    nodes.set_norm(matplotlib.colors.Normalize())
#    plt.colorbar(nodes)
#    plt.legend()
#    axes = plt.gca()
#    axes.set_xlim([-.2,.2])
#    axes.set_ylim([-.2,.2])
#    if args.figure_base:
#        eigencenter_hairball.savefig(args.figure_base+net_name+"_eigencenter_hairball.svg")
