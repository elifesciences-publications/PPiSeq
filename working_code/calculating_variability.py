#!/usr/bin/env python3

import argparse
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import re
import os
import csv

def csvprinter(obj,header,filename):
    """
    This prints obj to a file using the stats_base argument as a prefix then 
    the name_suffix to end the file. header is a tuple of column names.
    Prints a CSV to that file.
    """
    if filename:
        with open(filename,"w") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            if type(obj) == type(list()):
                writer.writerows(obj)
            if type(obj) == type(dict()):
                writer.writerows([ (i,obj[i]) for i in obj])

parser = argparse.ArgumentParser()
parser.add_argument("input_nxp",type=str,
    help="I expect a networkx network as a pickle, should be a MultiGraph "+
        "made from all the condition graphs, merged."
    )
parser.add_argument("--output_base",required=True,
    help="basename for output stats"
    )
#parser.add_argument("--stats_base",default=False,
#    help="basename for output stats files"
#    )
args = parser.parse_args()

multi_net = nx.read_gpickle(args.input_nxp)

# calculate variability for an interaction across environments
# calculate variability of interactions within environments
# sum these up in different ways, and spit them out

entropy_by_node = dict()
entropy_by_edge = dict()

conditions = set()
for e in multi_net.edges():
    [ conditions.add(val['condition']) for i,val in multi_net.get_edge_data(e[0],e[1]).items()]
print("I see conditions: ", conditions)

def normed_by_sum(inz):
    return(inz/np.sum(inz))

def normed_by_max(inz):
    return(inz/np.max(inz))

for u in multi_net.adj:
    u_array = pd.DataFrame(0,index=multi_net.adj[u].keys(),columns=conditions)
    for v in multi_net.adj[u]:
        for i,edge in multi_net.adj[u][v].items():
            u_array.loc[v,edge['condition']] = float(edge['fitness'])
    u_array_by_node = u_array.copy()
    u_array_by_node = u_array_by_node / u_array_by_node.values.sum()
    frac_active = u_array_by_node.apply(np.sum,axis=0)
    entropy_by_node[u] = frac_active.apply(lambda x: -x*np.log2(x) ).sum()

csvprinter(entropy_by_node,("YORF","entropy_across_conditions"),args.output_base+"_node_entropy.csv")

#    by_condition = {}
#    all_partners = []
#    for v in g.adj[u]:
#        all_partners.append(v)
#        datar = g.get_edge_data(u,v)
#        for id,each_edge in datar.items():
#            try:
#                by_condition[each_edge['condition']].append(float(each_edge['fitness']))
#            except:
#                by_condition[each_edge['condition']] = [float(each_edge['fitness'])]
#    hubiness[u] = len(all_partners) >= 3
#    max_interactions = 0
#    max_fitness_sum = 0
#    for condition in by_condition:
#        if len(by_condition[condition]) > max_interactions:
#            max_interactions = len(by_condition[condition])
#        if np.sum(by_condition[condition]) > max_fitness_sum:
#            max_fitness_sum = np.sum(by_condition[condition])
#    interaction_entropy = 0
#    for condition in by_condition:
#        frac = len(by_condition[condition]) / max_interactions
#        interaction_entropy += -frac*np.log2(frac)
#    interaction_entropies[u] = interaction_entropy
#    fitness_entropy = 0
#    for condition in by_condition:
#        frac = np.sum(by_condition[condition]) / max_fitness_sum
#        fitness_entropy += -frac*np.log2(frac)
#    fitness_entropies[u] = fitness_entropy
#    properties[u] = {'degree':len(all_partners),'hubiness':hubiness[u],'interaction_entropy':interaction_entropies[u],'fitness_entropy':fitness_entropies[u]}
#nx.set_node_attributes(g,properties)
#nx.write_graphml(g,"../tmp/annotated_multigraph.graphml")
#
#hist = plt.figure()
#plt.hist([i for k,i in interaction_entropies.items()])
##plt.semilogx()
#plt.title("Barplot of node degrees")
#plt.ylabel("Count of nodes")
#plt.xlabel("Degree of node")
#hist.show()
#hist.savefig("../tmp/interaction_entropies_hist.png")
#
#hist = plt.figure()
#plt.hist([i for k,i in fitness_entropies.items()])
##plt.semilogx()
#plt.title("Barplot of node degrees")
#plt.ylabel("Count of nodes")
#plt.xlabel("Degree of node")
#hist.show()
#hist.savefig("../tmp/interaction_entropies_hist.png")
#
#scatter = plt.figure()
#plt.scatter([i for k,i in nx.get_node_attributes(g,'interaction_entropy').items()],
#    [i for k,i in nx.get_node_attributes(g,'fitness_entropy').items()])
#plt.xlabel("interaction entropy")
#plt.ylabel("fitness entropy")
#plt.show()
#
#pos = nx.layout.spring_layout(g)
#
#np.max([g.nodes[k]['interaction_entropy'] for k in g.nodes])
#
#plt.figure(figsize=(16,8))
#nodes = nx.draw_networkx_nodes(g, pos,
#    node_size=np.log2(de)*2,
#    node_color=[g.nodes[k]['interaction_entropy'] for k in g.nodes])
#edges = nx.draw_networkx_edges(g, pos, width=0.1, alpha=0.2)
#nodes.set_norm(matplotlib.colors.NoNorm())
#plt.colorbar(nodes)
#plt.legend()
#axes = plt.gca()
#axes.set_xlim([-.2,.2])
#axes.set_ylim([-.2,.2])
#plt.show()
#
##comms = nx.algorithms.community.asyn_lpa_communities(g)
#de = [ i[1] for i in this_net.degree()]
#
##csvprinter(this_net.degree(),("YORF","degree"),"_node_degrees.csv")
#
#degree_hist = plt.figure()
#plt.hist(de)
#plt.title(net_name+": barplot of node degrees")
#plt.ylabel("Count of nodes")
#plt.xlabel("Degree of node")
#if args.figure_base:
#    degree_hist.savefig(args.figure_base+net_name+"_degree_histogram.svg")
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
#
#



