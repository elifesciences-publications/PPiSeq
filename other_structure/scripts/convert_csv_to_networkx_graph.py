#!/usr/bin/env python3

# This expects this schema:
# PPI,Number_of_Barcodes,Barcode_sequences,Fitness,Fitness_estimaion_error,Counts_G0,Counts_G6,Counts_G12,Counts_G18
# or 
# PPI,Number_of_barcodes,Mean_fitness,SD,P_value,FDR_adjusted_P_value,Positive

import argparse
import csv
import networkx
import pickle
import re

parser = argparse.ArgumentParser()
parser.add_argument("input_csv",type=str,
    help="gimme a csv like those that Zhimin makes named "+
        "'PPI_barcodes_fitness_counts'"
    )
parser.add_argument("output_file",type=str,
    help="output path"
    )
parser.add_argument("--type",type=str,
    help="raw or mean?"
    )
args = parser.parse_args()

g = networkx.MultiDiGraph()

if args.type == "raw":
    with open(args.input_csv,'r') as f:
        c = csv.reader(f,delimiter=',')
        next(c)
        for i in c:
            idz = i[0].split("_")
            g.add_edge(idz[0],idz[1],fitness=i[3],error=i[4])
elif args.type == "mean":
    with open(args.input_csv,'r') as f:
        c = csv.reader(f,delimiter=',')
        next(c)
        for i in c:
            if i[6] == 0:
                continue
            idz = i[0].split("_")
            g.add_edge(idz[0],idz[1],fitness=i[2],error=i[3],fdr=i[5])
else:
    raise()

networkx.write_gpickle(g,args.output_file)
