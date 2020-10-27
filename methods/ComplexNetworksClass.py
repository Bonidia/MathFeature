#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import random
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import argparse
from igraph import *
from Bio import SeqIO

 
def header(file):
    """
    AB = avg_betweenness
    AEB = avg_edge_betweenness
    AD = avg_degree
    ASSD = assortativity_degree
    MAXD = maxdegree
    MIND = mindegree
    AS = avg_strength
    APL = average_path_length
    ASPL =average_shortest_paths
    TALU = transitivity_avglocal_undirected
    TU = transitivity_undirected
    AEC = avg_eigenvector_centrality
    NE = number_of_edges
    MOT3 = motifs_randesu_no_3 
    MOT4 = motifs_randesu_no_4
    """
    file = open(file, 'a')
    file.write("nameseq,")
    for i in range(1, threshold):
	    file.write("AB.{0},AD.{0},ASSD.{0},MAXD.{0},MIND.{0},AS.{0},APL.{0},TALU.{0},TU.{0},NE.{0},MOT3.{0},MOT4.{0},".format(i))
    file.write("label")
    file.write("\n")
    

def feature_extraction():
    metrics.append(mean(thresholdCN.betweenness(directed=False, weights=None, nobigint=True)))
    metrics.append(mean(thresholdCN.degree()))
    metrics.append(thresholdCN.assortativity_degree(directed=False)) # Returns the assortativity
    metrics.append(max(thresholdCN.degree()))
    metrics.append(min(thresholdCN.degree()))
    metrics.append(np.std(thresholdCN.degree())) # Returns the strength (weighted degree)
    metrics.append(thresholdCN.average_path_length(directed=False, unconn=False)) # Average path length
    metrics.append(thresholdCN.transitivity_avglocal_undirected()) # local transitivity (clustering coefficient) 
    metrics.append(thresholdCN.transitivity_undirected()) # global transitivity (clustering coefficient) 
    metrics.append(mean(cn.ecount())) # Counts the number of edges
    metrics.append(thresholdCN.motifs_randesu_no(size=3))
    metrics.append(thresholdCN.motifs_randesu_no(size=4))
    return


def patterns(seq, win):
    """
    Generate k-mers: subsequences of length k 
    contained in a biological sequence.
    """
    seqlen = len(seq)
    for i in range(seqlen):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return

            
def file_record(foutput):
    """Generates CSV with the extracted features"""
    file = open(foutput, 'a')
    file.write("%s," % (name_seq))
    metrics_preprocessing = np.nan_to_num(metrics)
    for x in metrics_preprocessing:
        # print (x)
        file.write("%s," % (x))
    file.write(label_dataset)
    file.write("\n")
    print ("Recorded Sequence!!!")
    return
    
    
def complex_network(finput,foutput,label,ksize,threshold):
    """Generates complex network"""
    global name_seq, cn, metrics, thresholdCN
    header(foutput) # Header
    for seq_record in SeqIO.parse(finput, "fasta"): # Reading Sequences
        perm = []
        metrics = []
        cn = Graph()
        # print(summary(cn))
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        kmer = []
        for codon in patterns(seq, ksize): # Generates every chosen k pattern
            kmer.append(str(codon))
        # print(kmer)
        for i in range(len(kmer)-1): # Position -1 -- Build the Network
            cn.add_vertices(kmer[i]) if kmer[i] not in perm else kmer[i]
            cn.add_vertices(kmer[i+1]) if kmer[i+1] not in perm else kmer[i+1]
            cn.add_edges([(kmer[i],kmer[i+1])])
            perm.append(kmer[i]) if kmer[i] not in perm else kmer[i]
            perm.append(kmer[i+1]) if kmer[i+1] not in perm else kmer[i+1]
        # print(summary(cn))
        for t in range(1, threshold):
            """
            Extract features using a threshold scheme.
            Eliminates each edge of size < t and extract features again.
            """
            matrix_adj = pd.DataFrame(cn.get_adjacency())
            thresholdCN = np.where(matrix_adj < t, 0, matrix_adj)
            thresholdCN = Graph.Adjacency(thresholdCN.astype(int).tolist(), mode=ADJ_UNDIRECTED)
            # print(t)
            if thresholdCN.ecount() < 1:
                for i in range(t, threshold):
                    for i in range(1,13):
                        metrics.append(0)
                break
            else:
                feature_extraction()
        file_record(foutput)
    # return cn.write_adjacency("cn")
    return


#############################################################################    
if __name__ == "__main__":
    print("\n")
    print("###################################################################################")
    print("###################### Feature Extraction: Complex Network ########################")
    print("######    Arguments: -i input -o output -l label -k pattern -t threshold   ########")
    print("##########               Author: Robson Parmezan Bonidia                ###########")
    print("###################################################################################")
    print("\n")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Fasta format file, E.g., test.fasta')
    parser.add_argument('-o', '--output', help='CSV format file, E.g., test.csv')
    parser.add_argument('-l', '--label', help='Dataset Label, E.g., lncRNA, mRNA, sncRNA ...')
    parser.add_argument('-k', '--pattern', help='k size, E.g., 2, 3 (default = codon), 4 ...')
    parser.add_argument('-t', '--threshold', help='threshold size, E.g., 2, 3 (default = 20) ...')
    args = parser.parse_args()
    finput = str(args.input)
    foutput = str(args.output)
    label_dataset = str(args.label)
    ksize = int(args.pattern)
    threshold = int(args.threshold) + 1 # always +1
    complex_network(finput,foutput,label_dataset,ksize,threshold)
#############################################################################