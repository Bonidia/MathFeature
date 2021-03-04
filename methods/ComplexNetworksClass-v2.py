#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import warnings
import argparse
from igraph import *
from Bio import SeqIO
warnings.filterwarnings('ignore')

 
def header(file):
    file = open(file, 'a')
    file.write('nameseq,')
    for i in range(1, ksize + 1):
        file.write('avg_betweenness.{0}-mer,avg_degree.{0}-mer,assortativity_degree.{0}-mer,'
                   'maxdegree.{0}-mer,mindegree.{0}-mer,std_degree.{0}-mer,average_path_length.{0}-mer,'
                   'transitivity_avglocal_undirected.{0}-mer,transitivity_undirected.{0}-mer,'
                   'number_of_edges.{0}-mer,motifs3.{0}-mer,motifs4.{0}-mer,authority_score.{0}-mer,'
                   'closeness.{0}-mer,Burts_constraint.{0}-mer,count_multiple.{0}-mer,density.{0}-mer,'
                   'diameter.{0}-mer,eccentricity.{0}-mer,edge_betweenness.{0}-mer,'
                   'hub_score.{0}-mer,maxdegreetwo.{0}-mer,neighborhood_size.{0}-mer,radius.{0}-mer,'
                   'avg_weighted_degree.{0}-mer,number_vertices.{0}-mer,'.format(i))
    file.write('label')
    file.write('\n')
    

def feature_extraction(thresholdCN):
    metrics.append(mean(thresholdCN.betweenness(directed=False, weights=None, nobigint=True)))
    metrics.append(mean(thresholdCN.degree()))
    metrics.append(thresholdCN.assortativity_degree(directed=False))  # Returns the assortativity
    metrics.append(max(thresholdCN.degree()))
    metrics.append(min(thresholdCN.degree()))
    metrics.append(np.std(thresholdCN.degree()))  # Returns the strength (weighted degree)
    metrics.append(thresholdCN.average_path_length(directed=False, unconn=False))  # Average path length
    metrics.append(thresholdCN.transitivity_avglocal_undirected())  # local transitivity (clustering coefficient)
    metrics.append(thresholdCN.transitivity_undirected())  # global transitivity (clustering coefficient)
    metrics.append(cn.ecount())  # Counts the number of edges
    metrics.append(thresholdCN.motifs_randesu_no(size=3))
    metrics.append(thresholdCN.motifs_randesu_no(size=4))
    metrics.append(mean(thresholdCN.authority_score()))
    metrics.append(mean(thresholdCN.closeness(vertices=None, mode=ALL, cutoff=None, weights=None, normalized=True)))  # Calculates the closeness centralities of given vertices in a graph
    metrics.append(mean(thresholdCN.constraint(vertices=None, weights=None)))  # Calculates Burt's constraint scores for given vertices in a graph.
    metrics.append(mean(thresholdCN.count_multiple(edges=None)))  # Counts the multiplicities of the given edges.
    metrics.append(thresholdCN.density(loops=False))  # Calculates the density of the graph.
    metrics.append(thresholdCN.diameter(directed=False, unconn=False, weights=None))  # Calculates the diameter of the graph.
    metrics.append(mean(thresholdCN.eccentricity(vertices=None, mode=ALL)))  # Calculates the eccentricities of given vertices in a graph.
    metrics.append(mean(thresholdCN.edge_betweenness(directed=False, cutoff=None, weights=None)))  # Calculates or estimates the edge betweennesses in a graph.
    metrics.append(mean(thresholdCN.hub_score()))  # Calculates Kleinberg's hub score for the vertices of the graph.
    metrics.append(thresholdCN.maxdegree())  # Returns the maximum degree of a vertex set in the graph.
    metrics.append(mean(thresholdCN.neighborhood_size()))  # For each vertex specified by vertices, returns the number of vertices reachable from that vertex in at most order steps
    metrics.append(thresholdCN.radius())  # Calculates the radius of the graph.
    metrics.append(mean(thresholdCN.strength()))  # Returns the strength (weighted degree) of some vertices from the graph.
    metrics.append(cn.vcount())  # Counts the number of vertices.
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

            
def file_record(name_seq, foutput, label):
    """Generates CSV with the extracted features"""
    file = open(foutput, 'a')
    file.write('%s,' % name_seq)
    metrics_preprocessing = np.nan_to_num(metrics)
    for x in metrics_preprocessing:
        # print (x)
        file.write('%s,' % x)
    file.write(label)
    file.write('\n')
    print('Recorded Sequence!!!')
    return
    
    
def complex_network(finput, foutput, label, ksize):
    """Generates complex network"""
    global name_seq, cn, metrics
    header(foutput)  # Header
    for seq_record in SeqIO.parse(finput, 'fasta'):  # Reading Sequences
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        metrics = []
        for k in range(1, ksize + 1):
            cn = Graph()
            kmer = []
            for subseq in patterns(seq, k):  # Generates k pattern
                kmer.append(str(subseq))
            # print(kmer)
            vertices = np.unique(kmer)
            for vert in vertices:
                cn.add_vertices(vert)
            for i in range(len(kmer)-1):  # Position -1 -- Build the Network
                cn.add_edges([(kmer[i], kmer[i+1])])
            # print(summary(cn))
            feature_extraction(cn)
        file_record(name_seq, foutput, label)
    # return cn.write_adjacency("cn")
    return


#############################################################################    
if __name__ == '__main__':
    print('\n')
    print('###################################################################################')
    print('###################### Feature Extraction: Complex Network ########################')
    print('############    Arguments: -i input -o output -l label -k pattern    ##############')
    print('##########               Author: Robson Parmezan Bonidia                ###########')
    print('###################################################################################')
    print('\n')
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Fasta format file, E.g., test.fasta')
    parser.add_argument('-o', '--output', help='CSV format file, E.g., test.csv')
    parser.add_argument('-l', '--label', help='Dataset Label, E.g., lncRNA, mRNA, sncRNA ...')
    parser.add_argument('-k', '--pattern', help='k size (k-mer), E.g., 2, 3, 4 ...')
    args = parser.parse_args()
    finput = str(args.input)
    foutput = str(args.output)
    label_dataset = str(args.label)
    ksize = int(args.pattern)
    complex_network(finput, foutput, label_dataset, ksize)
#############################################################################
