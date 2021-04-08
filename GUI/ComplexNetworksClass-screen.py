#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import numpy as np
import warnings
import pandas as pd
from gooey import Gooey, GooeyParser
from igraph import *
from Bio import SeqIO
warnings.filterwarnings("ignore")
path = os.path.dirname(__file__)
if path != '':
    path = path + '/'

 
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
    metrics.append(thresholdCN.assortativity_degree(directed=False))  # Returns the assortativity
    metrics.append(max(thresholdCN.degree()))
    metrics.append(min(thresholdCN.degree()))
    metrics.append(np.std(thresholdCN.degree()))  # Returns the strength (weighted degree)
    metrics.append(thresholdCN.average_path_length(directed=False, unconn=False))  # Average path length
    metrics.append(thresholdCN.transitivity_avglocal_undirected())  # local transitivity (clustering coefficient)
    metrics.append(thresholdCN.transitivity_undirected())  # global transitivity (clustering coefficient)
    metrics.append(mean(cn.ecount()))  # Counts the number of edges
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
    file.write("%s," % name_seq)
    metrics_preprocessing = np.nan_to_num(metrics)
    for x in metrics_preprocessing:
        # print (x)
        file.write("%s," % x)
    file.write(label_dataset)
    file.write("\n")
    print("Recorded Sequence!!!")
    return
    
    
def complex_network(finput, foutput, label, ksize, threshold):
    """Generates complex network"""
    global name_seq, cn, metrics, thresholdCN
    header(foutput)  # Header
    for seq_record in SeqIO.parse(finput, "fasta"):  # Reading Sequences
        perm = []
        metrics = []
        cn = Graph()
        # print(summary(cn))
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        kmer = []
        for codon in patterns(seq, ksize):  # Generates every chosen k pattern
            kmer.append(str(codon))
        # print(kmer)
        for i in range(len(kmer)-1):  # Position -1 -- Build the Network
            cn.add_vertices(kmer[i]) if kmer[i] not in perm else kmer[i]
            cn.add_vertices(kmer[i+1]) if kmer[i+1] not in perm else kmer[i+1]
            cn.add_edges([(kmer[i], kmer[i+1])])
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
                    for i in range(1, 13):
                        metrics.append(0)
                break
            else:
                feature_extraction()
        file_record(foutput)
    # return cn.write_adjacency("cn")
    return


@Gooey(program_name='MathFeature',
       default_size=(800, 600),
       language='english',
       tabbed_groups=True,
       image_dir=path + 'img',
       menu=[{
           'name': 'File',
           'items': [{
               'type': 'AboutDialog',
               'menuTitle': 'About',
               'name': 'MathFeature',
               'description': 'Feature Extraction Package for Biological Sequences Based on Mathematical Descriptors',
               'version': '1.2.1',
               'copyright': '2021',
               'website': 'https://bonidia.github.io/MathFeature/',
               'developer': 'https://bonidia.github.io/website/'},
               {
                   'type': 'Link',
                   'menuTitle': 'Visit Our Site',
                   'url': 'https://bonidia.github.io/MathFeature/'
               }]},
           {
               'name': 'Help',
               'items': [{
                   'type': 'Link',
                   'menuTitle': 'Documentation',
                   'url': 'https://bonidia.github.io/MathFeature/'
               }]
           }])
def main():
    parser = GooeyParser(description="Feature Extraction Package for Biological Sequences")

    screen = parser.add_argument_group('Complex Network')
    # test = parser.add_argument_group('ORF Descriptor')
    screen.add_argument('-i',
                        '--input_file',
                        metavar='Input File',
                        help='Fasta format file, e.g., test.fasta',
                        required=True,
                        widget='FileChooser')
    screen.add_argument('-o',
                        '--output_file',
                        metavar='Output File',
                        help='CSV format file, e.g., test.csv',
                        required=True,
                        widget='FileSaver')
    screen.add_argument('-l',
                        '--label',
                        metavar='Dataset Label',
                        help='Dataset Label, e.g., 0, 1, lncRNA, mRNA, sncRNA.',
                        required=True,
                        action='store')
    screen.add_argument('-k',
                        '--pattern',
                        metavar='Frequency of k-mer',
                        help='Frequency of k-mer, e.g., 2, 3 (default = codon), 4...',
                        default=3,
                        required=True,
                        action='store')
    screen.add_argument('-t',
                        '--threshold',
                        metavar='Threshold Size',
                        help='threshold size, e.g., 2, 3, 20...',
                        required=True,
                        action='store')
    return parser.parse_args()


#############################################################################    
if __name__ == "__main__":
    args = main()
    ffinput = str(args.input_file)
    ffoutput = str(args.output_file)
    label_dataset = str(args.label)
    ksize = int(args.pattern)
    threshold = int(args.threshold) + 1  # always +1
    start_time = time.time()
    complex_network(ffinput, ffoutput, label_dataset, ksize, threshold)
    print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
