#!/usr/bin/env python
# -*- coding: utf-8 -*-


#############################################################################

import sys
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import collections
import argparse
from Bio import SeqIO
from itertools import product


#############################################################################


def headerNAC(k):
    global listTotal, caracteres
    listTotal = []
    file = open(foutput, 'a')
    file.write('%s,' % ('nameseq'))
    permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
    # print (permsList)
    for perm in permsList:
        # print (perm)
        listTotal.append(perm)
        file.write('%s,' % (str(perm)))
    file.write('label')
    file.write('\n')
    return


def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0, seqlen, step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return        
    

def chunks_two(seq, win):
    seqlen = len(seq)
    for i in range(seqlen):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return

            
def file_record(name_seq):
    file = open(foutput, 'a')
    file.write('%s,' % (name_seq))
    for x in probabilities:
        # print (x)
        file.write('%s,' % (str(x[1])))
    file.write(labelDataset)
    file.write('\n')
    print ('Recorded Sequence: %s' % (name_seq))
    return
    

# Amino Acid Composition (AAC), Dipeptide composition (DPC), Tripeptide composition (TPC)
def nacSeq(k):
    global probabilities
    headerNAC(k)
    for seq_record in SeqIO.parse(finput, 'fasta'):
        seq = seq_record.seq
        seq = seq.upper()	
        name_seq = seq_record.name
        probabilities = []
        kmer = {}
        totalWindows = (len(seq) - k) + 1 # (L - k + 1)
        permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
        for key in permsList:
            kmer[key] = 0
        kmer = collections.OrderedDict(sorted(kmer.items()))
        for subseq in chunks_two(seq, k):
            if subseq in kmer:
                # print(subseq)
                kmer[subseq] = kmer[subseq] + 1
            else:
                kmer[subseq] = 1
        for key, value in kmer.items():
            # print (key)
            # print (value)
            probabilities.append([str(key), value/totalWindows])
        file_record(name_seq)
    return


#############################################################################


def header_kmer(ksize):
    global listTotal,caracteres
    listTotal = []
    file = open(foutput,'a')
    file.write('%s,' % ('nameseq'))
    for k in range(1, ksize + 1):
        permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
        # print (permsList)
        for perm in permsList:
            # print (perm)
            listTotal.append(perm)
            file.write('%s,' % (str(perm)))
    file.write('label')
    file.write('\n')
    return


# k-mer
def findKmers(ksize):
    global probabilities
    header_kmer(ksize)
    for seq_record in SeqIO.parse(finput,'fasta'):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        probabilities = []
        for k in range(1, ksize + 1):
            kmer = {}
            totalWindows = (len(seq) - k) + 1  # (L - k + 1)
            permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
            for key in permsList:
                kmer[key] = 0
            kmer = collections.OrderedDict(sorted(kmer.items()))
            for subseq in chunks_two(seq,k):
                if subseq in kmer:
                    # print(subseq)
                    kmer[subseq] = kmer[subseq] + 1
                else:
                    kmer[subseq] = 1
            for key, value in kmer.items():
                # print (key)
                # print (value)
                probabilities.append([str(key), value/totalWindows])
        file_record(name_seq)
    return


#############################################################################


def chunks_kgap(seq, win, step):
    seqlen = len(seq)
    for i in range(0, seqlen, step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        j += (step - 1)
        if j>=seqlen: break
    return


def kgap_record(name_seq):
    file = open(foutput, 'a')
    file.write('%s,' % (name_seq))
    for x in number_kmers:
        # print (x)
        file.write('%s,' % (str(x[1])))
    file.write(labelDataset)
    file.write('\n')
    print ('Recorded Sequence: %s' % (name_seq))
    return


# xxKGAP - Correct with step - new approach - my idea
def seqKGAP(window, ksize, step):
    global number_kmers
    header_kmer(ksize)
    for seq_record in SeqIO.parse(finput,'fasta'):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        number_kmers = []
        for k in range(1, ksize + 1):
            kmer = {}
            permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
            for key in permsList:
                kmer[key] = 0
            kmer = collections.OrderedDict(sorted(kmer.items()))
            for subseq in chunks_kgap(seq, window, step):
                for subsubseq in chunks_two(subseq, k):
                    if subsubseq in kmer:
                        kmer[subsubseq] = kmer[subsubseq] + 1
                    else:
                        kmer[subsubseq] = 1
            totalWindows = sum(kmer.values())
            for key, value in kmer.items():
                # print (key)
                # print (value)
                number_kmers.append([str(key), value/totalWindows])
        kgap_record(name_seq)
    return

        
#############################################################################    
if __name__ == '__main__':
    print('\n')
    print('###################################################################################')
    print('##########      Feature Extraction: Several Techniques - Protein        ###########')
    print('##########               Author: Robson Parmezan Bonidia                ###########')
    print('###################################################################################')
    print('\n')
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Fasta format file, E.g., test.fasta')
    parser.add_argument('-o', '--output', help='CSV format file, E.g., test.csv')
    parser.add_argument('-l', '--label', help='Dataset Label, E.g., lncRNA, mRNA, sncRNA ...')
    parser.add_argument('-t', '--type', help='type of Feature Extraction...')
    args = parser.parse_args()
    finput = str(args.input)
    foutput = str(args.output)
    labelDataset = str(args.label)
    tech = str(args.type)
    #############################################################################
    # Amino | Code | Molar mass | Van der Waals Volume | V | Polarity(0 = Non, 1 = Polar) |
    # Acidity(1 = Basic (weak), 2 = Basic, 3 = Basic (strong), 4 = Neutral, 5 = Acidic ) | Hydropathy | Isoelectric.
    table = {
        'GCA': ["Alanine","Ala/A",89.09,67,92,0,4,1.8,6.01],
        'GCC': ["Alanine","Ala/A",89.09,67,92,0,4,1.8,6.01],
        'GCG': ["Alanine","Ala/A",89.09,67,92,0,4,1.8,6.01],
        'GCT': ["Alanine","Ala/A",89.09,67,92,0,4,1.8,6.01],
        'CGA': ["Arginine","Arg/R",174.2,148,225,1,3,-4.5,10.76],
        'CGC': ["Arginine","Arg/R",174.2,148,225,1,3,-4.5,10.76],
        'CGG': ["Arginine","Arg/R",174.2,148,225,1,3,-4.5,10.76],
        'CGT': ["Arginine","Arg/R",174.2,148,225,1,3,-4.5,10.76],
        'AGA': ["Arginine","Arg/R",174.2,148,225,1,3,-4.5,10.76],
        'AGG': ["Arginine","Arg/R",174.2,148,225,1,3,-4.5,10.76],
        'AAC': ["Asparagine","Asn/N",132.11,96,135,1,4,-3.5,5.41],
        'AAT': ["Asparagine","Asn/N",132.11,96,135,1,4,-3.5,5.41],
        'GAC': ["Aspartic Acid","Asp/D",133.1,91,125,1,5,-3.5,2.85],
        'GAT': ["Aspartic Acid","Asp/D",133.1,91,125,1,5,-3.5,2.85],
        'TGC': ["Cysteine","Cys/C",121.15,86,106,1,4,2.5,5.05],
        'TGT': ["Cysteine","Cys/C",121.15,86,106,1,4,2.5,5.05],
        'GAA': ["Glutamic Acid","Glu/E",147.13,109,161,1,5,-3.5,3.15],
        'GAG': ["Glutamic Acid","Glu/E",147.13,109,161,1,5,-3.5,3.15],
        'CAA': ["Glutamine","Gln/Q",146.15,114,155,1,4,-3.5,5.65],
        'CAG': ["Glutamine","Gln/Q",146.15,114,155,1,4,-3.5,5.65],
        'GGA': ["Glycine","Gly/G",75.06,48,66,0,4,-0.4,6.06],
        'GGC': ["Glycine","Gly/G",75.06,48,66,0,4,-0.4,6.06],
        'GGG': ["Glycine","Gly/G",75.06,48,66,0,4,-0.4,6.06],
        'GGT': ["Glycine","Gly/G",75.06,48,66,0,4,-0.4,6.06],
        'CAC': ["Histidine","His/H",155.15,118,167,1,1,-3.2,7.6],
        'CAT': ["Histidine","His/H",155.15,118,167,1,1,-3.2,7.6],
        'ATA': ["Isoleucine","Ile/I",131.17,124,169,0,4,4.5,6.05],
        'ATC': ["Isoleucine","Ile/I",131.17,124,169,0,4,4.5,6.05],
        'ATT': ["Isoleucine","Ile/I",131.17,124,169,0,4,4.5,6.05],
        'CTA': ["Leucine","Leu/L",131.17,124,168,0,4,3.8,6.01],
        'CTC': ["Leucine","Leu/L",131.17,124,168,0,4,3.8,6.01],
        'CTG': ["Leucine","Leu/L",131.17,124,168,0,4,3.8,6.01],
        'CTT': ["Leucine","Leu/L",131.17,124,168,0,4,3.8,6.01],
        'TTA': ["Leucine","Leu/L",131.17,124,168,0,4,3.8,6.01],
        'TTG': ["Leucine","Leu/L",131.17,124,168,0,4,3.8,6.01],
        'AAA': ["Lysine","Lys/K",146.18,135,171,1,2,-3.9,9.6],
        'AAG': ["Lysine","Lys/K",146.18,135,171,1,2,-3.9,9.6],
        'ATG': ["Methionine","Met/M",149.2,124,171,0,4,1.9,5.74],
        'TTC': ["Phenylalanine","Phe/F",165.19,135,203,0,4,2.8,5.49],
        'TTT': ["Phenylalanine","Phe/F",165.19,135,203,0,4,2.8,5.49],
        'CCA': ["Proline","Pro/P",115.13,90,129,0,4,-1.6,6.3],
        'CCC': ["Proline","Pro/P",115.13,90,129,0,4,-1.6,6.3],
        'CCG': ["Proline","Pro/P",115.13,90,129,0,4,-1.6,6.3],
        'CCT': ["Proline","Pro/P",115.13,90,129,0,4,-1.6,6.3],
        'TCA': ["Serine","Ser/S",105.09,73,99,1,4,-0.8,5.68],
        'TCC': ["Serine","Ser/S",105.09,73,99,1,4,-0.8,5.68],
        'TCG': ["Serine","Ser/S",105.09,73,99,1,4,-0.8,5.68],
        'TCT': ["Serine","Ser/S",105.09,73,99,1,4,-0.8,5.68],
        'AGC': ["Serine","Ser/S",105.09,73,99,1,4,-0.8,5.68],
        'AGT': ["Serine","Ser/S",105.09,73,99,1,4,-0.8,5.68],
        'ACA': ["Threonine","Thr/T",119.12,93,122,1,4,-0.7,5.6],
        'ACC': ["Threonine","Thr/T",119.12,93,122,1,4,-0.7,5.6],
        'ACG': ["Threonine","Thr/T",119.12,93,122,1,4,-0.7,5.6],
        'ACT': ["Threonine","Thr/T",119.12,93,122,1,4,-0.7,5.6],
        'TGG': ["Tryptophan","Trp/W",204.22,163,240,0,4,-0.9,5.89],
        'TAC': ["Tyrosine","Tyr/Y",181.19,141,203,1,4,-1.3,5.64],
        'TAT': ["Tyrosine","Tyr/Y",181.19,141,203,1,4,-1.3,5.64],
        'GTC': ["Valine","Val/V",117.14,105,142,0,4,4.2,6],
        'GTG': ["Valine","Val/V",117.14,105,142,0,4,4.2,6],
        'GTT': ["Valine","Val/V",117.14,105,142,0,4,4.2,6],
        'GTA': ["Valine","Val/V",117.14,105,142,0,4,4.2,6],
    }
    #############################################################################
    caracteres = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    #############################################################################
    if tech == 'AAC' or tech == 'aac':
        nacSeq(1)
    elif tech == 'DPC' or tech == 'dpc':
        nacSeq(2)
    elif tech == 'TPC' or tech == 'tpc':
        nacSeq(3)
    elif tech == 'KMER' or tech == 'kmer':
        k = int(input('Frequency of k-mer (e.g., 3, 4): '))
        findKmers(k)
    elif tech == 'RCKmer' or tech == 'rckmer':
        k = int(input('Frequency of k-mer (e.g., 3, 4): '))
        findKmers(k)
    elif tech == 'Kstep' or tech == 'kstep':
        window = int(input('Sliding Window (e.g., 30, 40): '))
        step = int(input('Window Step (e.g., 3, 4): '))
        k = int(input('k-mer Frequency in Window (e.g., 3, 4): '))
        seqKGAP(window, k, step)
#############################################################################