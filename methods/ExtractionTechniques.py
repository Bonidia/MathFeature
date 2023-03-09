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
    

# Nucleic acid composition (NAC), Di-nucleotide composition (DNC), Tri-nucleotide composition (TNC)
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


def reverse_complement(seq):
    return ''.join([complement[base] for base in seq[::-1]])


# k-mer, RCKmer
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

            # exemplo: se seq=TT e ksize=3
            # entao nao precisa fazer a contagem
            if k <= len(seq):
                for subseq in chunks_two(seq,k):
                    subseq = reverse_complement(subseq) if tech == 'RCKmer' or tech == 'rckmer' else subseq
                    if subseq in kmer:
                        # print(subseq)
                        kmer[subseq] = kmer[subseq] + 1
                    else:
                        kmer[subseq] = 1

            for key, value in kmer.items():
                # print (key)
                # print (value)
                try:
                    v = value/totalWindows
                except ZeroDivisionError:
                    v = 0.0
                probabilities.append([str(key), v])
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
    print('################### Feature Extraction: Several Techniques ########################')
    print('##########               Author: Robson Parmezan Bonidia                ###########')
    print('###################################################################################')
    print('\n')
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Fasta format file, E.g., test.fasta')
    parser.add_argument('-o', '--output', help='CSV format file, E.g., test.csv')
    parser.add_argument('-l', '--label', help='Dataset Label, E.g., lncRNA, mRNA, sncRNA ...')
    parser.add_argument('-t', '--type', help='type of Feature Extraction...')
    parser.add_argument('-seq', '--seq', help='type of sequence, 1 = DNA and 2 = RNA')
    args = parser.parse_args()
    finput = str(args.input)
    foutput = str(args.output)
    labelDataset = str(args.label)
    tech = str(args.type)
    seq = int(args.seq)
    stepw = 1
    if seq == 1:
        caracteres = ['A', 'C', 'G', 'T']
        complement = {'A': 'T','C': 'G','G': 'C','T': 'A'}
    else:
        caracteres = ['A', 'C', 'G', 'U']
        complement = {'A': 'U','C': 'G','G': 'C','U': 'A'}
    if tech == 'NAC' or tech == 'nac':
        nacSeq(1)
    elif tech == 'DNC' or tech == 'dnc':
        nacSeq(2)
    elif tech == 'TNC' or tech == 'tnc':
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