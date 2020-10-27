#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import collections
import argparse
from Bio import SeqIO
from itertools import product


def perms():
    global listTotal, caracteres
    listTotal = []
    file = open(foutput, 'a')
    file.write("%s," % ("nameseq"))
    for k in range(1, ksize+1):
        permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
        # print (permsList)
        for perm in permsList:
            # print (perm)
            listTotal.append(perm)
            file.write("%s," % (str(perm)))
    file.write("label")
    file.write("\n")
    return


def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0,seqlen,step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return        
    

def chunksTwo(seq, win):
    seqlen = len(seq)
    for i in range(seqlen):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return

            
def fileRecord(name_seq):
    file = open(foutput, 'a')
    file.write("%s," % (name_seq))
    for x in probabilities:
        # print (x)
        file.write("%s," % (str(x[1])))
    file.write(labelDataset)
    file.write("\n")
    print ("Recorded Sequence: %s" % (name_seq))
    return
    

def findKmers():
    global probabilities
    perms()
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()	
        name_seq = seq_record.name
        probabilities = []
        for k in range(1, ksize+1):
            kmer = {}
            totalWindows = (len(seq) - k) + 1 # (L - k + 1)
            permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
            for key in permsList:
                kmer[key] = 0
            kmer = collections.OrderedDict(sorted(kmer.items()))
            for subseq in chunksTwo(seq, k):
                if subseq in kmer:
                    # print(subseq)
                    kmer[subseq] = kmer[subseq] + 1
                else:
                    kmer[subseq] = 1
            for key, value in kmer.items():
                # print (key)
                # print (value)
                probabilities.append([str(key), value/totalWindows])
        fileRecord(name_seq)
    return

        
#############################################################################    
if __name__ == "__main__":
    print("\n")
    print("###################################################################################")
    print("######################## Feature Extraction: k-mer scheme #########################")
    print("##########   Arguments: python3.5 -i input -o output -l label -k kmer   ###########")
    print("##########               Author: Robson Parmezan Bonidia                ###########")
    print("###################################################################################")
    print("\n")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Fasta format file, E.g., test.fasta')
    parser.add_argument('-o', '--output', help='CSV format file, E.g., test.csv')
    parser.add_argument('-l', '--label', help='Dataset Label, E.g., lncRNA, mRNA, sncRNA ...')
    parser.add_argument('-k', '--kmer', help='Range of k-mer, E.g., 1-mer (1) or 2-mer (1, 2) ...')
    parser.add_argument('-seq', '--seq', help='type of sequence, 1 = DNA and 2 = RNA')
    args = parser.parse_args()
    finput = str(args.input)
    foutput = str(args.output)
    labelDataset = str(args.label)
    ksize = int(args.kmer)
    seq = int(args.seq)
    stepw = 1
    if seq == 1:
    	caracteres = ["A", "C", "G", "T"]
    else:
    	caracteres = ["A", "C", "G", "U"]
    findKmers()
#############################################################################