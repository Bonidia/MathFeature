#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import warnings
warnings.filterwarnings("ignore")
from Bio import SeqIO


def count_sequences(arq):
    i = 0
    for seq_record in SeqIO.parse(arq, "fasta"):
        name_seq = seq_record.name
        i += 1
    print ("Number of sequences: %s" % (i))
    print("Finished")


#############################################################################    
if __name__ == "__main__":
    print("\n")
    print("###################################################################################")
    print("#################### Feature Extraction: Sequences Counting #######################")
    print("########################   Arguments: python3.5 -i input  #########################")
    print("##########               Author: Robson Parmezan Bonidia                ###########")
    print("###################################################################################")
    print("\n")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Fasta format file, E.g., test.fasta')
    args = parser.parse_args()
    finput = str(args.input)
    count_sequences(finput)
#############################################################################