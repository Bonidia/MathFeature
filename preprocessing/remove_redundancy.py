#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import warnings
warnings.filterwarnings("ignore")
from Bio import SeqIO


def remove_equal_sequences(data):
    dataset = {}
    arq = open("prep_dataset.fasta", 'a')
    print("Dataset: " + data)
    for seq_record in SeqIO.parse(data, "fasta"):
        name = seq_record.name
        seq = seq_record.seq
        if seq in dataset.values():
            print("Removed Equal Sequences: %s" % (name))
        else:
            dataset[name] = seq
            arq.write(">" + name)
            arq.write("\n")
            arq.write(str(seq))
            arq.write("\n")
    return


#############################################################################    
if __name__ == "__main__":
    print("\n")
    print("###################################################################################")
    print("################    Feature Extraction: Removing Redundancy   #####################")
    print("################         Arguments: python3.5 -i dataset1      ####################")
    print("##########               Author: Robson Parmezan Bonidia                ###########")
    print("###################################################################################")
    print("\n")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input1', help='Fasta format file 1, E.g., dataset1.fasta')
    args = parser.parse_args()
    finput_one = str(args.input1)
    remove_equal_sequences(finput_one)
#############################################################################
