#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import warnings
warnings.filterwarnings("ignore")
from Bio import SeqIO


def remove_equal_sequences(data1,data2):
    dataset1 = {}
    dataset2 = {}
    for seq_record in SeqIO.parse(data1, "fasta"):
        name = seq_record.name
        seq = seq_record.seq
        dataset1[name] = seq
    for seq_record in SeqIO.parse(data2, "fasta"):
        name = seq_record.name
        seq = seq_record.seq
        dataset2[name] = seq
    for i in range(1, 3):
        if i == 1:
            arq = open("new" + data1, 'a')
            print("Dataset: " + data1)
            for key, value in dataset1.items():
                if value in dataset2.values():
                    print("Removed Equal Sequences: %s" % (key))
                else:
                    arq.write(">" + key)
                    arq.write("\n")
                    arq.write(str(value))
                    arq.write("\n")
        else:
            arq = open("new" + data2, 'a')
            print("Dataset: " + data2)
            for key, value in dataset2.items():
                if value in dataset1.values():
                    print("Removed Equal Sequences: %s" % (key))
                else:
                    arq.write(">" + key)
                    arq.write("\n")
                    arq.write(str(value))
                    arq.write("\n")
    return


#############################################################################    
if __name__ == "__main__":
    print("\n")
    print("###################################################################################")
    print("################ Feature Extraction: Removing Equal Sequences #####################")
    print("################  Arguments: python3.5 -i dataset1 -t dataset2 ####################")
    print("##########               Author: Robson Parmezan Bonidia                ###########")
    print("###################################################################################")
    print("\n")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input1', help='Fasta format file 1, E.g., dataset1.fasta')
    parser.add_argument('-t', '--input2', help='Fasta format file 2, E.g., dataset2.fasta')
    args = parser.parse_args()
    finput_one = str(args.input1)
    finput_two = str(args.input2)
    remove_equal_sequences(finput_one,finput_two)
#############################################################################