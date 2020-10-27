#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import re
import argparse
warnings.filterwarnings("ignore")
from Bio import SeqIO


def preprocessing(finput,foutput):
    alphabet = ("B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|V|W|X|Y|Z")
    file = open(foutput, 'a')
    for seq_record in SeqIO.parse(finput, "fasta"):
        name_seq = seq_record.name
        seq = seq_record.seq
        if re.search(alphabet, str(seq)) is not None:
            print(name_seq)
            print("Removed Sequence")
        else:
            file.write(">%s" % (str(name_seq)))
            file.write("\n")
            file.write(str(seq))
            file.write("\n")
            print(name_seq)
            print("Included Sequence")
    print("Finished")


#############################################################################    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Fasta format file, E.g., dataset.fasta')
    parser.add_argument('-o', '--output', help='Fasta format file, E.g., preprocessing.fasta')
    args = parser.parse_args()
    finput = str(args.input)
    foutput = str(args.output)
    preprocessing(finput,foutput)
#############################################################################]