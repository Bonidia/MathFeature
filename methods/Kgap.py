#!/usr/bin/env python
# -*- coding: utf-8 -*-


#############################################################################

import argparse
from Bio import SeqIO
from itertools import product
import warnings
warnings.filterwarnings('ignore')

#############################################################################


def kgap_record(name_seq):
    file = open(foutput, 'a')
    file.write('%s,' % name_seq)
    for x in number_kmers:
        # print (x)
        file.write('%s,' % (str(x[1])))
    file.write(labelDataset)
    file.write('\n')
    print('Recorded Sequence: %s' % name_seq)
    return


def chunks_kgap(seq, gap, before, after):
    seqlen = len(seq)
    chunks = []
    for i in range(0, seqlen):
        j = i + after + gap + before
        if j < seqlen:
            chunks.append(str(seq[i:j]))
            # print(str(seq[i:j]))
    return chunks


def header_kgap(permsList):
    file = open(foutput, 'a')
    file.write('%s,' % 'nameseq')
    for perm in permsList:
        file.write('%s,' % str(perm))
    file.write('label')
    file.write('\n')
    return


def compare(chunks, label):
    table = {}
    # total = len(chunks)
    for i in label:
        count = 0
        for j in chunks:
            aux = j[:before]
            aux2 = j[before + k:]
            aux3 = i[:before]
            aux4 = i[before + k:]
            if aux == aux3 and aux2 == aux4:
                count = count + 1

        # prob = count/total
        prob = count
        dic = {i: prob}
        table.update(dic)
    
    return table


def kgap():
    global number_kmers
    i = 0
    label = []
    gap = ''
    while i < k:
        gap = gap + '_'
        i = i + 1

    permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=before)]
    permsList2 = [''.join(str(i) for i in x) for x in product(caracteres, repeat=after)]

    for i in permsList:
        for j in permsList2:
            char = i + gap + j
            label.append(char)
    header_kgap(label)

    for seq_record in SeqIO.parse(finput, 'fasta'):
        number_kmers = []
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        # total = len(seq)
        chunks = chunks_kgap(seq, k, before, after)
        # print(chunks)
        table = compare(chunks, label)
        for key, value in table.items():
            number_kmers.append([str(key), value])
        kgap_record(name_seq)
    return


#############################################################################    
if __name__ == "__main__":
    print("\n")
    print("###################################################################################")
    print("######################## Feature Extraction: k-gap scheme #########################")
    print("##########   Arguments: python3.5 -i input -o output -l label -k kmer   ###########")
    print("#####   -seq type sequence -bef number before kgap -aft number after kgap   #######")
    print("###########               Author: Alvaro Pedroso Queiroz                ###########")
    print("###################################################################################")
    print("\n")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Fasta format file, E.g., test.fasta')
    parser.add_argument('-o', '--output', help='CSV format file, E.g., test.csv')
    parser.add_argument('-l', '--label', help='Dataset Label, E.g., lncRNA, mRNA, sncRNA ...')
    parser.add_argument('-k', '--kgap', help='Frequency of k-gap, E.g., 1, 2, 3, 4...')
    parser.add_argument('-bef', '--bef', help='Number before k-gap, E.g., 1, 2): ')
    parser.add_argument('-aft', '--aft', help='Number after k-gap, E.g., 1, 2): ')
    parser.add_argument('-seq', '--seq', help='type of sequence, 1 = DNA, 2 = RNA and 3 = Protein')
    args = parser.parse_args()
    finput = str(args.input)
    foutput = str(args.output)
    labelDataset = str(args.label)
    k = int(args.kgap)
    before = int(args.bef)
    after = int(args.aft)
    seq = int(args.seq)
    if seq == 1:
        caracteres = ["A", "C", "G", "T"]
        kgap()
    elif seq == 2:
        caracteres = ["A", "C", "G", "U"]
        kgap()
    elif seq == 3:
        caracteres = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        kgap()
    else:
        print('Check parameter: -seq, options: 1 = DNA, 2 = RNA and 3 = Protein')
#############################################################################
