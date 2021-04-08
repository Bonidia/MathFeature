#!/usr/bin/env python
# -*- coding: utf-8 -*-


#############################################################################

import time
import os
import warnings
from Bio import SeqIO
from itertools import product
from gooey import Gooey, GooeyParser
warnings.filterwarnings('ignore')
path = os.path.dirname(__file__)
if path != '':
    path = path + '/'

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
    screen = parser.add_argument_group('Xmer k-Spaced Ymer composition frequency (kGap) - DNA/RNA/Protein')
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
    screen.add_argument('-k',
                        '--kgap',
                        metavar='Frequency of k-gap',
                        help='Frequency of k-gap, e.g., 1 = A_A, 2 = A__A, 3 = A___A...',
                        required=True,
                        action='store')
    screen.add_argument('-bef',
                        '--bef',
                        metavar='Before k-gap',
                        help='before - e.g., 1 = A_A, 2 = AA_A, 3 = AAA_A...',
                        required=True,
                        action='store')
    screen.add_argument('-aft',
                        '--aft',
                        metavar='After k-gap',
                        help='after - e.g., 1 = A_A, 2 = A_AA, 3 = A_AAA...',
                        required=True,
                        action='store')
    screen.add_argument('-s',
                        '--seq',
                        metavar='Type of sequence',
                        help='Type of sequence, e.g., DNA, RNA, Protein',
                        required=True,
                        choices=['DNA', 'RNA', 'Protein'],
                        nargs='+')
    screen.add_argument('-l',
                        '--label',
                        metavar='Dataset Label',
                        help='Dataset Label, e.g., 0, 1, lncRNA, mRNA, sncRNA.',
                        required=True,
                        action='store')
    return parser.parse_args()


#############################################################################    
if __name__ == "__main__":
    args = main()
    finput = str(args.input_file)
    foutput = str(args.output_file)
    labelDataset = str(args.label)
    k = int(args.kgap)
    before = int(args.bef)
    after = int(args.aft)
    seq = args.seq
    seq = ''.join(seq)
    start_time = time.time()
    if seq == 'DNA':
        caracteres = ["A", "C", "G", "T"]
        kgap()
    elif seq == 'RNA':
        caracteres = ["A", "C", "G", "U"]
        kgap()
    elif seq == 'Protein':
        caracteres = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        kgap()
    else:
        print('Check parameter: -seq, options: 1 = DNA, 2 = RNA and 3 = Protein')
    print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
