#!/usr/bin/env python
# -*- coding: utf-8 -*-


#############################################################################

import os
import time
import warnings
import collections
from Bio import SeqIO
from itertools import product
from gooey import Gooey, GooeyParser
warnings.filterwarnings('ignore')
path = os.path.dirname(__file__)
if path != '':
    path = path + '/'


#############################################################################


def headerNAC(k):
    global listTotal, caracteres
    listTotal = []
    file = open(foutput, 'a')
    file.write('%s,' % 'nameseq')
    permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
    # print (permsList)
    for perm in permsList:
        # print (perm)
        listTotal.append(perm)
        file.write('%s,' % str(perm))
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
    print('Recorded Sequence: %s' % name_seq)
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
        totalWindows = (len(seq) - k) + 1  # (L - k + 1)
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
    global listTotal, caracteres
    listTotal = []
    file = open(foutput, 'a')
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
    for seq_record in SeqIO.parse(finput, 'fasta'):
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

    screen = parser.add_argument_group('Other Descriptors - Protein')
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
    screen.add_argument('-t',
                        '--type',
                        metavar='Type of Feature Extraction',
                        required=True,
                        choices=['AAC', 'DPC', 'TPC', 'kmer'],
                        nargs='+')
    screen.add_argument('-k',
                        '--pattern',
                        metavar='Frequency of k-mer - Only for kmer option',
                        help='Frequency of k-mer, e.g., 2, 3, 20...',
                        default=3,
                        action='store')
    return parser.parse_args()

        
#############################################################################    
if __name__ == '__main__':
    args = main()
    finput = str(args.input_file)
    foutput = str(args.output_file)
    labelDataset = str(args.label)
    tech = args.type
    tech = ''.join(tech)
    ksize = int(args.pattern)
    stepw = 1
    start_time = time.time()
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
        findKmers(ksize)
    elif tech == 'RCKmer' or tech == 'rckmer':
        findKmers(ksize)
    print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
