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
    file.write('%s,' % name_seq)
    for x in probabilities:
        # print (x)
        file.write('%s,' % str(x[1]))
    file.write(labelDataset)
    file.write('\n')
    print('Recorded Sequence: %s' % name_seq)
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
    file.write('%s,' % 'nameseq')
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
                subseq = reverse_complement(subseq) if tech == 'RCKmer' or tech == 'rckmer' else subseq
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
    print('Recorded Sequence: %s' % name_seq)
    return


# xxKGAP - Correct with step - new approach - my idea
def seqKGAP(window, ksize, step):
    global number_kmers
    header_kmer(ksize)
    for seq_record in SeqIO.parse(finput, 'fasta'):
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

    screen = parser.add_argument_group('Other Descriptors - DNA/RNA')
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
                        choices=['NAC', 'DNC', 'TNC', 'kmer'],
                        nargs='+')
    screen.add_argument('-k',
                        '--pattern',
                        metavar='Frequency of k-mer - Only for kmer option',
                        help='Frequency of k-mer, e.g., 2, 3, 20...',
                        default=3,
                        action='store')
    screen.add_argument('-s',
                        '--seq',
                        metavar='Type of sequence',
                        help='Type of sequence, e.g., DNA or RNA',
                        required=True,
                        choices=['DNA', 'RNA'],
                        nargs='+')
    return parser.parse_args()


#############################################################################    
if __name__ == '__main__':
    args = main()
    finput = str(args.input_file)
    foutput = str(args.output_file)
    labelDataset = str(args.label)
    tech = args.type
    tech = ''.join(tech)
    seq = args.seq
    seq = ''.join(seq)
    ksize = int(args.pattern)
    stepw = 1
    start_time = time.time()
    if seq == 'DNA':
        caracteres = ['A', 'C', 'G', 'T']
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    else:
        caracteres = ['A', 'C', 'G', 'U']
        complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    if tech == 'NAC' or tech == 'nac':
        nacSeq(1)
    elif tech == 'DNC' or tech == 'dnc':
        nacSeq(2)
    elif tech == 'TNC' or tech == 'tnc':
        nacSeq(3)
    elif tech == 'KMER' or tech == 'kmer':
        findKmers(ksize)
    elif tech == 'RCKmer' or tech == 'rckmer':
        findKmers(ksize)
    elif tech == 'Kstep' or tech == 'kstep':
        window = int(input('Sliding Window (e.g., 30, 40): '))
        step = int(input('Window Step (e.g., 3, 4): '))
        k = int(input('k-mer Frequency in Window (e.g., 3, 4): '))
        seqKGAP(window, k, step)
    print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
