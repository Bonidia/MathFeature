#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import time
import warnings
import math
from gooey import Gooey, GooeyParser
from Bio import SeqIO
warnings.filterwarnings("ignore")
path = os.path.dirname(__file__)
if path != '':
    path = path + '/'


def header(ksize):
    file = open(foutput, 'a')
    file.write("nameseq,")
    for i in range(1, ksize+1):
        file.write("k" + str(i) + ",")
    file.write("label")
    file.write("\n")
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

            
def file_record():
    file = open(foutput, 'a')
    file.write("%s," % name_seq)
    for data in information_entropy:
        file.write("%s," % str(data))
    file.write(label_dataset)
    file.write("\n")
    print("Recorded Sequence!!!")
    return
    

def entropy_equation():
    header(ksize)
    global name_seq, information_entropy
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        information_entropy = []
        for k in range(1, ksize+1):
            probabilities = []
            kmer = {}
            total_windows = (len(seq) - k) + 1 # (L - k + 1)
            for subseq in chunks_two(seq, k):
                if subseq in kmer:
                    # print(subseq)
                    kmer[subseq] = kmer[subseq] + 1
                else:
                    kmer[subseq] = 1
            for key, value in kmer.items():
                # print(key)
                # print(value)
                probabilities.append(value/total_windows)
            if e == "Shannon" or e == "shannon":
                entropy_equation = [(p * math.log(p, 2)) for p in probabilities]
                entropy = -(sum(entropy_equation))
                information_entropy.append(entropy)
            else:
                q = 2
                entropy_equation = [(p ** q) for p in probabilities]
                entropy =  (1/(q - 1)) * (1 - sum(entropy_equation))
                information_entropy.append(entropy)
        file_record()
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

    screen = parser.add_argument_group('Entropy Descriptor')
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
    screen.add_argument('-k',
                        '--pattern',
                        metavar='Frequency of k-mer',
                        help='Frequency of k-mer, e.g., 2, 3, 20...',
                        default=10,
                        required=True,
                        action='store')
    screen.add_argument('-e',
                        '--entropy',
                        metavar='Type of Entropy',
                        help='Type of Entropy, E.g., Shannon or Tsallis.',
                        required=True,
                        choices=['Shannon', 'Tsallis'],
                        nargs='+',
                        default='Shannon')
    return parser.parse_args()


#############################################################################    
if __name__ == "__main__":
    args = main()
    finput = str(args.input_file)
    foutput = str(args.output_file)
    label_dataset = str(args.label)
    ksize = int(args.pattern)
    stepw = 1
    e = args.entropy
    e = ''.join(e)
    start_time = time.time()
    entropy_equation()
    print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
