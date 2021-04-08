#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import time
import numpy as np
import warnings
import statistics
import os
from gooey import Gooey, GooeyParser
from Bio import SeqIO
from scipy.fftpack import fft
warnings.filterwarnings('ignore')
path = os.path.dirname(__file__)
if path != '':
    path = path + '/'


#############################################################################
#############################################################################


def check_files():
    for name, label in dataset_labels.items():
        if os.path.exists(name):
            print('Dataset %s: Found File' % name)
            run = 1
        else:
            print('Dataset %s: File not exists' % name)
            run = 0
            break
    return run


def sequence_length():
    dlength = []
    for name, label in dataset_labels.items():
        length = 0
        for seq_record in SeqIO.parse(name, 'fasta'):
            seq = seq_record.seq
            if len(seq) > length:
                length = len(seq)
        dlength.append(length)
    # max_length = max(dlength)
    return max(dlength)

        
def file_record(name_seq, mapping, label):
    dataset = open(foutput, 'a')
    dataset.write('%s,' % (str(name_seq)))
    for map in mapping:
        dataset.write('%s,' % map)
        # dataset.write('{0:.4f},'.format(metric))
    dataset.write(label)
    dataset.write('\n')
    print('Recorded Sequence: %s' % name_seq)
    return


def chunksTwo(seq, win):
    seqlen = len(seq)
    for i in range(seqlen):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return


def accumulated_amino_frequency():
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, 'fasta'):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            mapping = []
            protein = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0,
                       'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}
            for i in range(len(seq)):
                if seq[i] in protein:
                    protein[seq[i]] += 1
                    mapping.append(protein[seq[i]] / (i + 1))
            padding = (max_length - len(mapping))
            mapping = np.pad(mapping, (0, padding), 'constant')
            file_record(name_seq, mapping, label)
    return


def kmer_frequency_protein(k):
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, "fasta"):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            mapping = []
            kmer = {}
            totalWindows = (len(seq) - k) + 1  # (L - k + 1)
            for subseq in chunksTwo(seq, k):
                # print(subseq)
                if subseq in kmer:
                    kmer[subseq] = kmer[subseq] + 1
                else:
                    kmer[subseq] = 1
            for subseq in chunksTwo(seq, k):
                # print(kmer[subseq])
                # print(kmer[subseq]/totalWindows)
                mapping.append(kmer[subseq]/totalWindows)
            padding = ((max_length - k + 1) - len(mapping))
            mapping = np.pad(mapping, (0, padding), 'constant')
            file_record(name_seq, mapping, label)
    return


def integer_mapping_protein():
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, 'fasta'):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            mapping = []
            protein = {'A': 1, 'C': 5, 'D': 4, 'E': 7, 'F': 14, 'G': 8, 'H': 9, 'I': 10, 'K': 12, 'L': 11,
                       'M': 13, 'N': 3, 'P': 15, 'Q': 6, 'R': 2, 'S': 16, 'T': 17, 'V': 20, 'W': 18, 'Y': 19}
            for i in range(len(seq)):
                if seq[i] in protein:
                    mapping.append(protein[seq[i]])
            padding = (max_length - len(mapping))
            mapping = np.pad(mapping, (0, padding), 'constant')
            file_record(name_seq, mapping, label)
    return


def eiip_mapping_protein():
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, 'fasta'):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            mapping = []
            protein = {'A': 0.3710, 'C': 0.08292, 'D': 0.12630, 'E': 0.00580, 'F': 0.09460,
                       'G': 0.0049, 'H': 0.02415, 'I': 0.0000, 'K': 0.37100, 'L': 0.0000,
                       'M': 0.08226, 'N': 0.00359, 'P': 0.01979, 'Q': 0.07606, 'R': 0.95930, 'S': 0.08292,
                       'T': 0.09408, 'V': 0.00569, 'W': 0.05481, 'Y': 0.05159}
            for i in range(len(seq)):
                if seq[i] in protein:
                    mapping.append(protein[seq[i]])
            padding = (max_length - len(mapping))
            mapping = np.pad(mapping, (0, padding), 'constant')
            file_record(name_seq, mapping, label)
    return


#############################################################################
#############################################################################


def header_fourier():
    dataset = open(foutput, 'a')
    dataset.write('nameseq,average,median,maximum,minimum,peak,'
                  + 'none_levated_peak,sample_standard_deviation,population_standard_deviation,'
                  + 'percentile15,percentile25,percentile50,percentile75,amplitude,'
                  + 'variance,interquartile_range,semi_interquartile_range,'
                  + 'coefficient_of_variation,skewness,kurtosis,label')
    dataset.write('\n')
    return


def file_record_fourier(features, name_seq, label_dataset):
    dataset = open(foutput, 'a')
    dataset.write('%s,' % (str(name_seq)))
    for metric in features:
        dataset.write('%s,' % metric)
        # dataset.write('{0:.4f},'.format(metric))
    dataset.write(label_dataset)
    dataset.write('\n')
    print('Recorded Sequence: %s' % name_seq)
    return


def feature_extraction(features, spectrum, spectrumTwo):
    average = sum(spectrum)/len(spectrum)
    features.append(average)
    ###################################
    median = np.median(spectrum)
    features.append(median)
    ###################################
    maximum = np.max(spectrum)
    features.append(maximum)
    ###################################
    minimum = np.min(spectrum)
    features.append(minimum)
    ###################################
    peak = (len(spectrum)/3)/average
    features.append(peak)
    ###################################
    peak_two = (len(spectrumTwo)/3)/(np.mean(spectrumTwo))
    features.append(peak_two)
    ###################################
    standard_deviation = np.std(spectrum)  # standard deviation
    features.append(standard_deviation)
    ###################################
    standard_deviation_pop = statistics.stdev(spectrum)  # population sample standard deviation
    features.append(standard_deviation_pop)
    ###################################
    percentile15 = np.percentile(spectrum, 15)
    features.append(percentile15)
    ###################################
    percentile25 = np.percentile(spectrum, 25)
    features.append(percentile25)
    ###################################
    percentile50 = np.percentile(spectrum, 50)
    features.append(percentile50)
    ###################################
    percentile75 = np.percentile(spectrum, 75)
    features.append(percentile75)
    ###################################
    amplitude = maximum - minimum
    features.append(amplitude)
    ###################################
    # mode = statistics.mode(spectrum)
    ###################################
    variance = statistics.variance(spectrum)
    features.append(variance)
    ###################################
    interquartile_range = np.percentile(spectrum, 75) - np.percentile(spectrum, 25)
    features.append(interquartile_range)
    ###################################
    semi_interquartile_range = (np.percentile(spectrum, 75) - np.percentile(spectrum, 25))/2 
    features.append(semi_interquartile_range)
    ###################################
    coefficient_of_variation = standard_deviation/average
    features.append(coefficient_of_variation)
    ###################################
    skewness = (3 * (average - median))/standard_deviation
    features.append(skewness)   
    ###################################
    kurtosis = (np.percentile(spectrum, 75) - np.percentile(spectrum, 25)) / (2 * (np.percentile(spectrum, 90) - np.percentile(spectrum, 10))) 
    features.append(kurtosis)
    ###################################
    return


def accumulated_amino_frequency_fourier():
    header_fourier()
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, 'fasta'):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            features = []
            spectrum = []
            spectrumTwo = []
            mapping = []
            protein = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0,
                       'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}
            for i in range(len(seq)):
                if seq[i] in protein:
                    protein[seq[i]] += 1
                    mapping.append(protein[seq[i]] / (i + 1))
            Fmap = fft(mapping)
            for i in range(len(mapping)):
                specTotal = (abs(Fmap[i])**2)
                specTwo = (abs(Fmap[i]))
                spectrum.append(specTotal)
                spectrumTwo.append(specTwo)
            feature_extraction(features, spectrum, spectrumTwo)
            file_record_fourier(features, name_seq, label)
    return


def kmer_frequency_protein_fourier(k):
    header_fourier()
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, "fasta"):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            features = []
            spectrum = []
            spectrumTwo = []
            mapping = []
            kmer = {}
            totalWindows = (len(seq) - k) + 1  # (L - k + 1)
            for subseq in chunksTwo(seq, k):
                # print(subseq)
                if subseq in kmer:
                    kmer[subseq] = kmer[subseq] + 1
                else:
                    kmer[subseq] = 1
            for subseq in chunksTwo(seq, k):
                # print(kmer[subseq])
                # print(kmer[subseq]/totalWindows)
                mapping.append(kmer[subseq]/totalWindows)
            Fmap = fft(mapping)
            for i in range(len(mapping)):
                specTotal = (abs(Fmap[i])**2)
                specTwo = (abs(Fmap[i]))
                spectrum.append(specTotal)
                spectrumTwo.append(specTwo)
            feature_extraction(features, spectrum, spectrumTwo)
            file_record_fourier(features, name_seq, label)
    return


def integer_mapping_protein_fourier():
    header_fourier()
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, 'fasta'):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            features = []
            spectrum = []
            spectrumTwo = []
            mapping = []
            protein = {'A': 1, 'C': 5, 'D': 4, 'E': 7, 'F': 14, 'G': 8, 'H': 9, 'I': 10, 'K': 12, 'L': 11,
                       'M': 13, 'N': 3, 'P': 15, 'Q': 6, 'R': 2, 'S': 16, 'T': 17, 'V': 20, 'W': 18, 'Y': 19}
            for i in range(len(seq)):
                if seq[i] in protein:
                    mapping.append(protein[seq[i]])
            Fmap = fft(mapping)
            for i in range(len(mapping)):
                specTotal = (abs(Fmap[i])**2)
                specTwo = (abs(Fmap[i]))
                spectrum.append(specTotal)
                spectrumTwo.append(specTwo)
            feature_extraction(features, spectrum, spectrumTwo)
            file_record_fourier(features, name_seq, label)
    return


def eiip_mapping_protein_fourier():
    header_fourier()
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, 'fasta'):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            features = []
            spectrum = []
            spectrumTwo = []
            mapping = []
            protein = {'A': 0.3710, 'C': 0.08292, 'D': 0.12630, 'E': 0.00580, 'F': 0.09460,
                       'G': 0.0049, 'H': 0.02415, 'I': 0.0000, 'K': 0.37100, 'L': 0.0000,
                       'M': 0.08226, 'N': 0.00359, 'P': 0.01979, 'Q': 0.07606, 'R': 0.95930, 'S': 0.08292,
                       'T': 0.09408, 'V': 0.00569, 'W': 0.05481, 'Y': 0.05159}
            for i in range(len(seq)):
                if seq[i] in protein:
                    mapping.append(protein[seq[i]])
            Fmap = fft(mapping)
            for i in range(len(mapping)):
                specTotal = (abs(Fmap[i])**2)
                specTwo = (abs(Fmap[i]))
                spectrum.append(specTotal)
                spectrumTwo.append(specTwo)
            feature_extraction(features, spectrum, spectrumTwo)
            file_record_fourier(features, name_seq, label)
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

    screen = parser.add_argument_group('Numerical Mapping with Fourier - Protein')
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
                        '--representation',
                        metavar='Numerical Mapping with Fourier Transform',
                        required=True,
                        choices=['ANF', 'k-mer frequency', 'Integer', 'EIIP'],
                        nargs='+')
    screen.add_argument('-k',
                        '--pattern',
                        metavar='Frequency of k-mer - Only for k-mer frequency option',
                        help='Frequency of k-mer, e.g., 2, 3, 20...',
                        default=3,
                        action='store')
    return parser.parse_args()


#############################################################################
#############################################################################   
if __name__ == '__main__':
    args = main()
    finput = str(args.input_file)
    foutput = str(args.output_file)
    label_dataset = str(args.label)
    kmer = int(args.pattern)
    representation = args.representation
    representation = ''.join(representation)
    dataset_labels = {}
    dataset_labels[finput] = label_dataset
    start_time = time.time()
    if check_files() == 1:
        max_length = sequence_length()
        if representation == 'ANF':
            accumulated_amino_frequency_fourier()
        elif representation == 'k-mer frequency':
            kmer_frequency_protein_fourier(kmer)
        elif representation == 'Integer':
            integer_mapping_protein_fourier()
        elif representation == 'EIIP':
            eiip_mapping_protein_fourier()
        else:
            print('This package does not contain this approach - Check parameter -r')
    else:
        print('Some file does not exist')
    print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
#############################################################################
