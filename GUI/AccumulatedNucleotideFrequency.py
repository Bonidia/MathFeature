#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import numpy as np
from Bio import SeqIO
from itertools import product
from scipy.fftpack import fft, ifft
import warnings
import sys
import scipy.stats
import statistics
import os
import collections
warnings.filterwarnings('ignore')


#############################################################################
#############################################################################


def check_files():
    for name, label in dataset_labels.items():
        if os.path.exists(name):
            print('Dataset %s: Found File' % (name))
            run = 1
        else:
            print('Dataset %s: File not exists' % (name))
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
        dataset.write('%s,' % (map))
        # dataset.write('{0:.4f},'.format(metric))
    dataset.write(label)
    dataset.write('\n')
    print('Recorded Sequence: %s' % (name_seq))
    return


def accumulated_nucle_frequency():
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, 'fasta'):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            mapping = []
            A = 0
            C = 0
            T = 0
            G = 0
            for i in range(len(seq)):
                if seq[i] == 'A':
                    A += 1
                    mapping.append(A / (i + 1))
                elif seq[i] == 'C':
                    C += 1
                    mapping.append(C / (i + 1))
                elif seq[i] == 'T' or seq[i] == 'U':
                    T += 1
                    mapping.append(T / (i + 1))
                else:
                    G += 1
                    mapping.append(G / (i + 1))
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
        dataset.write('%s,' % (metric))
        # dataset.write('{0:.4f},'.format(metric))
    dataset.write(label_dataset)
    dataset.write('\n')
    print('Recorded Sequence: %s' % (name_seq))
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
    peak = (len(spectrum)/3)/(average)
    features.append(peak)
    ###################################
    peak_two = (len(spectrumTwo)/3)/(np.mean(spectrumTwo))
    features.append(peak_two)
    ###################################
    standard_deviation = np.std(spectrum) # standard deviation
    features.append(standard_deviation)
    ###################################
    standard_deviation_pop = statistics.stdev(spectrum) # population sample standard deviation 
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


def accumulated_nucle_frequency_fourier():
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
            A = 0
            C = 0
            T = 0
            G = 0
            for i in range(len(seq)):
                if seq[i] == 'A':
                    A += 1
                    mapping.append(A / (i + 1))
                elif seq[i] == 'C':
                    C += 1
                    mapping.append(C / (i + 1))
                elif seq[i] == 'T' or seq[i] == 'U':
                    T += 1
                    mapping.append(T / (i + 1))
                else:
                    G += 1
                    mapping.append(G / (i + 1))
            Fmap = fft(mapping)
            for i in range(len(mapping)):
                specTotal = (abs(Fmap[i])**2)
                specTwo = (abs(Fmap[i]))
                spectrum.append(specTotal)
                spectrumTwo.append(specTwo)
            feature_extraction(features, spectrum, spectrumTwo)
            file_record_fourier(features, name_seq, label)
    return


#############################################################################
#############################################################################   
if __name__ == '__main__':
    print('\n')
    print('###################################################################################')
    print('##########            Feature Extraction: Accumulated Nucle Frequency   ###########')
    print('##########  Arguments: -i number of datasets -o output -r representation  #########')
    print('##########            -r:  1 = Accumulated Nucle Frequency                #########')
    print('##########     -r:  2 = Accumulated Nucle Frequency with Fourier          #########')
    print('##########               Author: Robson Parmezan Bonidia                ###########')
    print('###################################################################################')
    print('\n')
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number', help='Fasta format file | Number of dataset or labels')
    parser.add_argument('-o', '--output', help='Csv format file | E.g., train.csv')
    parser.add_argument('-r', '--approach', help='1 = Accumulated Nucle Frequency, 2 = Accumulated Nucle Frequency with Fourier')
    args = parser.parse_args()
    n = int(args.number)
    foutput = str(args.output)
    representation = int(args.approach)
    dataset_labels = {}
    for i in range(1, n + 1):
        name = input('Dataset %s: ' % (i))
        label = input('Label for %s: ' % (name))
        print('\n')
        dataset_labels[name] = label
    if check_files() == 1:
        max_length = sequence_length()
        if representation == 1:
            accumulated_nucle_frequency()
        elif representation == 2:
            accumulated_nucle_frequency_fourier()
        else:
            print('This package does not contain this approach - Check parameter -r')
    else:
        print('Some file not exists')
#############################################################################
#############################################################################
