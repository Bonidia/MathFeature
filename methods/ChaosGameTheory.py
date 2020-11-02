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
warnings.filterwarnings("ignore")


#############################################################################
#############################################################################


def check_files():
    for name, label in dataset_labels.items():
        if os.path.exists(name):
            print("Dataset %s: Found File" % (name))
            run = 1
        else:
            print("Dataset %s: File not exists" % (name))
            run = 0
            break
    return run


def sequence_length():
    dlength = []
    for name, label in dataset_labels.items():
        length = 0
        for seq_record in SeqIO.parse(name, "fasta"):
            seq = seq_record.seq
            if len(seq) > length:
                length = len(seq)
        dlength.append(length)
    # max_length = max(dlength)
    return max(dlength)

        
def file_record(name_seq, mapping, label):
    dataset = open(foutput, 'a')
    dataset.write("%s," % (str(name_seq)))
    for map in mapping:
        dataset.write("%s," % (map))
        # dataset.write("{0:.4f},".format(metric))
    dataset.write(label)
    dataset.write("\n")
    print("Recorded Sequence: %s" % (name_seq))
    return


def chunksTwo(seq, win):
    seqlen = len(seq)
    for i in range(seqlen):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return


def frequency_chaos(k):
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, "fasta"):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            mapping = []
            kmer = {}
            totalWindows = (len(seq) - k) + 1 # (L - k + 1)
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


def classifical_chaos():
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, "fasta"):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            Sx = []
            Sy = []
            for nucle in seq:
                if nucle == "A":
                    Sx.append(1)
                    Sy.append(1)
                elif nucle == "C":
                    Sx.append(-1)
                    Sy.append(-1)
                elif nucle == "T" or nucle =="U":
                    Sx.append(-1)
                    Sy.append(1)
                else:
                    Sx.append(1)
                    Sy.append(-1)
            CGR_x = [] 
            CGR_y = []
            for i in range(0, len(Sx)):
            	if i == 0:
                    CGR_x.append(0.5 * Sx[i])
                    CGR_y.append(0.5 * Sy[i])
            	else:
                    CGR_x.append(0.5 * Sx[i] + 0.5 * CGR_x[i - 1])
                    CGR_y.append(0.5 * Sy[i] + 0.5 * CGR_y[i - 1])           
            concat = CGR_x + CGR_y
            padding = (max_length - len(Sx)) * 2
            mapping = np.pad(concat, (0, padding), 'constant')
            file_record(name_seq, mapping, label)
    return


#############################################################################
#############################################################################


def header_fourier():
    dataset = open(foutput, 'a')
    dataset.write("nameseq,average,median,maximum,minimum,peak,"
                  + "none_levated_peak,sample_standard_deviation,population_standard_deviation,"
                  + "percentile15,percentile25,percentile50,percentile75,amplitude,"
                  + "variance,interquartile_range,semi_interquartile_range,"
                  + "coefficient_of_variation,skewness,kurtosis,label")
    dataset.write("\n")
    return


def file_record_fourier(features, name_seq, label_dataset):
    dataset = open(foutput, 'a')
    dataset.write("%s," % (str(name_seq)))
    for metric in features:
        dataset.write("%s," % (metric))
        # dataset.write("{0:.4f},".format(metric))
    dataset.write(label_dataset)
    dataset.write("\n")
    print("Recorded Sequence: %s" % (name_seq))
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


def classifical_chaos_fourier():
    header_fourier()
    for finput, label in dataset_labels.items():
        for seq_record in SeqIO.parse(finput, "fasta"):
            seq = seq_record.seq
            seq = seq.upper()
            name_seq = seq_record.name
            features = []
            spectrum = []
            spectrumTwo = []
            Sx = []
            Sy = []
            for nucle in seq:
                if nucle == "A":
                    Sx.append(1)
                    Sy.append(1)
                elif nucle == "C":
                    Sx.append(-1)
                    Sy.append(-1)
                elif nucle == "T" or nucle =="U":
                    Sx.append(-1)
                    Sy.append(1)
                else:
                    Sx.append(1)
                    Sy.append(-1)
            CGR_x = [] 
            CGR_y = []
            for i in range(0, len(Sx)):
                if i == 0:
                    CGR_x.append(0.5 * Sx[i])
                    CGR_y.append(0.5 * Sy[i])
                else:
                    CGR_x.append(0.5 * Sx[i] + 0.5 * CGR_x[i - 1])
                    CGR_y.append(0.5 * Sy[i] + 0.5 * CGR_y[i - 1])           
            Fx = fft(CGR_x)
            Fy = fft(CGR_y)
            for i in range(len(Fx)):
                specTotal = (abs(Fx[i])**2) + (abs(Fy[i])**2)
                specTwo = (abs(Fx[i])) + (abs(Fy[i]))
                spectrum.append(specTotal)
                spectrumTwo.append(specTwo)
            feature_extraction(features, spectrum, spectrumTwo)
            file_record_fourier(features, name_seq, label)
    return


def frequency_chaos_fourier(k):
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
            totalWindows = (len(seq) - k) + 1 # (L - k + 1)
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


#############################################################################
#############################################################################   
if __name__ == "__main__":
    print("\n")
    print("###################################################################################")
    print("##########            Feature Extraction: Chaos Game Theory             ###########")
    print("##########  Arguments: -i number of datasets -o output -r representation  #########")
    print("##########      -r:  1 = Classifical Chaos Game Representation            #########")
    print("##########      -r:  2 = Frequency Chaos Game Representation              #########")
    print("##########      -r:  3 = Chaos Game Signal with Classifical Representation  #######")
    print("##########      -r:  4 = Chaos Game Signal with Frequency                 #########")
    print("##########               Author: Robson Parmezan Bonidia                ###########")
    print("###################################################################################")
    print("\n")
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--number', help='Fasta format file | Number of dataset or labels')
    parser.add_argument('-o', '--output', help='Csv format file | E.g., train.csv')
    parser.add_argument('-r', '--approach', help='1 = Classifical Chaos Game Representation, 2 = Frequency Chaos Game Representation, 3 = Chaos Game Signal with Classifical Representation, 4 = Chaos Game Signal with Frequency')
    args = parser.parse_args()
    n = int(args.number)
    foutput = str(args.output)
    representation = int(args.approach)
    dataset_labels = {}
    for i in range(1, n + 1):
        name = input("Dataset %s: " % (i))
        label = input("Label for %s: " % (name))
        print("\n")
        dataset_labels[name] = label
    if check_files() == 1:
        max_length = sequence_length()
        if representation == 1:
            classifical_chaos()
        elif representation == 2:
            k = int(input("Frequency of k-mer (e.g., 3, 4): "))
            frequency_chaos(k)
        elif representation == 3:
            classifical_chaos_fourier()
        elif representation == 4:
            k = int(input("Frequency of k-mer (e.g., 3, 4): "))
            frequency_chaos_fourier(k)
        else:
            print("This package does not contain this approach - Check parameter -r")
    else:
        print("Some file not exists")
#############################################################################
#############################################################################
