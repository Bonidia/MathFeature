#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import os
import time
import statistics
from Bio import SeqIO
from scipy.fftpack import fft
from gooey import Gooey, GooeyParser
import warnings
warnings.filterwarnings("ignore")
path = os.path.dirname(__file__)
if path != '':
    path = path + '/'


#############################################################################
#############################################################################

def header():
    dataset = open(foutput, 'a')
    dataset.write("nameseq,average,median,maximum,minimum,peak,"
                  + "none_levated_peak,sample_standard_deviation,population_standard_deviation,"
                  + "percentile15,percentile25,percentile50,percentile75,amplitude,"
                  + "variance,interquartile_range,semi_interquartile_range,"
                  + "coefficient_of_variation,skewness,kurtosis,label")
    dataset.write("\n")
    return


def file_record():
    dataset = open(foutput, 'a')
    dataset.write("%s," % (str(name_seq)))
    for metric in features:
        dataset.write("%s," % metric)
        # dataset.write("{0:.4f},".format(metric))
    dataset.write(label_dataset)
    dataset.write("\n")
    print('Recorded Sequence: %s' % name_seq)
    return


def feature_extraction():
    global features
    features = []
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


def binary_fourier():
    header()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        A = []
        C = []
        T = []
        G = []
        for nucle in seq:
            if nucle == "A":
                A.append(1)
            else:
                A.append(0)
            if nucle == "C":
                C.append(1)
            else:
                C.append(0)
            if nucle == "T" or nucle == "U":
                T.append(1)
            else:
                T.append(0)
            if nucle == "G":
                G.append(1)
            else:
                G.append(0)
        FA = fft(A)
        FC = fft(C)
        FT = fft(T)
        FG = fft(G)
        for i in range(len(seq)):
            specTotal = (abs(FA[i])**2) + (abs(FC[i])**2) + (abs(FT[i])**2) + (abs(FG[i])**2)
            specTwo = (abs(FA[i])) + (abs(FC[i])) + (abs(FT[i])) + (abs(FG[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        feature_extraction()
        file_record()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("A -- %s" % (A))
        # print("C -- %s" % (C))
        # print("T -- %s" % (T))
        # print("G -- %s" % (G))
        # print("\n")
        # print("A -- %s" % (abs(FA)))
        # print("C -- %s" % (abs(FC)))
        # print("T -- %s" % (abs(FT)))
        # print("G -- %s" % (abs(FG)))
        # print("\n")
        # print(spectrumTwo)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def zcurve_fourier():
    header()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        ###################################
        ###################################
        R = 0  # x[n] = (An + Gn) − (Cn + Tn) ≡ Rn − Yn
        Y = 0
        M = 0  # y[n] = (An + Cn) − (Gn + Tn) ≡ Mn − Kn
        K = 0
        W = 0  # z[n] = (An + Tn) − (Cn + Gn) ≡ Wn − Sn
        S = 0
        ###################################
        ###################################
        x = []
        y = []
        z = []
        for nucle in seq:
            if nucle == "A" or nucle == "G":
                R += 1
                x.append(R-Y)
            else:
                Y += 1
                x.append(R-Y)
            if nucle == "A" or nucle == "C":
                M += 1
                y.append(M-K)
            else:
                K += 1
                y.append(M-K)
            if nucle == "A" or nucle == "T" or nucle == "U":
                W += 1
                z.append(W-S)
            else:
                S += 1
                z.append(W-S)
        FX = fft(x)
        FY = fft(y)
        FZ = fft(z)
        for i in range(len(seq)):
            specTotal = (abs(FX[i])**2) + (abs(FY[i])**2) + (abs(FZ[i])**2)
            specTwo = (abs(FX[i])) + (abs(FY[i])) + (abs(FZ[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        feature_extraction()
        file_record()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print (seq)
        # print("\n")
        # print("X -- %s" % (x))
        # print("Y -- %s" % (y))
        # print("Z -- %s" % (z))
        # print("\n")
        # print("X -- %s" % (abs(FX)))
        # print("Y -- %s" % (abs(FY)))
        # print("Z -- %s" % (abs(FZ)))
        # print("\n")
        # print(spectrum)
        # print("\n")
        # print(spectrumTwo)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def integer_fourier():
    header()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        integer = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
                integer.append(0)
            elif nucle == "C":
                integer.append(1)
            elif nucle == "A":
                integer.append(2)
            else:
                integer.append(3)
        FI = fft(integer)
        for i in range(len(seq)):
            specTotal = (abs(FI[i])**2)
            specTwo = (abs(FI[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        feature_extraction()
        file_record()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("I -- %s" % (integer))
        # print("\n")
        # print("I -- %s" % (abs(FI)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def real_fourier():
    header()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        real = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
                real.append(1.5)
            elif nucle == "C":
                real.append(0.5)
            elif nucle == "A":
                real.append(-1.5)
            else:
                real.append(-0.5)
        FR = fft(real)
        for i in range(len(seq)):
            specTotal = (abs(FR[i])**2)
            specTwo = (abs(FR[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        feature_extraction()
        file_record()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("R -- %s" % (real))
        # print("\n")
        # print("R -- %s" % (abs(FR)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def eiip_fourier():
    header()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        eiip = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
                eiip.append(0.1335)
            elif nucle == "C":
                eiip.append(0.1340)
            elif nucle == "A":
                eiip.append(0.1260)
            else:
                eiip.append(0.0806)
        Feiip = fft(eiip)
        for i in range(len(seq)):
            specTotal = (abs(Feiip[i])**2)
            specTwo = (abs(Feiip[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        feature_extraction()
        file_record()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("R -- %s" % (eiip))
        # print("\n")
        # print("R -- %s" % (abs(Feiip)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def complex_number():
    header()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        complexNumber = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
                complexNumber.append(1-1j)
            elif nucle == "C":
                complexNumber.append(-1+1j)
            elif nucle == "A":
                complexNumber.append(1+1j)
            else:
                complexNumber.append(-1-1j)
        FcomplexNumber = fft(complexNumber)
        for i in range(len(seq)):
            specTotal = (abs(FcomplexNumber[i])**2)
            specTwo = (abs(FcomplexNumber[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        feature_extraction()
        file_record()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("R -- %s" % (complexNumber))
        # print("\n")
        # print("R -- %s" % (abs(FcomplexNumber)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
    return


def atomic_number():
    header()
    global spectrum, spectrumTwo, name_seq
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        spectrum = []
        spectrumTwo = []
        atomicNumber = []
        for nucle in seq:
            if nucle == "T" or nucle == "U":
                atomicNumber.append(66)
            elif nucle == "C":
                atomicNumber.append(58)
            elif nucle == "A":
                atomicNumber.append(70)
            else:
                atomicNumber.append(78)
        FatomicNumber = fft(atomicNumber)
        for i in range(len(seq)):
            specTotal = (abs(FatomicNumber[i])**2)
            specTwo = (abs(FatomicNumber[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        feature_extraction()
        file_record()
        ###################################
        ########## Debug ##################
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # spectrum.pop(0)
        # spectrum.pop(len(spectrum)-1)
        # print(seq)
        # print("\n")
        # print("R -- %s" % (atomicNumber))
        # print("\n")
        # print("R -- %s" % (abs(FcomplexNumber)))
        # print("\n")
        # print(spectrum)
        # fig = plt.figure()
        # ax = plt.axes()
        # ax.plot(spectrum)
        ###################################
        ###################################
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

    screen = parser.add_argument_group('Numerical Mapping with Fourier - DNA/RNA')
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
                        choices=['Binary', 'Z-curve', 'Real', 'Integer', 'EIIP', 'Complex Number', 'Atomic Number'],
                        nargs='+')
    screen.add_argument('-s',
                        '--seq',
                        metavar='Type of sequence',
                        help='Type of sequence, e.g., DNA or RNA',
                        required=True,
                        choices=['DNA', 'RNA'],
                        nargs='+')
    return parser.parse_args()


#############################################################################    
if __name__ == "__main__":
    args = main()
    finput = str(args.input_file)
    foutput = str(args.output_file)
    label_dataset = str(args.label)
    seq = args.seq
    seq = ''.join(seq)
    representation = args.representation
    representation = ''.join(representation)
    start_time = time.time()
    if representation == 'Binary':
        binary_fourier()
    elif representation == 'Z-curve':
        zcurve_fourier()
    elif representation == 'Real':
        real_fourier()
    elif representation == 'Integer':
        integer_fourier()
    elif representation == 'EIIP':
        eiip_fourier()
    elif representation == 'Complex Number':
        complex_number()
    elif representation == 'Atomic Number':
        atomic_number()
    else:
        print("This package does not contain this representation")
    print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
