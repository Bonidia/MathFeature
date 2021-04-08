"""
Modified from source code of PseKNC: A flexible web server for generating pseudo K-tuple nucleotide composition
downloaded from lin-group.cn/server/pseknc
article: https://doi.org/10.1016/j.ab.2014.04.001
"""


#############################################################################

from __future__ import division
import time
import os
import math
import re
import sys
import warnings
from Bio import SeqIO
from gooey import Gooey, GooeyParser
warnings.filterwarnings('ignore')
path = os.path.dirname(__file__)
if path != '':
    path = path + '/'


#############################################################################


def pseknc_record(listy, name_seq):
    i = 0
    z = 1
    listz = ''
    # Generating output file format
    while i < len(listy):
        if i < 1:
            listz = listz + str(listy[i])
            i = i + 1
        else:
            listz = listz + "," + str(listy[i])
            i = i + 1
    # Generating output file format
    output = name_seq + ',' + listz + ',' + label_dataset + '\n'
    # Writing to output file
    out_file = open(foutput, 'a')
    out_file.write(output)
    out_file.close()
    print('Recorded Sequence: %s' % name_seq)
    return


def header_pseknc(perms_list):
    file = open(foutput, 'a')
    file.write('%s,' % 'nameseq')
    i = 0
    # print(permsList)
    while i < len(perms_list):
        # print (perm)
        file.write('pseknc-' + str(i) + ',')
        i = i + 1
    file.write('label')
    file.write('\n')
    return


def header_pseknc2(perms_list):
    file = open(foutput, 'a')
    file.write('%s,' % 'nameseq')
    i = 0
    # print (permsList)
    while i < len(perms_list):
        # print (perm)
        file.write('pseknc-' + str(i) + ',')
        i = i + 1
    file.write('label')
    file.write('\n')
    return


#  generate_permutations
# ----------------------------------------
# inputs: chars = int
# outputs: all possible oligonucleotides of length k


def generate_permutations(chars):
    status = []
    for tmp in range(chars):
        status.append(0)
    last_char = len(allowed_chars)
    rows = []
    for x in range(last_char ** chars):
        rows.append("")
        for y in range(chars - 1, -1, -1):
            key = status[y]
            rows[x] = allowed_chars[key] + rows[x]
        for pos in range(chars - 1, -1, -1):
            if status[pos] == last_char - 1:
                status[pos] = 0
            else:
                status[pos] += 1
                break
    return rows

#  _mean: Calculates mean value of a list of numbers
# -----------------------------------------
# input: listy, a list of physicochemical property values
# output: the mean value of listy


def _mean(listy):
    return float(sum(listy))/len(listy)


def _std(listy, ddof=1):
    mean = _mean(listy)
    temp = [math.pow(i-mean, 2) for i in listy]
    res = math.sqrt(float(sum(temp))/len(listy))
    return res


def sepSequence(seq, k):
    i = k-1
    seqq = []
    while i < len(seq):
        j = 0
        nuc = ''
        while j < k:
            nuc = seq[i-j] + nuc
            j = j + 1
        seqq.append(nuc)
        i += 1
    return seqq


#  simplePseKNC: Calculates the frequencies of all possible
#  oligonucleotides in the input sequence
# ----------------------------------------
# inputs: inputFile = a string, outputFile = a string, k = int, 
#         formatt = string
# output: nothing, writing to a file


def simplePseKNC(input_file, output_file, k, formatt, geneticMaterial):
    listy = ''
    # Need to generate all possible oligonucleotides
    olinucz = generate_permutations(k)
    header_pseknc2(olinucz)
    for seq_record in SeqIO.parse(input_file, 'fasta'):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name
        listy = sepSequence(seq, k)
        z = 1
        out_file_name = output_file
        out_file = open(out_file_name, 'a')
        aux = 0
        listz = ''
        for olinuc in olinucz:
            freq = listy.count(olinuc)
            freq = int(freq/len(listy)*1000)/1000.0
            if olinucz.index(olinuc) < 1:
                # print(olinucz.index(olinuc))
                # print('olinuc: ',olinuc)
                listz = listz + str(freq)
            else:
                listz = listz + "," + str(freq)
        out_file.write(name_seq + ',')
        out_file.write(listz)
        out_file.write(',' + label_dataset + '\n')
        print('Recorded Sequence: %s' % name_seq)
    out_file.close()
    return


#  getValues: Returns a line of values for one property
#             from physicochemical property files
# -----------------------------------------
# input: prop = string of one property and supInfo = string
# output: values = a string representing a list of property values


def getValues(prop, sup_info):
    values = ""
    name = re.search(prop, sup_info)
    if name:
        strr = prop + '\s*\,(.+)'
        b = re.search(strr, sup_info)
        if b:
            values = b.group(1)
    return values


#  getSpecificValue: Returns a property value for a specific di- or tri-
#                    nucleotide
# -----------------------------------------
# input: olinuc = string, prop = string, supInfo = string
# output: value = an int that is property value of an olinuc


def getSpecificValue(olinuc, olinucs, prop, values, sup_info):
    values = values.split(",")
    count = olinucs.index(olinuc)
    value = values[count]
    return float(value)


#  hn: Hn function 
# -----------------------------------------
# inputs: olinuc = string, prop = string, supInfo = string
# output: temp = int


def hn(olinuc, olinucs, prop, sup_info, values):
    h0 = float(getSpecificValue(olinuc, olinucs, prop, values, sup_info))
    valueS = [float(x) for x in values.split(",")]
    temp = float((h0 - _mean(valueS)) / _std(valueS))
    return temp


#  theta2: Theta(i,i+j) function, for Type I
# -----------------------------------------
# input: seq = string, props = string, i = int,
# Type = int, k = int, j = int, supInfo = string
# output: summ = int


def theta2(seq, olinucs, props, i, k, j, sup_info):
    summ = 0
    values = ''
    for prop in props:
        values = getValues(prop, sup_info).rstrip()
        hn1 = hn(seq[i-1], olinucs, prop, sup_info, values)
        hn2 = hn(seq[i+j-1], olinucs, prop, sup_info, values)
        subsqr = math.pow(float(hn1-hn2), 2)
        summ = summ + subsqr
    return float(summ)/len(props)


#  J: J(i,i+j) function, for Type II
# --------------------------------------
# inputs: seq = string, prop = string, i = int,
#         k = int, j = int, supInfo = string
# output: product = int


def J(seq, olinucs, prop, i, k, j, sup_info):
    values = getValues(prop, sup_info)
    hn1 = hn(seq[i-1], olinucs, prop, sup_info, values)
    hn2 = hn(seq[i+j-1], olinucs, prop, sup_info, values)
    return float(hn1*hn2)


#  theta1: Theta(j) and Tau(LamGam) function
# -----------------------------------------
# input: seq = string, props = string, Type = int
#        k = int, j = int, supInfo = string
# output: final = int


def theta1(seq, olinucs, props, Type, k, j, sup_info):
    k = int(k)
    gamma = len(props)
    seqq = sepSequence(seq, k)
    i = 1
    a = 0
    if Type == 1:
        var = len(seq) - int(k) - int(j) + 1
        while i <= var:
            b = 0
            b = theta2(seqq, olinucs, props, i, k, j, sup_info)
            a = a + b
            i = i + 1
    else:
        ii = 0
        var = len(seq) - int(k) - int(j/gamma)
        while i <= var:
            b = 0
            if ii == gamma:
                ii = 0
                b = J(seqq, olinucs, props[ii], i, k, int(j/gamma), sup_info)
                a = a + b
            else:
                b = J(seqq, olinucs, props[ii], i, k, int(j/gamma), sup_info)
                a = a + b
            ii = ii + 1
            i = i + 1
    final = float(a)/var
    return final


#  pseKNCHelper: Creates list of adjusted frequency values for 4^k 
#  oligonucleotides and the lambda terms
# -----------------------------------------
# input: seq = string, props = string, Type = int, k = int, j = int, w = int
# output: freqs = list of ints


def pseKNCHelper(seq, olinucs, props, Type, k, j, w, genetic_material, sup_info):
    gamma = len(props)
    freqs = []
    seqq = sepSequence(seq, k)
    olinucs = olinucs.split(",")
    for olinuc in olinucs:
        freq = seqq.count(olinuc)
        freq = float(freq/len(seqq))
        total = 0
        i = 1
        if Type == 2:
            j = j/gamma
        while i <= j:
            total = total + theta1(seq, olinucs, props, Type, k, i, sup_info)
            i = i + 1
        total = float(freq/(1 + (float(w)*total)))
        total = int(total*1000) / 1000.0
        freqs.append(total)
    # Computing Lambda terms...
    fourK = math.pow(4, k)
    mu = fourK + 1
    while (fourK+1) <= mu <= (fourK + j):
        top = float(w) * theta1(seq, olinucs, props, Type, k, int(mu-fourK), sup_info)
        bottomTheta = 0
        bottom = 0
        i = 1
        while 1 <= i <= j:
            bottomTheta += theta1(seq, olinucs, props, Type, k, i, sup_info)
            i += 1
        bottom = 1 + (float(w) * bottomTheta)
        term = float(top / bottom)
        term = int(term * 1000) / 1000.0
        freqs.append(term)
        mu += 1
    return freqs

#  pseKNC: Opens input and output files, calls functions to calculate
#  values and writes outputs to the output file in appropriate format
# -----------------------------------------
# input: inputFile = string, outputFile = string, propNames = string,
#        Type = int, k = int, j = int, w = int, formatt = string
# output: nothing, write to output file


def pseKNC(input_file, output_file, props, Type, k, j, w, formatt, genetic_material):
    j = float(j)
    gamma = len(props)
    # Getting supporting info from files
    sup_file_fame = data_folder
    sup_file = open(sup_file_fame, 'r')
    sup_info = sup_file.read()
    o = re.search('Physicochemical properties\,(.+)\n', sup_info)
    olinucs = ''
    if o:		
        olinucs = o.group(1).rstrip()
    sup_file.close()
    # Calculating frequencies
    aux = 0
    for seq_record in SeqIO.parse(input_file, 'fasta'):
        seq = seq_record.seq
        seq = seq.upper()
        name_seq = seq_record.name	
        if Type == 1:
            listy = pseKNCHelper(seq, olinucs, props, Type, k, j, w, genetic_material, sup_info)
        else:
            listy = pseKNCHelper(seq, olinucs, props, Type, k, j*gamma, w, genetic_material, sup_info)
        if aux == 0:
            header_pseknc(listy)
            aux = 1
        pseknc_record(listy, name_seq)
    return


def parameters(data_folder, x, t, k, j, w, s, finput, foutput, seq):
    inputFile = finput
    outputFile = foutput
    if seq == 2:
        geneticMaterial = 'RNA'
    else:
        geneticMaterial = 'DNA' 

    # Getting list of properties (ex. Tilt, shift)
    properties = ''
    props = []
    propFileName = x
    propFile = open(propFileName, 'r')
    for line in propFile:
        e = re.search('(.+)', line)
        if e:
            properties = properties + e.group(1) + ","
    propFile.close()
    props = properties.split(",")
    props = props[0:len(props)-1]
    
    # Checking if type argument is 1 or 2
    if t == 1 or t == 2:
        Type = t
    else:
        print('The Type argument must be 1 or 2.')
        sys.exit()

    # Checking if k argument is 2 or 3
    if k == 2 or k == 3:
        kay = k
    else:
        print('The K argument must be 2 or 3.')
        sys.exit()

    # Checking if weight argument is between 0 and 1.0
    if w > 0 and w <= 1:
        weight = w
    else:
        print('The weight factor argument must be between 0.1 and 1.0.')
        sys.exit()
    
    # Checking if the lambda argument is a whole number and smaller than L-k
    InFileName = inputFile
    InFile = open(InFileName, 'r')
    listy = ''
    for line in InFile:
        g = re.search('\>', line)
        if g:
            label = line
        else:
            listy = listy + line.rstrip()
    InFile.close()
    ell = len(listy)
    if (float(j) % 1) == 0 and (1 <= float(j) < (ell - kay)):
        lam = j
    else:
        print('Lambda must be a whole number and smaller than the length of the query sequence minus the k-tuple number (lambda < L-k).')
        print('Length of query sequence: ', ell)
        print('k: ', kay)
        sys.exit()
    
    formatt = 'csv'
    SupFileName = data_folder
    # SupFile = open(SupFileName, 'r')
    tt = 0
    for prop in props:
        with open(SupFileName, 'r')as SupFile:
            for line in SupFile:
                t = re.search(prop, line)
                if t:
                    tt = 1
            if not(tt == 1):
                print('"' + prop + ' was not found in the supporting information file."')
                print('Please check that the k value and the query properties correspond.')
                sys.exit()
        SupFile.close()

    if s == 1:
        simplePseKNC(inputFile, outputFile, int(kay), formatt, geneticMaterial)
    else:
        if kay == 2 or kay == 3:
            kay = int(kay)
        else:
            print('The k-tuple argument must be 2 or 3 (dinucleotides or trinucleotides).')
            sys.exit()
        pseKNC(inputFile, outputFile, props, Type, kay, lam, weight, formatt, geneticMaterial)
    return


@Gooey(program_name='MathFeature',
       default_size=(1000, 800),
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
    screen = parser.add_argument_group('PseKNC - DNA/RNA')
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
    screen.add_argument('-x',
                        '--prop',
                        metavar='Properties - Names',
                        help='File containing list of properties to be used in calculations, e.g., see File directory',
                        required=True,
                        widget='FileChooser')
    screen.add_argument('-xp',
                        '--prop_values',
                        metavar='Properties - Values',
                        help='File containing list of properties (values) to be used in calculations, e.g., see File directory',
                        required=True,
                        widget='FileChooser')
    screen.add_argument('-t',
                        '--type',
                        metavar='PseKNC Type',
                        required=True,
                        choices=['Type 1 PseKNC', 'Type 2 PseKNC'],
                        nargs='+')
    screen.add_argument('-k',
                        '--kind',
                        metavar='Kind of oligonucleotide',
                        help='Dinucleotide or Trinucleotide',
                        required=True,
                        choices=['Dinucleotide', 'Trinucleotide'],
                        nargs='+')
    screen.add_argument('-j',
                        '--plambda',
                        metavar='lambda Parameter',
                        help='Must be smaller than the length of any query sequence, e.g., 1',
                        default=1,
                        required=True,
                        action='store')
    screen.add_argument('-w',
                        '--weight',
                        metavar='Weight Parameter',
                        help='It can be a value between (0, 1], e.g, 1.0',
                        default=1.0,
                        required=True,
                        action='store')
    screen.add_argument('-s',
                        '--frequency',
                        metavar='Frequency',
                        help='Frequency of each oligonucleotide in the input sequence, e.g., 2',
                        default=2,
                        required=True,
                        action='store')
    screen.add_argument('-seq',
                        '--seq',
                        metavar='Type of sequence',
                        help='Type of sequence, e.g., DNA, RNA',
                        required=True,
                        choices=['DNA', 'RNA'],
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
    x = str(args.prop)
    data_folder = str(args.prop_values)
    t_pse = args.type
    t_pse = ''.join(t_pse)
    t = 2
    if t_pse == 'Type 1 PseKNC':
        t = 1
    else:
        t = 2
    k_olig = args.kind
    k_olig = ''.join(k_olig)
    k = 0
    if k_olig == 'Dinucleotide':
        k = 2
    else:
        k = 3
    j = int(args.plambda)
    w = float(args.weight)
    s = int(args.frequency)
    seq = args.seq
    seq = ''.join(seq)
    label_dataset = str(args.label)
    if seq == 'DNA':
        allowed_chars = ["A", "C", "G", "T"]
    else:
        allowed_chars = ["A", "C", "G", "U"]
    start_time = time.time()
    parameters(data_folder, x, t, k, j, w, s, finput, foutput, seq)
    print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
