#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import re
import os
import time
from Bio import SeqIO
from gooey import Gooey, GooeyParser
warnings.filterwarnings("ignore")
path = os.path.dirname(__file__)
if path != '':
    path = path + '/'


def preprocessing(finput, foutput):
    alphabet = "B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|V|W|X|Y|Z"
    file = open(foutput, 'a')
    for seq_record in SeqIO.parse(finput, "fasta"):
        name_seq = seq_record.name
        seq = seq_record.seq
        if re.search(alphabet, str(seq)) is not None:
            print(name_seq)
            print("Removed Sequence")
        else:
            file.write(">%s" % (str(name_seq)))
            file.write("\n")
            file.write(str(seq))
            file.write("\n")
            print(name_seq)
            print("Included Sequence")
    print("Finished")


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
    screen = parser.add_argument_group('Preprocessing (To eliminate noise from sequences) - DNA/RNA')
    screen.add_argument('-i',
                        '--input_file',
                        metavar='Input File',
                        help='Fasta format file, e.g., test.fasta',
                        required=True,
                        widget='FileChooser')
    screen.add_argument('-o',
                        '--output_file',
                        metavar='Output File',
                        help='Fasta format file, e.g., preprocessing.fasta',
                        required=True,
                        widget='FileSaver')
    return parser.parse_args()


#############################################################################    
if __name__ == "__main__":
    args = main()
    finput = str(args.input_file)
    foutput = str(args.output_file)
    start_time = time.time()
    preprocessing(finput, foutput)
    print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
