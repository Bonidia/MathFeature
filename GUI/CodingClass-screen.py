"""
Modified from source code of LGC 1.0,
downloaded from https://bigd.big.ac.cn/lgc/calculator
article: https://academic.oup.com/bioinformatics/article-abstract/35/17/2949/5288512?redirectedFrom=fulltext
"""


import operator
import warnings
import os
import time
import numpy as np
import scipy.stats
from gooey import Gooey, GooeyParser
from Bio import SeqIO
warnings.filterwarnings("ignore")
path = os.path.dirname(__file__)
if path != '':
	path = path + '/'


def des_start_code(codes):
	if codes in ('ATG', 'AUG'):
		return True
	return False


def des_end_code(codes):
	if codes in ('TAA', 'UAA', 'TAG', 'UAG', 'TGA', 'UGA'):
		return True
	return False


def read_by_three(string, offset):
	flag = True
	length = len(string)
	start = end = -1
	i = 0
	result = set()
	while i < length-2:
		codes = string[i:i+3]
		if des_start_code(codes) and flag:
			start = i
			flag = False
		if des_end_code(codes) and not flag:
			end = i + 2
			flag = True
		if (end > start) and (start != -1):
			result.add((start + offset, end + offset))
		i = i + 3
	return result


def get_gc(string):
	gc = ((string.count('G') + string.count('C')) / len(string)) * 100
	return gc


def get_info(string, pos):
	length = pos[1] - pos[0] + 1
	gc = get_gc(string[pos[0]:pos[1]+1])
	return str(pos[0]), str(pos[1]), str(length), str(gc)


def orf(seq):
	result_info = []
	strings = [seq, seq[1:], seq[2:]]
	for index, string in enumerate(strings):
		# print(index)
		# print(string)
		positions = read_by_three(string, index)
		positions = sorted(positions, key=operator.itemgetter(0))
		# print(positions)
		for pos in positions:
			result_info.append(get_info(seq, pos))
	# print(result_info)
	# print(len(result_info))
	return result_info


def run(finput, foutput, label_dataset):
	file = open(foutput, 'a')
	file.write('nameseq,maximum_ORF_length,minimum_ORF_length,std_ORF_length,average_ORF_length,cv_ORF_length,' + ''
			   + 'maximum_GC_content_ORF,minimum_GC_content_ORF,std_GC_content_ORF,' + ''
			   + 'average_GC_content_ORF,cv_GC_content_ORF,label')
	file.write('\n')
	for seq_record in SeqIO.parse(finput, 'fasta'):
		seq = seq_record.seq
		seq = seq.upper()
		name_seq = seq_record.name
		file.write('%s,' % name_seq)
		measures = orf(seq)
		if len(measures) > 0:
			length_orf = []
			gc_mea = []
			for values in measures:
				length_orf.append(int(values[2]))
				gc_mea.append(float(values[3]))
			# print(length_orf)
			# print(gc_mea)
			file.write('%s,' % max(length_orf))
			file.write('%s,' % min(length_orf))
			file.write('%s,' % np.std(length_orf))
			file.write('%s,' % np.mean(length_orf))
			file.write('%s,' % scipy.stats.variation(length_orf))
			file.write('%s,' % max(gc_mea))
			file.write('%s,' % min(gc_mea))
			file.write('%s,' % np.std(gc_mea))
			file.write('%s,' % np.mean(gc_mea))
			file.write('%s,' % scipy.stats.variation(gc_mea))
			file.write('%s' % label_dataset)
			file.write('\n')
		else:
			file.write('0,0,0,0,0,0,0,0,0,0,')
			file.write('%s' % label_dataset)
			file.write('\n')
		print('Recorded Sequence: %s' % name_seq)
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

	screen = parser.add_argument_group('ORF Descriptor')
	# test = parser.add_argument_group('ORF Descriptor')
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
	return parser.parse_args()


######################################################
######################################################
if __name__ == '__main__':
	args = main()
	ffinput = str(args.input_file)
	ffoutput = str(args.output_file)
	label = str(args.label)
	start_time = time.time()
	run(ffinput, ffoutput, label)
	print('Computation time %s seconds' % (time.time() - start_time))
######################################################
######################################################
