import time
import os
import numpy as np
import pandas as pd
import warnings
from gooey import Gooey, GooeyParser
from Bio import SeqIO
warnings.filterwarnings('ignore')
path = os.path.dirname(__file__)
if path != '':
	path = path + '/'


def header_features(table, foutput):
	file = open(foutput, 'a')
	file.write('nameseq,')
	for index, row in table.iterrows():
		file.write(str('{pname}.average,{pname}.sample_standard_deviation,{pname}.percentile25,'
				   '{pname}.percentile50,{pname}.percentile75,{pname}.interquartile_range,'
				   '{pname}.coefficient_of_variation,').format(pname = row['AccNo']))
	file.write('label')
	file.write('\n')
	return


def file_record(features, name_seq, label_dataset):
	dataset = open(foutput, 'a')
	dataset.write('%s,' % (str(name_seq)))
	for metric in features:
		dataset.write('%.5f,' % metric)
	dataset.write(label_dataset)
	dataset.write('\n')
	print('Recorded Sequence: %s' % name_seq)
	return


def feature_extraction(features, map_table):
	for index, row in map_table.iterrows():
		prop = np.array(row)

		average = sum(prop)/len(prop)
		features.append(average)

		standard_deviation = np.std(prop)  # standard deviation
		features.append(standard_deviation)

		percentile25 = np.percentile(prop, 25)
		features.append(percentile25)

		percentile50 = np.percentile(prop, 50)
		features.append(percentile50)

		percentile75 = np.percentile(prop, 75)
		features.append(percentile75)

		interquartile_range = np.percentile(prop, 75) - np.percentile(prop, 25)
		features.append(interquartile_range)

		if average != 0:
			coefficient_of_variation = standard_deviation/average
			features.append(coefficient_of_variation)
		else:
			features.append(0)

	return features


def mapping(finput, foutput, label_dataset, table):
	table = pd.read_csv(table, sep='\t')
	header_features(table, foutput)
	for seq_record in SeqIO.parse(finput, 'fasta'):
		seq = seq_record.seq
		seq = seq.upper()
		name_seq = seq_record.name
		map_table = pd.DataFrame()
		features = []
		i = 0
		for aa in seq:
			if aa in table:
				map_table[i] = table[aa]
				i += 1
		feature_extraction(features, map_table)
		file_record(features, name_seq, label_dataset)
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
	screen = parser.add_argument_group('AAindex Table - Protein')
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
	screen.add_argument('-t',
						'--AAindex',
						metavar='AAindex Table',
						help='AAindex Table, e.g., see Files directory',
						required=True,
						widget='FileChooser')
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
	finput = str(args.input_file)
	foutput = str(args.output_file)
	label_dataset = str(args.label)
	table = str(args.AAindex)
	start_time = time.time()
	mapping(finput, foutput, label_dataset, table)
	print('Computation time %s seconds' % (time.time() - start_time))
######################################################
######################################################
