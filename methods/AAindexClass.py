import argparse
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
import warnings
warnings.filterwarnings('ignore')


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


######################################################
######################################################
if __name__ == '__main__':
	print("\n")
	print("###################################################################################")
	print("#####################  Feature Extraction: AAindex Table   ########################")
	print("#####################    Author: Robson Parmezan Bonidia   ########################")
	print("###################################################################################")
	print("\n")
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', help='Fasta format file, E.g., test.fasta')
	parser.add_argument('-o', '--output', help='CSV format file, E.g., test.csv')
	parser.add_argument('-l', '--label', help='Dataset Label, E.g., lncRNA, mRNA, sncRNA ...')
	parser.add_argument('-t', '--table_aaindex', help='Table AAindex - e.g., see files directory...')
	args = parser.parse_args()
	finput = str(args.input)
	foutput = str(args.output)
	label_dataset = str(args.label)
	table = str(args.table_aaindex)
	start_time = time.time()
	mapping(finput, foutput, label_dataset, table)
	print('Computation time %s seconds' % (time.time() - start_time))
######################################################
######################################################
