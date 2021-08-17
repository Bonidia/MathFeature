#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import warnings
import sys
import os
warnings.filterwarnings("ignore")


def check_files():
	for name in dataset_names:
		if os.path.exists(name):
			print("Dataset %s: Found File" % (name))
			run = 1
		else:
			print("Dataset %s: File not exists" % (name))
			run = 0
			break
	return run


def concatenate():
	dataframes = pd.concat([pd.read_csv(f) for f in dataset_names], axis=1)
	dataframes = dataframes.loc[:, ~dataframes.columns.duplicated()]
	label = dataframes.pop('label')
	dataframes.insert(len(dataframes.columns), 'label', label)
	dataframes.to_csv(foutput, index=False)
	print("The files were concatenate")


#############################################################################    
if __name__ == "__main__":
	print("\n")
	print("###################################################################################")
	print("##########                    MathFeature: Concatenate                  ###########")
	print("##########                Arguments: -n number of datasets              ###########")
	print("##########                 Author: Robson Parmezan Bonidia              ###########")
	print("###################################################################################")
	print("\n")
	print("Files in the directory") 
	print(os.listdir(os.getcwd()))
	print("\n")
	#########################################################################
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--number', help='number of datasets')
	parser.add_argument('-o', '--output', help='Csv format file | E.g., combined.csv')
	args = parser.parse_args()
	n = int(args.number)
	foutput = str(args.output)
	dataset_names = []
	for i in range(1, n+1):
		name = input("Dataset %s: " % (i))
		dataset_names.append(name)
	if check_files() == 1:
		concatenate()
	else:
		print("Some file not exists")
#############################################################################
