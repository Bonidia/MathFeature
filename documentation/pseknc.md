![Python](https://img.shields.io/badge/python-v3.7-blue)
![Dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)
![Status](https://img.shields.io/badge/status-up-brightgreen)

<h1 align="center">
  <img src="https://github.com/Bonidia/MathFeature/blob/master/img/MathFeature.png" alt="MathFeature" width="350">
</h1>

<h4 align="center">Feature Extraction Package for Biological Sequences Based on Mathematical Descriptors</h4>
	
<p align="center">
  <a href="https://github.com/Bonidia/MathFeature">Home</a> •
  <a href="#authors">Key Features</a> •
  <a href="#list-of-files">List of files</a> •
  <a href="#dependencies">Dependencies</a> •
  <a href="#installing-dependencies-and-package">Installing</a> •
  <a href="#how-to-use">How To Use</a> •
  <a href="#citation">Citation</a> 
</p>

<h1 align="center"></h1>

## Preprocessing

Before executing any method in this package, it is necessary to run a pre-processing script, to eliminate any noise from the sequences (e.g., other letters as: N, K ...,). To use this script, follow the example below:

**Important:** This package only accepts sequence files in *Fasta* format as input to the methods.

```sh
To run the tool (Example): $ python3.7 preprocessing/preprocessing.py -i input -o output


Where:

-h = help

-i = Input - Fasta format file, e.g., test.fasta

-o = output - Fasta format file, e.g., output.fasta
```

**Running:**

```sh
$ python3.7 preprocessing/preprocessing.py -i dataset.fasta -o preprocessing.fasta 
```


## Pseudo k-tuple nucleotide composition - PseKNC

To use this model, follow the example below:

```sh 
To run the code (Example): $ python3.7 methods/PseKNC.py -h


Where:

-i = Input - Fasta format file, E.g., test.fasta

-o = Output - CSV format file, E.g., test.csv.

-l = label - lncRNA, circRNA...

-x = prop - e.g., Name of file containing list of properties to be used in calculations.

-xp = prop_values - e.g., Name of file containing list of properties (values) to be used in calculations.

-seq = type of sequence, e.g., 1 = DNA, 2 = RNA

-t = type, e.g., 1 - Type 1 PseKNC, 2 - Type 2 PseKNC

-k = Kind of oligonucleotide: 2 - Dinucleotide, 3 - Trinucleotide

-j = lambda, e.g., Set the value of lambda parameter in the PseKNC algorithm. Must be smaller than the length of any query sequence, E.g., 1

-w = weight e.g., Set the value of weight parameter in the PseKNC algorithm. It can be a value between (0,1], E.g, 1.0

-s = frequency, e.g., Calculate only the frequency of each oligonucleotide in the input sequence. Unless otherwise specified, E.g., 2.
```

**Running:**

```sh
Example 1 - When seq = DNA and k = 2

$ python3.7 methods/PseKNC.py -i sequence.fasta -o sequence.csv -l 1 -x files/propNames-DNA-k2.txt -xp files/propValues-DNA-k2.txt -seq 1 -t 2 -k 2 -j 1 -w 1.0 -s 2

Example 2 - When seq = DNA and k = 3

$ python3.7 methods/PseKNC.py -i sequence.fasta -o sequence.csv -l 1 -x files/propNames-DNA-k3.txt -xp files/propValues-DNA-k3.txt -seq 1 -t 2 -k 3 -j 1 -w 1.0 -s 2

Example 3 - When seq = RNA and k = 2

$ python3.7 methods/PseKNC.py -i sequence.fasta -o sequence.csv -l 1 -x files/propNames-RNA-k2.txt -xp files/propValues-RNA-k2.txt -seq 2 -t 2 -k 2 -j 1 -w 1.0 -s 2
```

**Note** Input sequences for feature extraction must be in fasta format.

**Note** This example will generate a csv file with the extracted features.
