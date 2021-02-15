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


## Xmer k-Spaced Ymer Composition Frequency (kGap).

To use this model, follow the example below:

```sh 
To run the code (Example): $ python3.7 methods/Kgap.py -i input -o output -l label -k kgap -bef before -aft after -seq type


Where:

-i = Input - Fasta format file, E.g., test.fasta

-o = Output - CSV format file, E.g., test.csv.

-l = label - lncRNA, circRNA...

-k = gap - e.g., Frequency of kgap, E.g. 1 = A_A, 2 = A__A, 3 = A___A...

-bef = before - e.g., 1 = A_A, 2 = AA_A, 3 = AAA_A...

-aft = after - e.g., 1 = A_A, 2 = A_AA, 3 = A_AAA...

-seq = type of sequence, e.g., 1 = DNA, 2 = RNA and 3 = Protein
```

**Running:**

```sh
$ python3.7 methods/Kgap.py -i sequences.fasta -o sequences.csv -l test -k 1 -bef 1 -aft 2 -seq 1
```

**Note** Input sequences for feature extraction must be in fasta format.

**Note** This example will generate a csv file with the extracted features.
