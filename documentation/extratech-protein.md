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

## Other techniques

MathFeature also provides other techniques known in the literature: k-mer (for protein), Amino acid composition (AAC), Dipeptide composition (DPC), Tripeptide composition (TPC).

**Important:** This package only accepts sequence files in *Fasta* format as input to the methods.

## Customizable k-mer, AAC, DPC, TPC

To use this model, follow the example below:

```sh
To run the code (Example): $ python3.7 methods/ExtractionTechniques-Protein.py -i input -o output -l label -t technique


Where:

-h = help

-i = Input - Fasta format file, e.g., test.fasta

-o = output - CSV format file, e.g., test.csv

-l = Label - Dataset Label, e.g., lncRNA, mRNA, sncRNA, DNA, 0, 1

-t = type of Feature Extraction - e.g., AAC or DPC or TPC or kmer or kstep
```

**Running:**

```sh
$ python3.7 methods/ExtractionTechniques-Protein.py -i protein.fasta -o dataset.csv -l DNA -t AAC
```

**Note:** Kstep =  Customizable k-mer - You can configure sliding window, step and k-mer.

**Note:** Input sequences for feature extraction must be in fasta format.

**Note:** This example will generate a csv file with the extracted features.
