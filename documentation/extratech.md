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

Before executing any method in this package, it is necessary to run a pre-processing script, to eliminate any noise from the sequences (e.g., other letters as: N, K ...). To use this script, follow the example below:

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


## Customizable k-mer, NAC, DNC, TNC

MathFeature also provides other techniques known in the literature: k-mer, Nucleic acid composition (NAC), Di-nucleotide composition (DNC), Tri-nucleotide composition (TNC). To use this model, follow the example below:

```sh
To run the code (Example): $ python3.7 methods/ExtractionTechniques.py -i input -o output -l label -t technique -seq DNA/RNA


Where:

-h = help

-i = Input - Fasta format file, e.g., test.fasta

-o = output - CSV format file, e.g., test.csv

-l = Label - Dataset Label, e.g., lncRNA, mRNA, sncRNA, DNA, 0, 1

-t = type of Feature Extraction - e.g., NAC or DNC or TNC or kmer or kstep

-seq = type of sequence, 1 = DNA and 2 = RNA'
```

**Running:**

```sh
$ python3.7 methods/ExtractionTechniques.py -i sequence.fasta -o dataset.csv -l DNA -t NAC -seq 1
```

**Note:** Kstep =  Customizable k-mer - You can configure sliding window, step and k-mer.

**Note:** Input sequences for feature extraction must be in fasta format.

**Note:** This example will generate a csv file with the extracted features.
