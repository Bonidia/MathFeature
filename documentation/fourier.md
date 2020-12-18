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

## Numerical Mapping and Fourier Transform

To generate features based in a Fourier approach, we apply the Discrete Fourier Transform (DFT), widely used for digital image and signal processing (here GSP). To calculate DFT, we used the Fast Fourier Transform (FFT). However, to use GSP techniques, it is necessary to apply a numeric representation for the transformation or mapping of genomic data. Thereby, we study 7 numerical mapping techniques (see our pipeline in this [article](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2)). To use this model, follow the example below:

```sh
To run the code (Example): $ python3.7 methods/FourierClass.py -i input -o output -l label -r representation


Where:

-h = help

-i = Input - Fasta format file, e.g., test.fasta

-o = output - CSV format file, e.g., test.csv

-l = Label - Dataset Label, e.g., lncRNA, mRNA, sncRNA

-r = representation/mappings, e.g., 1 = Binary, 2 = Z-curve, 3 = Real, 4 = Integer, 5 = EIIP, 6 = Complex Number, 7 = Atomic Number.
```

**Running:**

```sh
$ python3.7 methods/FourierClass.py -i sequences.fasta -o sequences.csv -l mRNA -r 2
```

**Note** Input sequences for feature extraction must be in fasta format.

**Note** This example will generate a csv file with the extracted features.
