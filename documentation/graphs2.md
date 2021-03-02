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

## Complex Networks - v2 - New Features (Faster and without threshold)

This method generates features based on complex networks, see our pipeline in this [article](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2). To use this model, follow the example below:

```sh
To run the tool (Example): $ python3.7 methods/ComplexNetworksClass-v2.py -i input -o output -l label -k kmer


Where:

-h = help

-i = Input - Fasta format file, e.g., test.fasta

-o = output - CSV format file, e.g., test.csv

-l = Label - Dataset Label, e.g., lncRNA, mRNA, sncRNA

-k = k size, e.g., 2 (k1, k2), 3 (k1, k2, k3).
```

**Running:**

```sh
$ python3.7 methods/ComplexNetworksClass-v2.py -i sequences.fasta -o sequences.csv -l mRNA -k 3
```

**Note** Input sequences for feature extraction must be in fasta format.

**Note** This example will generate a csv file with the extracted features.
