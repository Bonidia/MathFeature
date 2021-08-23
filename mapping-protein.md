![Python](https://img.shields.io/badge/python-v3.7-blue)
![Dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)
![Status](https://img.shields.io/badge/status-up-brightgreen)

<h1 align="center">
  <img src="img/MathFeature.png" alt="MathFeature" width="350">
</h1>

<h4 align="center">Feature Extraction Package for Biological Sequences Based on Mathematical Descriptors</h4>
	
<p align="center">
  <a href="https://bonidia.github.io/MathFeature/">Home</a> •
  <a href="#authors">Key Features</a> •
  <a href="#list-of-files">List of files</a> •
  <a href="#dependencies">Dependencies</a> •
  <a href="#installing-dependencies-and-package">Installing</a> •
  <a href="#how-to-use">How To Use</a> •
  <a href="#citation">Citation</a> 
</p>

<h1 align="center"></h1>

## Numerical Mapping - Protein

This method generates a numerical mapping of all sequence. Essentially, we provide 4 mappings (Protein). Nevertheless, this method will generate a vector with the size of the largest sequence. We developed a code that applies everything automatically. Therefore, it is necessary to pass all the classes/labels that will form the dataset. Thereby. to use this model, follow the example below:

```sh 
To run the code (Example): $ python3.7 methods/Mappings-Protein.py -n number of datasets/labels -o output -r representation


Where:

-h = help

-n = number of datasets/labels

-o = output - CSV format file, e.g., test.csv

-r = representation/mappings, e.g., 1 = Accumulated Amino Acid Frequency, 3 = Kmer Frequency Mapping, 5 = Integer Mapping, 7 = EIIP Mapping
```

**Running:**

```sh
$ python3.7 methods/Mappings-Protein.py -n 2 -o dataset.csv -r 2
```

**Note** Input sequences for feature extraction must be in fasta format.

**Note** This example will generate a csv file with the extracted features.