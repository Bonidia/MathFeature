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
  <a href="https://bonidia.github.io/MathFeature/">Documentation</a> •
  <a href="#list-of-files">List of files</a> •
  <a href="#dependencies">Dependencies</a> •
  <a href="#installing-dependencies-and-package">Installing</a> •
  <a href="#list-of-descriptors">Descriptors</a> •
  <a href="#how-to-use">How To Use</a> •
  <a href="#citation">Citation</a> 
</p>

<h1 align="center"></h1>

Machine learning algorithms have been very successfully applied to extract new and relevant knowledge from biological sequences. However, the predictive performance of these algorithms is largely affected by how the sequences are represented. Thereby, the main challenge is how to numerically represent a biological sequence in a numeric vector with an efficient mathematical expression. Several feature extraction techniques have been proposed for biological sequences, where most of them are available in feature extraction packages. However, there are relevant approaches that are not available in existing packages, techniques based on mathematical descriptors, e.g., Fourier, entropy, and graphs. Therefore, this paper presents a new package, named MathFeature, which implements mathematical descriptors able to extract relevant information from biological sequences. MathFeature provides 20 approaches based on several studies found in the literature, e.g., multiple numeric mappings, genomic signal processing, chaos game theory, entropy, and complex networks. MathFeature also allows the extraction of alternative features, complementing the existing packages.


## Authors

* Robson Parmezan Bonidia, Danilo Sipoli Sanches, and André Carlos Ponce de Leon Ferreira de Carvalho.

* **Correspondence:** rpbonidia@gmail.com or bonidia@usp.br


## Publication

Submitted


## List of files

 - **examples:** Files of Example;
 - **methods:** Main Files - Feature Extraction Models, e.g., Fourier, Numerical Mapping, Entropy, Complex Networks;
 - **preprocessing:** Preprocessing Files;
 - **README:** Documentation;
 - **requirements:** Dependencies.


## Dependencies

- Python (>=3.7.3)
- Biopython
- Igraph
- NumPy 
- Pandas
- SciPy


## Installing dependencies and package

It is important to note that we consider that the Python language is installed. Otherwise, access [here](https://www.python.org/downloads/release/python-375/).

```sh
$ git clone https://github.com/Bonidia/MathFeature.git MathFeature

$ cd MathFeature

$ pip3 install -r requirements.txt

$ apt-get -y install python3-igraph
```

## List of Descriptors

Descriptors calculated by MathFeature for DNA, RNA, and Protein sequences: [Click here.](documentation/descriptors.md)

## How to use

We proposed an open-source Python package called MathFeature, that implements feature extraction approaches using mathematical features, including 20 descriptors organized into five categories. To our best knowledge, MathFeature is the first package that computes biological sequence features based on various mathematical descriptors. In this section, 5 feature extraction groups are available: **(1)** numerical mapping techniques, **(2)** numerical mapping techniques with Fourier transform, **(3)** techniques with game chaos, **(4)** techniques with Entropy, **(5)** techniques with complex networks. Moreover, we provide some additional scripts for feature extraction and preprocessing.

See our [documentation](https://bonidia.github.io/MathFeature).

## Feature Selection

If you want to apply feature selection techniques, visit our [repository](https://github.com/Bonidia/FeatureSelection-FSRV).

## Citation

If you use this code in a scientific publication, we would appreciate citations to the following paper:

Submitted.
