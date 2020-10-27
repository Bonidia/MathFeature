![Python](https://img.shields.io/badge/python-v3.7-blue)
![Dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

<h1 align="center">
  <img src="https://github.com/Bonidia/MathFeature/blob/master/img/MathFeature.png" alt="MathFeature" width="350">
</h1>

<h4 align="center">Feature Extraction Package for Biological Sequences Based on Mathematical Approaches</h4>

<p align="center">
  <a href="https://github.com/Bonidia/MathFeature">Home</a> •
  <a href="#authors">Key Features</a> •
  <a href="#list-of-files">List of files</a> •
  <a href="#dependencies">Dependencies</a> •
  <a href="#installing-dependencies-and-package">Installing</a> •
  <a href="#how-to-use">How To Use</a> •
  <a href="#license">License</a> •
  <a href="#citation">Citation</a> 
</p>

<h1 align="center"></h1>

Lorem ipsum dolor sit amet, consectetur adipiscing elit. In vitae lorem luctus, mattis ligula eget, dignissim ex. Nunc rutrum, mauris sit amet lacinia pellentesque, felis sem egestas nunc, sit amet mattis lacus velit vitae purus. Sed scelerisque bibendum metus, non aliquam neque dictum non. Aenean iaculis lobortis tempor. Donec pretium sem accumsan tellus tincidunt commodo. Proin egestas ante ligula, eget luctus erat scelerisque quis. Aenean venenatis gravida massa, id elementum lectus tempus quis. Duis maximus odio sit amet quam commodo sodales. 


## Authors

* Robson Parmezan Bonidia, André Carlos Ponce de Leon Ferreira de Carvalho, and Danilo Sipoli Sanches.

* **Correspondence:** robservidor@gmail.com


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

It is important to note that we consider that the Python language is installed. Otherwise, access: https://www.python.org/downloads/release/python-375/.

```sh
$ git clone https://github.com/Bonidia/MathFeature.git MathFeature

$ cd MathFeature

$ pip3 install -r requirements.txt

$ apt-get -y install python3-igraph
```

## How to use

We proposed an open-source Python package called MathFeature, that implements feature extraction approaches using mathematical models, including 20 approaches/descriptors organized into five categories. To our best knowledge, MathFeature is the first package that computes biological sequence features based on various mathematical approaches. In this section, 5 feature extraction approaches are available: **7** numerical mapping techniques, **7** numerical mapping techniques with Fourier transform, **3** techniques with game chaos, **2** techniques with Entropy, **1** with Complex Networks. Moreover, we provide some scripts for preprocessing.

* [Preprocessing](https://github.com/Bonidia/MathFeature/blob/master/documentation/preprocessing.md)
* [Mapping Numerical](https://github.com/Bonidia/MathFeature/blob/master/documentation/mapping.md)
* [Chaos Game Representation](https://github.com/Bonidia/MathFeature/blob/master/documentation/chaos.md)
* [Mapping Numerical and Fourier Transform](https://github.com/Bonidia/MathFeature/blob/master/documentation/fourier.md)
* [Shannon and Tsallis Entropy](https://github.com/Bonidia/MathFeature/blob/master/documentation/entropy.md)
* [Complex Networks](https://github.com/Bonidia/MathFeature/blob/master/documentation/graphs.md)
* [Only k-mer](https://github.com/Bonidia/MathFeature/blob/master/documentation/kmer.md)


## License

This project is licensed under the MIT License - see the LICENSE.md file for details


## Citation

If you use this code in a scientific publication, we would appreciate citations to the following paper:

Submitted.
