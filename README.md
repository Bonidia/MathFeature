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
  <a href="http://mathfeature.icmc.usp.br/">Web Server</a> •
  <a href="https://github.com/Bonidia/MathFeature-WebServer">Web Server - LocalHost</a> •
  <a href="#list-of-files">List of files</a> •
  <a href="#dependencies">Dependencies</a> •
  <a href="#installing-dependencies-and-package">Installing</a> •
  <a href="#list-of-descriptors">Descriptors</a> •
  <a href="#how-to-use">How To Use</a> •
  <a href="#GUI">GUI</a> •
  <a href="#citation">Citation</a> 
</p>

<h1 align="center"></h1>

## Update news!!!

**MathFeature** - Now available on Web Server - Access on http://mathfeature.icmc.usp.br/ - 2021-09-30

**Web Server** - To work on your machine or network - Access on https://github.com/Bonidia/MathFeature-WebServer - 2021-11-05

## Abstract

One of the main challenges in the application of Machine Learning (ML) algorithms to biological sequence data is how to numerically represent a sequence in a numeric input vector. Feature extraction techniques capable of extracting numerical information from biological sequences have been reported in the literature. However, many of these techniques are not available in existing packages, such as mathematical descriptors. This paper presents a new package, MathFeature, which implements mathematical descriptors able to extract relevant numerical information from biological sequences, i.e., DNA, RNA, and Proteins (prediction of structural features along the primary sequence of amino acids). MathFeature makes available 20 numerical feature extraction descriptors based on approaches found in the literature, e.g., multiple numeric mappings, genomic signal processing, chaos game theory, entropy, and complex networks. MathFeature also allows the extraction of alternative features, complementing the existing packages. To ensure that our descriptors are robust and to assess their relevance, experimental results are presented in nine case studies. According to these results, the features extracted by MathFeature shown high performance (0.6350-0.9897, accuracy), both applying only mathematical descriptors, but also hybridization with well-known descriptors in the literature. Finally, through MathFeature, we overcome several studies in eight benchmark datasets, exemplifying the robustness and viability of the proposed package. MathFeature advances in the area by bringing descriptors not available in other packages, as well as allowing non-experts to use feature extraction techniques.

## Authors

* Robson Parmezan Bonidia, Douglas S. Domingues, Danilo Sipoli Sanches, and André Carlos Ponce de Leon Ferreira de Carvalho.

* **Correspondence:** rpbonidia@gmail.com or bonidia@usp.br


## Publication

Robson P Bonidia, Douglas S Domingues, Danilo S Sanches, André C P L F de Carvalho, MathFeature: feature extraction package for DNA, RNA and protein sequences based on mathematical descriptors, Briefings in Bioinformatics, 2021; bbab434, https://doi.org/10.1093/bib/bbab434.

## List of files

 - **case studies:** case studies used in our article;
 - **GUI:** GUI (Graphical User Interface)-based platform;
 - **examples:** Files of Example;
 - **files:** files used in some methods;
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

MathFeature can be run on the console, but we also provide a GUI-based platform.

## Docker Image - Terminal - MathFeature v1.0

It is important to note that we consider that the Docker is installed.

Docker commands - Examples

1 - https://www.docker.com/sites/default/files/d8/2019-09/docker-cheat-sheet.pdf

2 - https://dockerlabs.collabnix.com/docker/cheatsheet/

3 - https://github.com/wsargent/docker-cheat-sheet

```sh
$ docker pull bio21061993/mathfeature:latest

$ docker run -it --name mathfeature-terminal bio21061993/mathfeature bash

$ git clone https://github.com/Bonidia/MathFeature.git MathFeature-Terminal

$ cd MathFeature-Terminal

$ conda activate mathfeature-terminal

$ Run the desired scripts - See Documentation
```

## Docker Image - GUI - MathFeature - Option 1

Docker commands - Examples - Running GUI applications using Docker 

1 - https://cuneyt.aliustaoglu.biz/en/running-gui-applications-in-docker-on-windows-linux-mac-hosts/

2 - https://betterprogramming.pub/running-desktop-apps-in-docker-43a70a5265c4

3 - https://medium.com/@SaravSun/running-gui-applications-inside-docker-containers-83d65c0db110

4 - https://sourabhbajaj.com/blog/2017/02/07/gui-applications-docker-mac/

5 - https://www.cloudsavvyit.com/10520/how-to-run-gui-applications-in-a-docker-container/

```sh
$ docker pull bio21061993/mathfeature-gui

$ For Linux: sudo docker run -it --name mathfeature-gui --net=host --env="DISPLAY" -v /root:/root --volume="$HOME/.Xauthority:/root/.Xauthority:rw" bio21061993/mathfeature-gui

$ For Mac: docker run -it --name mathfeature-gui -e DISPLAY=$IP:0 -v /tmp/.X11-unix:/tmp/.X11-unix bio21061993/mathfeature-gui

$ git clone https://github.com/Bonidia/MathFeature.git MathFeature

$ cd MathFeature/

$ export PATH=/home/miniconda3/bin:$PATH

$ conda activate mathfeature-gui or source activate mathfeature-gui

$ python GUI/main.py

$ Run the desired scripts - See Documentation
```

## Terminal - Option 2

It is important to note that we consider that the Python language is installed. Otherwise, access [here](https://www.python.org/downloads/release/python-375/).

```sh
$ git clone https://github.com/Bonidia/MathFeature.git MathFeature

$ cd MathFeature

$ pip3 install -r requirements.txt

$ apt-get -y install python3-igraph
```

## Conda - Terminal - Option 3

Another way to install MathFeature is by using miniconda, e.g.:

```sh
$ git clone https://github.com/Bonidia/MathFeature.git MathFeature

$ cd MathFeature
```

**1 - Install Miniconda:** 

```sh

See documentation: https://docs.conda.io/en/latest/miniconda.html

$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

$ chmod +x Miniconda3-latest-Linux-x86_64.sh

$ ./Miniconda3-latest-Linux-x86_64.sh

$ export PATH=~/miniconda3/bin:$PATH

```

**2 - Create environment:**

```sh

conda env create -f mathfeature-terminal.yml -n mathfeature-terminal

```

**3 - Activate environment:**

```sh

conda activate mathfeature-terminal

```

**4 - You can deactivate the environment, using:**

```sh

conda deactivate

```

## Conda - GUI - Option 4 - Linux

Another way to install MathFeature is by using GUI, e.g.:

```sh
$ git clone https://github.com/Bonidia/MathFeature.git MathFeature

$ cd MathFeature
```

**1 - Install Miniconda:** 

```sh

See documentation: https://docs.conda.io/en/latest/miniconda.html

$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

$ chmod +x Miniconda3-latest-Linux-x86_64.sh

$ ./Miniconda3-latest-Linux-x86_64.sh

$ export PATH=~/miniconda3/bin:$PATH

$ apt-get install libgtk-3-dev

```

**2 - Create environment:**

```sh

conda env create -f mathfeature-gui.yml -n mathfeature-gui

```

**3 - Activate environment:**

```sh

conda activate mathfeature-gui

python GUI/main.py

```

**4 - You can deactivate the environment, using:**

```sh

conda deactivate

```

## List of Descriptors

Descriptors calculated by MathFeature for DNA, RNA, and Protein sequences: [Click here.](documentation/descriptors.md)

## How to use

We proposed an open-source Python package called MathFeature, that implements feature extraction approaches using mathematical features, including 20 descriptors organized into five categories. To our best knowledge, MathFeature is the first package that computes biological sequence features based on various mathematical descriptors. In this section, 5 feature extraction groups are available: **(1)** numerical mapping techniques, **(2)** numerical mapping techniques with Fourier transform, **(3)** techniques with game chaos, **(4)** techniques with Entropy, **(5)** techniques with complex networks. Moreover, we provide some additional scripts for feature extraction and preprocessing.

See our [documentation](https://bonidia.github.io/MathFeature).

## GUI

<h1 align="left">
  <img src="https://github.com/Bonidia/MathFeature/blob/master/img/math1.png" alt="MathFeature" width="650">
</h1>

<h1 align="left">
  <img src="https://github.com/Bonidia/MathFeature/blob/master/img/math2.png" alt="MathFeature" width="650">
</h1>


## Feature Selection

If you want to apply feature selection techniques, visit our [repository](https://github.com/Bonidia/FeatureSelection-FSRV).

## Citation

If you use this code in a scientific publication, we would appreciate citations to the following paper:

Robson P Bonidia, Douglas S Domingues, Danilo S Sanches, André C P L F de Carvalho, MathFeature: feature extraction package for DNA, RNA and protein sequences based on mathematical descriptors, Briefings in Bioinformatics, 2021; bbab434, https://doi.org/10.1093/bib/bbab434.

```sh

@article{10.1093/bib/bbab434,
    author = {Bonidia, Robson P and Domingues, Douglas S and Sanches, Danilo S and de Carvalho, André C P L F},
    title = "{MathFeature: feature extraction package for DNA, RNA and protein sequences based on mathematical descriptors}",
    journal = {Briefings in Bioinformatics},
    year = {2021},
    month = {11},
    issn = {1477-4054},
    doi = {10.1093/bib/bbab434},
    url = {https://doi.org/10.1093/bib/bbab434},
    note = {bbab434},
    eprint = {https://academic.oup.com/bib/advance-article-pdf/doi/10.1093/bib/bbab434/41108442/bbab434.pdf},
}

```
