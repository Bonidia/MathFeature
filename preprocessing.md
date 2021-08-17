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



## Script: Concatenate datasets/descriptors

Concatenate datasets/descriptors. To use this script, follow the example below:

**Important:** *CSV* format as input to the method.

```sh
To run the tool (Example): $ python3.7 preprocessing/concatenate.py -n number of datasets/files -o output


Where:

-h = help

-n = Input - number of datasets/files, e.g., orf.csv, entropy.csv

-o = output - Fasta format file, e.g., output.csv
```

**Running:**

```sh
$ python3.7 preprocessing/concatenate.py -n 2 -o hybrid_dataset.csv 
```


## Script: Preprocessing

To eliminate any noise from the sequences (e.g., other letters as: N, K ...,). To use this script, follow the example below:

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
$ python3.7 preprocessing/preprocessing.py -i sequences.fasta -o preprocessing.fasta 
```


## Script: Sampling

Sample sequences, select a quantity. To use this script, follow the example below:

**Important:** This package only accepts sequence files in *Fasta* format as input to the methods.

```sh
To run the tool (Example): $ python3.7 preprocessing/sampling.py -i input -o output -p samples


Where:

-h = help

-i = Input - Fasta format file, e.g., test.fasta

-o = output - Fasta format file, e.g., output.fasta

-p = samples - Number of samples
```

**Running:**

```sh
$ python3.7 preprocessing/sampling.py -i sequences.fasta -o sampling.fasta -p 1000
```

## Script: Sequences Count

Returns the number of sequences/samples. To use this script, follow the example below:

**Important:** This package only accepts sequence files in *Fasta* format as input to the methods.

```sh
To run the tool (Example): $ python3.7 preprocessing/count_sequences.py -i input


Where:

-h = help

-i = Input - Fasta format file, e.g., test.fasta
```

**Running:**

```sh
$ python3.7 preprocessing/count_sequences.py -i sequences.fasta
```

## Script: Bases/Nucleotide Count by Sequences

Bases/Nucleotide Count by Sequence. To use this script, follow the example below:

**Important:** This package only accepts sequence files in *Fasta* format as input to the methods.

```sh
To run the tool (Example): $ python3.7 preprocessing/count_bases.py -i input


Where:

-h = help

-i = Input - Fasta format file, e.g., test.fasta
```

**Running:**

```sh
$ python3.7 preprocessing/count_bases.py -i sequences.fasta
```

## Script: Redundancy 1

Removes repeated sequences in a file. To use this script, follow the example below:

**Important:** This package only accepts sequence files in *Fasta* format as input to the methods.

```sh
To run the tool (Example): $ python3.7 preprocessing/remove_redundancy.py -i input


Where:

-h = help

-i = Input - Fasta format file, e.g., test.fasta
```

**Running:**

```sh
$ python3.7 preprocessing/remove_redundancy.py -i sequences.fasta
```

## Script: Redundancy 2

Removes repeated sequences between two files. To use this script, follow the example below:

**Important:** This package only accepts sequence files in *Fasta* format as input to the methods.

```sh
To run the tool (Example): $ python3.7 preprocessing/remove_equal_sequences.py -i input1 -t input2


Where:

-h = help

-i = Input 1 - Fasta format file, e.g., test.fasta

-t = Input 2 - Fasta format file, e.g., test.fasta
```

**Running:**

```sh
$ python3.7 preprocessing/remove_equal_sequences.py -i sequences1.fasta -t sequences2.fasta
```
