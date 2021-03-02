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

# List of Descriptors

Descriptors calculated by MathFeature for DNA, RNA, and Protein sequences.


| Descriptor groups     | Descriptor                           | Dimension  | Sequence        | Example (Study with Application or Theory)                                     |
|    :---:              |  :---:                               |  :---:     |  :---:          |  :---:                                                                         |
|                                                                                                                                                                              |
|                       | Binary                               | L * 4      | DNA/RNA         | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2) - [Ref 2](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.68.3805) |
|                       | Z-curve                              | L * 3      | DNA/RNA         | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2) - [Ref 2](https://www.tandfonline.com/doi/abs/10.1080/07391102.1994.10508031)|
|                       | Real                                 | L          | DNA/RNA         | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2) - [Ref 2](https://link.springer.com/article/10.1155/S111086570430925X)       |
| **Numerical Mapping** | Integer                              | L          | DNA/RNA/Protein      | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2) - [Ref 2](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1582-4934.2002.tb00196.x) |
|                       | EIIP                                 | L          | DNA/RNA/Protein      | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2) - [Ref 2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1891688/) |
|                       | Complex Number                       | L          | DNA/RNA         | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2) - [Ref 2](https://ieeexplore.ieee.org/abstract/document/8361572)     |
|                       | Atomic Number                        | L          | DNA/RNA         | [Ref 1](https://doi.org/10.1117/12.732283) - [Ref 2]( https://doi.org/10.1371/journal.pone.0173288)                                       |
|                                                                                                                                                                              |
|                       | Chaos Game Representation            | L * 2      | DNA/RNA         | [Ref 1](https://doi.org/10.1093/nar/18.8.2163)                                                                                            |                                                                             |
| **Chaos Game**        | Frequency Chaos Game Representation  | L - k + 1  | DNA/RNA         |                                                                                |
|                       | Chaos Game Signal (with Fourier)     | 19         | DNA/RNA         | [Ref 1](https://doi.org/10.1016/j.ygeno.2016.08.002)                                                                                      |                                                                                                                                                                              
| **Fourier Transform** | Numerical Mapping + Fourier          | 19         | DNA/RNA/Protein | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2)           |
|                                                                                                                                                                              |
| **Entropy**           | Shannon                              | k          | DNA/RNA/Protein | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2)           |
|                       | Tsallis                              | k          | DNA/RNA/Protein | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2)           |
|                                                                                                                                                                              |
| **Graphs**            | Complex Networks (with threshold)         | 12 * t     | DNA/RNA/Protein | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2) - [Ref 2](https://doi.org/10.1093/nar/gky462) |
|                       | Complex Networks (without threshold - v2) | 27 * k | DNA/RNA/Protein | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2) - [Ref 2](https://doi.org/10.1093/nar/gky462) |
|                                                                                                                                                                              |
|                       | Basic k-mer                          | 4^k        | DNA/RNA         | [Ref 1](https://www.biorxiv.org/content/10.1101/2020.06.08.140368v2)   | 
|                       | Customizable k-mer                   | 4^k        | DNA/RNA         |                                                                                
|                       | Nucleic acid composition (NAC)       | 4          | DNA/RNA         | [Ref 1](https://doi.org/10.1016/j.cmpb.2017.05.008)                            
|                       | Di-nucleotide composition (DNC)      | 16         | DNA/RNA         | [Ref 1](https://doi.org/10.1016/j.cmpb.2017.05.008)                            
|                       | Tri-nucleotide composition (TNC)     | 64         | DNA/RNA         | [Ref 1](https://doi.org/10.1016/j.cmpb.2017.05.008) 
|                       | ORF Features or Coding Features      | 10         | DNA/RNA         | [Ref 1](https://doi.org/10.1093/bib/bbab011) - [Ref 2](https://www.nature.com/articles/srep34838) |
|                       | Fickett score                        | 2          | DNA/RNA         | [Ref 1](https://academic.oup.com/nar/article/41/6/e74/2902455)
|                       | Pseudo K-tuple nucleotide composition | -          | DNA/RNA         | [Ref 1](https://doi.org/10.1016/j.ab.2014.04.001)
| **Other techniques**  | Accumulated Nucleotide Frequency-ANF | L          | DNA/RNA/Protein | [Ref 1](https://www.nature.com/articles/srep13859)                           
|                       | ANF with Fourier                     | 19         | DNA/RNA/Protein |   
| 			| Xmer k-Spaced Ymer Composition Frequency (kGap)   | 4^X * 4^Y or 20^X * 20^Y  | DNA/RNA/Protein | [Ref 1](https://doi.org/10.3389/fbioe.2020.00134) - [Ref 2](https://doi.org/10.1093/bioinformatics/btz165) |  
|                       | Amino acid composition (AAC)         | 20         | Protein         | [Ref 1](https://doi.org/10.3389/fcell.2020.578901)                             
|                       | Dipeptide composition (DPC)          | 400        | Protein         | [Ref 1](https://doi.org/10.3389/fcell.2020.578901)                             
|                       | Tripeptide composition (TPC)         | 8000       | Protein         | [Ref 1](https://doi.org/10.3389/fcell.2020.578901)                            
|                       | Basic k-mer                          | 20^k       | Protein         |                                                                                
|                       | Customizable k-mer                   | 20^k       | Protein         |   
|                       | Kmer Frequency Mapping               | L - k + 1  | Protein         |
|                       | Kmer Frequency Mapping with Fourier  | 19         | Protein         |



To use any descriptor, see our [documentation](https://bonidia.github.io/MathFeature/).

**Note 1:** L = length of the longest sequence.

**Note 2:** k = frequencies of k-mer.

**Note 3:** t = threshold: number of subgraphs.

**Note 4:** The reference column represents some studies that apply the descriptor (Similar approach). Other references are cited in our article.
