<p align="center">
<img src="logo.svg" alt="DNA Spider-Web" title="DNASpiderWeb" width="20%"/>
</p>


# DNA SPIDER-WEB

[![CircleCI](https://circleci.com/gh/circleci/circleci-docs.svg?style=shield)](https://circleci.com/gh/HaolingZHANG/DNASpiderWeb)
[![Build Status](https://travis-ci.org/HaolingZHANG/DNASpiderWeb.svg)](https://travis-ci.org/HaolingZHANG/DNASpiderWeb)
[![Coverage Status](https://coveralls.io/repos/HaolingZHANG/DNASpiderWeb/badge.svg?branch=master&service=github)](https://coveralls.io/github/HaolingZhang/DNASpiderWeb?branch=master)
[![PythonVersion](https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9-blue)](https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9-blue)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?maxAge=259200)](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/LICENSE)


As a genetic material, DNA has become an attractive medium for storing digital information gradually.
In recent years, a growing number of functional biochemical operations and storage environments were introduced, 
biochemical constraints are not confined to long single-nucleotide repeats and abnormal GC content.
However, trade-offs between information density and compatibility to biochemical operations and algorithms require in-depth investigation, 
to create transcoding algorithms considering novel biochemical constraints and their combinations.
Here, we design an automatic generator, named SPIDER-WEB, 
to create corresponding graph-based algorithms 
which could be used directly or served as a benchmark for artificial algorithm design.
The generated coding algorithms provide 
[high code rate](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/experiments/step_1_compatibility.py),
[low parameter sensitivity](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/experiments/step_3_stability.py),
[certain self-repair ability](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/experiments/step_4_repairability.py), and
[valuable privacy protection potentiality](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/experiments/step_5_privacy.py).
Through this work, applicable algorithms with appropriate density-compatibility trade-off under arbitrary local biochemical constraints could be generated in an automated way. 
It is also suggested that more kinds of biochemical constraints can be further investigated as more complex operations would be needed in future DNA storage systems.

## Installation
You can install this package using pip:
```sh
pip install dsw
```

The packages requires a python version >=3.7, 
as well as some basic libraries listed in [requirements file](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/requirements.txt).
That is, [numpy>=1.17.1](https://pypi.org/project/numpy/) and [scipy>=1.3.1](https://pypi.org/project/scipy/).

Furthermore, if you want to try / repeat the completed experiments in this work.
Some additional libraries need to be installed, that is, 
[Chamaeleo>=1.34](https://pypi.org/project/Chamaeleo/) and [matplotlib>=3.1.1](https://pypi.org/project/matplotlib/).
These experimental Python scripts [here](https://github.com/HaolingZHANG/DNASpiderWeb/tree/main/experiments) are single threaded. 
It may take about a month to complete all experiments on a conventional laptop (reference: Intel i7-4710MQ @ 2.50GHz).
In order to further understand the experimental situation, like getting raw data (about 730MB), 
please do not hesitate to contact zhanghaoling/at/genomics.cn.

## Basic evaluation
Some investigations of 
[SPIDER-WEB](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/spiderweb.py) are recorded 
[here](https://github.com/HaolingZHANG/DNASpiderWeb/tree/main/experiments/results/figures/), 
which are generated directly through the process Python scripts 
[here](https://github.com/HaolingZHANG/DNASpiderWeb/tree/main/experiments/).
Compared with 
[HEDGES](https://www.pnas.org/content/117/31/18489.full), 
[DNA Fountain](https://www.science.org/doi/abs/10.1126/science.aaj2038) and 
[Yin-Yang Code](https://www.biorxiv.org/content/10.1101/829721v3), 
the rough (**not academic**) conclusions are as follows:

<p align="center">
<img src="./experiments/results/figures/[6-0] result overviews.svg" title="result overview" width="100%"/>
</p>

where 
* [C-R] refers to the code rate, the statistical result of which is the average number of binary bits encoded per nucleotide; 
* [C-A] refers to the codability, the statistical result of which is the success rate of encoding process;
* [P-S] refers to the parameter sensibility, the statistical result of which is the normalized standard deviation of code rate;
* [R-A] refers to the repairability, the statistical result of which is the normalized average DNA length of lossless decoding.


## Customization
### local biochemical constraints set
You can create your customized local biochemical constraint filter by inheriting [DefaultBioFilter](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/biofilter.py#L4), such as
Here is a simple example in this work named [LocalBioFilter](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/biofilter.py#L30).

### capacity approximation
Through our customized approximater, 
see [here](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/graphized.py#L537), 
you can easily obtain the capacity under the specific biochemical constraint set.
```python
from dsw import LocalBioFilter, find_vertices, connect_valid_graph, approximate_capacity

bio_filter = LocalBioFilter(max_homopolymer_runs=2, gc_range=[0.4, 0.6], undesired_motifs=["GCC"])
vertices = find_vertices(length=10, bio_filter=bio_filter)
graph = connect_valid_graph(length=10, vertices=vertices)
capacity = approximate_capacity(graph=graph)
print(capacity)
```

The capacity approximater is based on the Perron–Frobenius theorem,
we have proved its applicability for variable-length coding algorithms in the publication.

Unfortunately, when the observed length (or window length) reaches 8, 
the data capacity will achieve 32.0GB with the data type (double-precision floating-point format), 
which is unable to be allocated (both MATLAB and Python platforms).
Therefore, see [here](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/graphized.py),
we provide a series of special data structures and approximation algorithms 
to complete the representation of graphs and the solution of largest eigenvalues.
 
To verify its [reliability](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/experiments/step_2_reliability.py), 
we compared our proposed method with Numpy ``linalg.eig'' function in some small-scale matrices.
As part of the experimental results, the median value of relative error analysis is 2.97E-11, 
which represents our proposed method is high reliability and negligible error. 

<p align="center">
<img src="./experiments/results/figures/[2-2] reliability detailed.svg" title="reliability" width="100%"/>
</p>

### coding algorithm generation
According to SPIDER-WEB, you can obtain the corresponding variable-length algorithms:

```python
from dsw import LocalBioFilter, find_vertices, connect_coding_graph

bio_filter = LocalBioFilter(max_homopolymer_runs=2, gc_range=[0.4, 0.6], undesired_motifs=["ACA", "CAC", "GTG", "TGT"])
vertices = find_vertices(length=10, bio_filter=bio_filter)
graph = connect_coding_graph(length=10, vertices=vertices, threshold=1)
```
In [experiments folder](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/experiments/__init__.py), 
you can easily find 12 examples in our experiments, the local biochemical constraint set are:

| set index | homopolymer run-length constraint | regionalized GC content constraint | undesired motif constraint | capacity |
| ---- | ---- | ---- | ---- | ----|
| 1 |  2  | 50%~50% | restriction sites | 1.0000 |
| 2 |  1  |   N/A   | N/A | 1.5850 |
| 3 | N/A | 10%~30% | N/A | 1.6302 |
| 4 |  2  | 40%~60% | similar structures | 1.6698 |
| 5 |  2  | 40%~60% | N/A | 1.7761 |
| 6 | N/A | 50%~70% | N/A | 1.7958 |
| 7 |  3  | 40%~60% | N/A | 1.8114 |
| 8 |  4  | 40%~60% | N/A | 1.8152 |
| 9 |  3  |   N/A   | N/A | 1.9824 |
| 10 |  4  |   N/A   | N/A | 1.9957 |
| 11 |  5  |   N/A   | N/A | 1.9989 |
| 12 |  6  |   N/A   | N/A | 1.9997 |

where the restriction sites represent AGCT, GACGC, CAGCAG, GATATC, GGTACC, CTGCAG, GAGCTC, GTCGAC, AGTACT, ACTAGT, GCATGC, AGGCCT, and TCTAGA;
and the similar structures in Nanopore sequencing are AGA, GAG, CTC, and TCT.

There simple transcoding performances of generated coding algorithms are as shown below:

<p align="center">
<img src="./experiments/results/figures/[1-1] compatibility code rates.svg" title="compatibility" width="100%"/>
</p>

where the blue line is the approximated capacity and 
red violin plot (with median point) is tens of thousands of groups of (bit/nucleotide) results.


## Citation
If you think this repo helps or being used in your research, please consider refer this paper.
Here is a Bibtex entry:

````
@article{zhang2021spider,
  title={SPIDER-WEB is all your need for creating coding algorithms in DNA-based storage},
  author={Zhang, Haoling and Lan, Zhaojun and Ping, Zhi and Zhang, Yiwei and Shen, Yue},
  journal={XXX},
  year={XXX},
  pages={XXX},
}
````