<p align="center">
<img src="./docs/source/_static/logo.svg" alt="DNA Spider-Web" title="DNASpiderWeb" width="100%"/>
</p>

[![CircleCI](https://circleci.com/gh/HaolingZHANG/DNASpiderWeb/tree/main.svg?style=shield&circle-token=6aeac22720c5828585591fa5e2f4917bcaae9a72)](https://circleci.com/gh/HaolingZHANG/DNASpiderWeb/tree/main)
[![PythonVersion](https://img.shields.io/badge/python-3.7-blue)](https://img.shields.io/badge/python-3.7-blue)
[![License](https://img.shields.io/badge/License-BGI_Research-orange.svg)](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/LICENSE.pdf)


DNA has been considered a promising medium for storing digital information.
As an essential step in the DNA-based data storage workflow, 
coding scheme is responsible to implement functions including bit-to-base transcoding, error correction, etc. 
In previous studies, these functions are normally realized by introducing independent coding algorithms. 
Here, we report a graph-based architecture, named SPIDER-WEB, 
providing an all-in-one coding solution by generating customized algorithms automatically. 
SPIDER-WEB supports correcting at maximum 4% edit errors in the DNA sequences 
including substitution and insertion/deletion (indel), 
while the required logical redundancy at only 5.5%. 
Since no DNA sequence pretreatment is required for correcting and decoding processes, 
SPIDER-WEB offers the function of real-time information retrieval, 
which is 305.08 times faster than the speed of single-molecule sequencing techniques. 
Our retrieval process can improve 2 orders of magnitude faster compared to conventional one 
under megabyte-level data and can be scalable to fit exabyte-level data. 
Therefore, SPIDER-WEB holds the potential to improve the practicability in large-scale data storage applications.

## Installation
You can install this package using pip:
```sh
pip install DNASpiderWeb
```

The packages requires a python version >=3.7, 
as well as some basic libraries 
(only [numpy 1.17.1](https://pypi.org/project/numpy/) and [networkx 2.6.3](https://pypi.org/project/networkx/).
The license is customized by the BGI-Research, see 
[here](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/LICENSE.pdf).

Furthermore, if you want to try / repeat the completed experiments in this work.
Some additional libraries need to be installed, that is, 
[Chamaeleo 1.34](https://pypi.org/project/Chamaeleo/), [matplotlib 3.1.1](https://pypi.org/project/matplotlib/), 
and [biopython 1.78](https://pypi.org/project/biopython/).
These experimental Python scripts 
[here](https://github.com/HaolingZHANG/DNASpiderWeb/tree/main/experiments) are single threaded. 
It may take about several months to complete all experiments on a conventional laptop 
(reference: Intel i7-4710MQ @ 2.50GHz).
In order to further understand the experimental situation, the core raw data are saved
[here](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/experiments/raw).
If you want the whole raw data, please do not hesitate to contact us.

In addition, the module usage and customization information are shown on the 
[ReadtheDocs website](https://dnaspiderweb.readthedocs.io/en/latest/index.html).

If you are interested in detailed design, evaluations, conclusions, mathematical proofs, and implementations, 
please refer to our [publication](https://arxiv.org/abs/2204.02855v3).

## Repository Structure
The structure of this library is shown below:
```html
├── dsw                                     // Source codes of SPIDER-WEB.
│    ├── __init__.py                        // Exhibition of class and method calls.
│    ├── biofilter.py                       // Biochemical constraint filter to judge whether the candidate DNA string is valid or invalid.
│    │    ├── DefaultBioFilter              // Default biochemical constraint filter inherited by all related filters.
│    │    ├── LocalBioFilter                // Local biochemical constraint filter in our work.
│    ├── graphized.py                       // Special data structures and functions related to graph theory.
│    │    ├── get_complete_accessor         // Get a complete accessor with the required observed length.
│    │    ├── adjacency_matrix_to_accessor  // Convert the adjacency matrix to the equivalent accessor (compressed matrix).
│    │    ├── accessor_to_adjacency_matrix  // Convert the accessor to its equivalent adjacency matrix.
│    │    ├── latter_map_to_accessor        // Convert the latter map (linked storage structure of graph) to its equivalent accessor.
│    │    ├── accessor_to_latter_map        // Convert the accessor to its equivalent latter map.
│    │    ├── remove_useless                // Remove useless vertices (the out-degree of witch less than threshold) in the latter map.
│    │    ├── obtain_formers                // Obtain in-degree vertex indices based on the current vertex index.
│    │    ├── obtain_latters                // Obtain out-degree vertex indices based on the current vertex index.
│    │    ├── obtain_leaf_vertices          // Obtain leaf vertex indices based on the current vertex index and the depth.
│    │    ├── approximate_capacity          // Approximate the capacity of the specific graph through Perron–Frobenius theorem.
│    │    ├── path_matching                 // Perform saturation repair by matching the path of the accessor.
│    │    ├── calculate_intersection_score  // Calculate the intersection score based on the breach-first search (further version).
│    ├── operation.py                       // Progress monitor and digital calculation operation.
│    │    ├── Monitor                       // Monitor which outputting the progress based on current state and total state.
│    │    ├── calculus_addition             // Do huge number addition calculus with a small base value, as number + base.
│    │    ├── calculus_subtraction          // Do huge number subtraction calculus with a small base value, as number - base.
│    │    ├── calculus_multiplication       // Do huge number multiplication calculus with a small base value, as number * base.
│    │    ├── calculus_division             // Do huge number division calculus with a small base value, as number / base and number % base.
│    │    ├── bit_to_number                 // Convert a bit array to its equivalent decimal number.
│    │    ├── number_to_bit                 // Convert a decimal number to its equivalent bit array with specific length.
│    │    ├── dna_to_number                 // Convert a DNA string to its equivalent decimal number.
│    │    ├── number_to_dna                 // Convert a decimal number to its equivalent DNA string with specific length.
│    ├── spiderweb.py                       // Generating, transcoding, repairing pipelines of SPIDER-WEB.
│    │    ├── encode                        // Encode a bit array by the specific accessor.
│    │    ├── decode                        // Decode a DNA string by the specific accessor.
│    │    ├── set_vt                        // Set (or calculate) Varshamov-Tenengolts-based path check for DNA string.
│    │    ├── repair_dna                    // Repair the DNA string containing one or more errors.
│    │    ├── find_vertices                 // Find valid vertices based on the given the biochemical constraints.
│    │    ├── connect_valid_graph           // Connect a valid graph by valid vertices.
│    │    ├── connect_coding_graph          // Connect a coding algorithm by valid vertices and the threshold for minimum out-degree.
│    │    ├── remove_nasty_arc              // Remove the nasty arc based on the intersection scores (further version).
│    │    ├── create_random_shuffles        // Create the shuffles for accessor through the random mechanism.
├── experiments                             // Experiment module of SPIDER-WEB.
│    ├── __init__.py                        // Preset parameters in the simulation experiment.
│    ├── code_encode.py                     // Script in the encoding simulation process.
│    ├── code_repair.py                     // Script in the correcting simulation process.
│    ├── evaluations.py                     // Script of all the evaluation experiments.
│    ├── show_main.py                       // Script showing data in the main text.
│    ├── show_supp.py                       // Script showing data in the supplementary.
│    ├── sort_data.py                       // Script arranging the core raw data.
├── tests                                   // Test module of source codes.
│    ├── test_accessor_vs_latter_map.py     // Unit test for the conversion between the accessor and the latter map.
│    ├── test_accessor_vs_matrix.py         // Unit test for the conversion between the accessor and the adjacency matrix.
│    ├── test_bio_filters.py                // Unit test for the correctness of the biochemical constraint filter.
│    ├── test_capacities.py                 // Unit test for the reliability if the capacity approximation.
│    ├── test_coding.py                     // Unit test for the default or faster encoding/decoding correctness.
│    ├── test_generating.py                 // Unit test for the generating correctness.
│    ├── test_number_vs_binary_message.py   // Unit test for the conversion between the decimal number and binary message.
│    ├── test_number_vs_dna_string.py       // Unit test for the conversion between the decimal number and DNA string.
│    ├── test_operations.py                 // Unit test for the correctness of large number basic operations.
│    ├── test_repair.py                     // Unit test for the correcting process.
│    ├── test_shuffles.py                   // Unit test for the encoding/decoding correctness when using the shuffle strategy.
├── README.md                               // Description document of library.
```
The installation process only includes folder 'dsw' and  'tests'.


## Citation
If you think this repository helps or being used in your research, please consider refer this work.
Here is a Bibtex entry:

````
@article{zhang2021spider,
  title={SPIDER-WEB generates coding algorithms with superior error tolerance and real-time information retrieval capacity},
  author={Zhang, Haoling and Lan, Zhaojun and Zhang, Wenwei and Xu, Xun and Ping, Zhi and Zhang, Yiwei and Shen, Yue},
  journal={arXiv preprint arXiv:2204.02855},
  year={2022}
}
````