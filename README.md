<p align="center">
<img src="logo.svg" alt="DNA Spider-Web" title="DNASpiderWeb" width="40%"/>
</p>

---

# SPECTRA and DNA Spider-Web
As a genetic material, DNA has become an attractive medium for storing digital information gradually.
In recent years, a growing number of functional biochemical operations and storage environments were introduced, 
biochemical constraints are not confined to long single-nucleotide repeats and abnormal GC content.
However, trade-offs between information density and compatibility to biochemical operations and algorithms require in-depth investigation, 
to create transcoding algorithms considering novel biochemical constraints and their combinations.
Here, we construct a simple calculator, called [**SPECTRA**](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/spectra.py), 
to make use of the spectral radius to compute the density-compatibility trade-off under arbitrary local biochemical constraints.
Based on [**SPECTRA**](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/pipelines/step_3_practicability.py), 
we design an automatic generator, 
named [**SPIDER-WEB**](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/spiderweb.py), 
to create corresponding graph-based algorithms 
which could be used directly or served as a benchmark for artificial algorithm design.
Through this work, applicable algorithms with appropriate density-compatibility trade-off under arbitrary local biochemical constraints could be generated in an automated way. 
It is also suggested that more kinds of biochemical constraints can be further investigated as more complex operations would be needed in future DNA storage systems.

## Usage and Customization
### Biochemical Constraints Setting
You can create your customized local biochemical constraint filter by inheriting [DefaultBioFilter](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/biofilter.py#L4), such as

```python
from dsw.biofilter import DefaultBioFilter

class CustomizedBioFilter(DefaultBioFilter):

    def __init__(self, **requirements):
        super().__init__(screen_name="Custom")
        self.requirements = requirements
        
    def valid(self, oligo):
        # TODO do validity screening through self.requirements.
        return True
```
A simple example in this work is [LocalBioFilter](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/biofilter.py#L30).

### Directed Graph Obtained from Valid DNA Sequences
You can generate your directed graph based on the above-mentioned customized filter or local biochemical constraint group:
```python
from numpy import load
from dsw.biofilter import LocalBioFilter
from dsw.spiderweb import find_vertices, connect_default_graph

# bio_filter = CustomizedBioFilter(requirements=your requirements)
bio_filter = LocalBioFilter(max_homopolymer_runs=2, gc_range=[0.4, 0.6], undesired_motifs=["ACA", "CAC", "GTG", "TGT"])
find_vertices(length=10, bio_filter=bio_filter, save_path="local")
vertices = load(file="local" + "[v].npy")
connect_default_graph(length=10, vertices=vertices, save_path="local")
```
In [pipelines folder](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/pipelines/step_1_generate.py), 
you can easily find 8 examples.

### Trade-off Calculator
Through **SPECTRA**, you can obtain the theoretical trade-off under your requirements, that is,
```python
from numpy import load
from dsw.spectra import calculate_capacity

graph = load("local[g].npy")
print(calculate_capacity(graph=graph, replay=10))
```

### Automatic Algorithm Generator
According to **SPIDER-WEB**, you can obtain algorithms (fixed-length or variable-length):

```python
from numpy import load
from dsw.spiderweb import connect_fixed_graph, connect_variable_graph

vertices = load(file="local" + "[v].npy")
connect_fixed_graph(length=10, vertices=vertices, maximum_stride=5, save_path="../entities/FLC")
connect_variable_graph(length=10, vertices=vertices, threshold=1, save_path="../entities/VLC")
```
In [pipelines folder](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/pipelines/step_1_generate.py), 
you can easily find 8 examples.


## Basic Evaluation

Some basic evaluation results of **SPECTRA**, and **SPIDER-WEB** are recorded [here](https://github.com/HaolingZHANG/DNASpiderWeb/tree/main/results).

## Citation
If you think this repo helps or being used in your research, please consider refer this paper. 

````
@article{zhang2021general,
  title={General Trade-Off Calculator and Automatic Algorithm Generator for Arbitrary Local Biochemical Constraints in DNA Digital Data Storage},
  author={Zhang, Haoling and Lan, Zhaojun and Ping, Zhi and Zhang, Yiwei and Shen, Yue and Ge, Gennian},
  journal={XXX},
  year={XXX},
  pages={XXX},
}
````

Thank you!