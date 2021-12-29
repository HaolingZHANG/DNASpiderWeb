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
### local biochemical constraints set
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

### Directed graph generated from valid DNA string
You can generate your directed graph based on the above-mentioned customized filter or local biochemical constraint group:
```python
from dsw import LocalBioFilter, find_vertices, connect_valid_graph

# bio_filter = CustomizedBioFilter(requirements=your requirements)
bio_filter = LocalBioFilter(max_homopolymer_runs=2, gc_range=[0.4, 0.6], undesired_motifs=["ACA", "CAC", "GTG", "TGT"])
vertices = find_vertices(length=10, bio_filter=bio_filter, save_path="local")
graph = connect_valid_graph(length=10, vertices=vertices, save_path="local")
```
Through our customized approximater, 
see [here](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/graphoids.py#L476), 
you can easily obtain the capacity under the specific biochemical constraint set.
```python
from dsw import LocalBioFilter, find_vertices, connect_valid_graph, approximate_capacity

bio_filter = LocalBioFilter(max_homopolymer_runs=2, gc_range=[0.4, 0.6], undesired_motifs=["ACA", "CAC", "GTG", "TGT"])
vertices = find_vertices(length=10, bio_filter=bio_filter, save_path="local")
graph = connect_valid_graph(length=10, vertices=vertices, save_path="local")
print(approximate_capacity(graph=graph, replay=10))
```


### Automatic algorithm generator named SPIDER-WEB
According to **SPIDER-WEB**, you can obtain the corresponding variable-length algorithms:

```python
from numpy import load
from dsw import connect_coding_graph

vertices = load(file="local" + "[v].npy")
connect_coding_graph(length=10, vertices=vertices, threshold=1, save_path="../results/data/")
```
In [experiments folder](https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/experiments/__init__.py), 
you can easily find 12 examples in our experiments.

## Basic evaluations

Some basic evaluation results of **SPIDER-WEB** are recorded 
[here](https://github.com/HaolingZHANG/DNASpiderWeb/tree/main/experiments/results/figures/).

These figures are generated directly through the process Python scripts 
[here](https://github.com/HaolingZHANG/DNASpiderWeb/tree/main/experiments/).

Compared with 
[**HEDGES**](https://www.pnas.org/content/117/31/18489.full), 
[**DNA Fountain**](https://www.science.org/doi/abs/10.1126/science.aaj2038) and 
[**Yin-Yang Code**](https://www.biorxiv.org/content/10.1101/829721v3), 
the rough (not academic) conclusions are as follows:

<p align="center">
<img src="./experiments/results/figures/[6-0] result overviews.png" title="result overview" width="100%"/>
</p>


## Availability of process data

Please do not hesitate to contact zhanghaoling[at]genomics.cn.

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