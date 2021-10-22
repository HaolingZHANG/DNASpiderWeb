__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from numpy import load
from dsw import find_vertices
from dsw import connect_default_graph, connect_fixed_graph, connect_variable_graph
from dsw import LocalBioFilter


"""
Examples of restriction enzyme sites.
|cite| Richard J. Roberts (1983) Nucleic Acids Research
"""
cut_segments = [
    "AGCT",    # AluI*:    Arthrobacter luteus
    "GACGC",   # HgaI:     Haemophilus gallinarum
    "CAGCAG",  # EcoP15I:  Escherichia coli
    "GATATC",  # EcoRV*:   Escherichia coli
    "GGTACC",  # KpnI:     Klebsiella pneumoniae
    "CTGCAG",  # PstI:     Providencia stuartii
    "GAGCTC",  # SacI:     Streptomyces achromogenes
    "GTCGAC",  # SalI:     Streptomyces albus
    "AGTACT",  # ScaI*:    Streptomyces caespitosus
    "ACTAGT",  # SpeI:     Sphaerotilus natans
    "GCATGC",  # SphI:     Streptomyces phaeochromogenes
    "AGGCCT",  # StuI*:    Streptomyces tubercidicus
    "TCTAGA"   # XbaI:     Xanthomonas badrii
]


"""
Similar structure of nucleotide affects Nanopore Sequencing's recognition of it.
|cite| Ian M. Derrington et al. (2010) Proceedings of the National Academy of Sciences
|cite| Ryan R. Wick et al. (2019) Genome biology
"""
nanopore_segments = [
    "AGA", "GAG", "CTC", "TCT"
]


"""
Local GC stability.
"""
balanced_gc_range = [0.5, 0.5]

"""
Well-accepted GC bias.
|cite| Yaniv Erlich and Dina Zielinski (2017) Science
|cite| Zhi Ping et al. (2021) bioRxiv
"""
accepted_gc_bias = [0.4, 0.6]

"""
High GC content range.
|cite| Jake L. Weissman et al. (2019) PLoS Genetics
"""
high_gc_bias = [0.5, 0.7]

"""
Low GC content range.
|cite| Anne-Kathrin Dietel et al. (2019) PLoS Genetics
"""
low_gc_bias = [0.1, 0.3]

"""
When the stride exceeds 5, it takes at least 20 minutes to complete a calculation iteration.
However, it usually takes 10 iteration to calculate whether the graph with specific out-degrees is feasible.
"""
maximum_stride = 5


def obtain_filters():
    return {
        1: LocalBioFilter(max_homopolymer_runs=2, gc_range=balanced_gc_range, undesired_motifs=cut_segments),
        2: LocalBioFilter(max_homopolymer_runs=1),  # Goldman method
        3: LocalBioFilter(gc_range=low_gc_bias),
        4: LocalBioFilter(max_homopolymer_runs=2, gc_range=accepted_gc_bias, undesired_motifs=nanopore_segments),
        5: LocalBioFilter(max_homopolymer_runs=2, gc_range=accepted_gc_bias),
        6: LocalBioFilter(gc_range=high_gc_bias),
        7: LocalBioFilter(max_homopolymer_runs=4, gc_range=accepted_gc_bias),  # DNA Fountain and Yin-Yang Code
        8: LocalBioFilter(max_homopolymer_runs=3)  # Church method
    }


if __name__ == "__main__":
    filters = obtain_filters()
    for index, bio_filter in filters.items():
        print("Calculate graph " + str(index) + ".")
        find_vertices(length=10, bio_filter=bio_filter, save_path="../entities/BM" + str(index))
        vertices = load(file="../entities/BM" + str(index) + "[v].npy")
        connect_default_graph(length=10, vertices=vertices, save_path="../entities/BM" + str(index))
        connect_fixed_graph(length=10, vertices=vertices, maximum_stride=5, save_path="../entities/FLC" + str(index))
        connect_variable_graph(length=10, vertices=vertices, threshold=1, save_path="../entities/VLC" + str(index))
        print()
