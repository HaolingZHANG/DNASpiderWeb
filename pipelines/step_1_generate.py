__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import numpy
from dsw.generator import find_vertices
from dsw.generator import connect_default_graph, connect_fixed_graph, connect_variable_graph
from dsw.biofilter import LocalBioFilter
from pipelines import draw_graph_information, draw_fixed_graph_evolution


# Examples of restriction enzymes.
# Richard J. Roberts (1983) Nucleic Acids Research
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


# Similar structure of nucleotide affects Nanopore Sequencing's recognition of it.
# Ian M. Derrington et al. (2010) Proceedings of the National Academy of Sciences
# Ryan R. Wick et al. (2019) Genome biology
nanopore_segments = [
    "AGA", "GAG", "CTC", "TCT"
]


# When the stride exceeds 5, it takes at least 20 minutes to complete a calculation iteration.
# However, it usually takes 10 iteration to calculate whether the graph with specific out-degrees is feasible.
maximum_stride = 5


if __name__ == "__main__":
    bio_filters = {
        1: LocalBioFilter(max_homopolymer_runs=2, max_gc_bias=0.0, ignore_motifs=cut_segments),
        2: LocalBioFilter(max_homopolymer_runs=1),  # Goldman method
        3: LocalBioFilter(max_homopolymer_runs=2, max_gc_bias=0.1, ignore_motifs=nanopore_segments),
        4: LocalBioFilter(max_homopolymer_runs=2, max_gc_bias=0.1),
        5: LocalBioFilter(max_homopolymer_runs=4, max_gc_bias=0.1),
        6: LocalBioFilter(max_homopolymer_runs=3)  # Church method
    }

    path_groups = {}
    for index, bio_filter in bio_filters.items():
        print("Calculate graph " + str(index) + ".")
        find_vertices(length=10, bio_filter=bio_filter, save_path="../entities/BM" + str(index))
        vertices = numpy.load(file="../entities/BM" + str(index) + "[v].npy")
        connect_default_graph(length=10, vertices=vertices, save_path="../entities/BM" + str(index))
        connect_fixed_graph(length=10, vertices=vertices, maximum_stride=5, save_path="../entities/FLC" + str(index))
        connect_variable_graph(length=10, vertices=vertices, threshold=1, save_path="../entities/VLC" + str(index))
        print()

    draw_graph_information()
    draw_fixed_graph_evolution(maximum_stride)
