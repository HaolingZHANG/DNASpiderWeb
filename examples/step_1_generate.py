__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from dsw.generator import find_vertices, connect_graph
from dsw.biofilter import LocalBioFilter


# Examples of restriction enzymes.
# Roberts, R. J. (1983). Nucleic acids research.
cut_segments = [
    "AGCT",  # AluI*： Arthrobacter luteus
    "GACGC",  # HgaI: Haemophilus gallinarum
    "CAGCAG",  # EcoP15I: Escherichia coli
    "GATATC",  # EcoRV*: Escherichia coli
    "GGTACC",  # KpnI: Klebsiella pneumoniae
    "CTGCAG",  # PstI: Providencia stuartii
    "GAGCTC",  # SacI: Streptomyces achromogenes
    "GTCGAC",  # SalI: Streptomyces albus
    "AGTACT",  # ScaI*: Streptomyces caespitosus
    "ACTAGT",  # SpeI: Sphaerotilus natans
    "GCATGC",  # SphI: Streptomyces phaeochromogenes
    "AGGCCT",  # StuI*: Streptomyces tubercidicus
    "TCTAGA"  # XbaI: Xanthomonas badrii
]


# Similar structure of nucleotide affects Nanopore Sequencing's recognition of it.
# Derrington, I. M., et al. (2010). Proceedings of the National Academy of Sciences.
# Wick, R. R., et al. (2019). Genome biology.
nanopore_segments = [
    "AGA", "GAG", "CTC", "TCT"
]


if __name__ == "__main__":
    bio_filters = {"1": LocalBioFilter(max_homopolymer_runs=2, max_gc_bias=0.0, ignore_motifs=cut_segments),
                   "2": LocalBioFilter(max_homopolymer_runs=1),
                   "3": LocalBioFilter(max_homopolymer_runs=2, max_gc_bias=0.1, ignore_motifs=nanopore_segments),
                   "4": LocalBioFilter(max_homopolymer_runs=2, max_gc_bias=0.1),
                   "5": LocalBioFilter(max_homopolymer_runs=4, max_gc_bias=0.1),
                   "6": LocalBioFilter(max_homopolymer_runs=3)}

    for screen_index, bio_filter in bio_filters.items():
        print("Calculate graph " + str(screen_index) + ".")
        vertices = find_vertices(length=10, bio_filter=bio_filter, save_path="../outputs/constraint" + screen_index)
        connect_graph(length=10, vertices=vertices, threshold=1, save_path="../outputs/constraint" + screen_index)
        print()
