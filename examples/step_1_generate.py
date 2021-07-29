__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from dsw.generate import find_vertices, connect_graph
from dsw.screen import LocalScreen


# Examples of restriction enzymes.
# Roberts, R. J. (1983). Nucleic acids research.
cut_segments = [
    "GACGC",  # HgaI: Haemophilus gallinarum
    "AGCT",  # AluI*： Arthrobacter luteus
    "GATATC",  # EcoRV*: Escherichia coli
    "CAGCAG",  # EcoP15I: Escherichia coli
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

# similar
# Derrington, I. M., et al. (2010). Proceedings of the National Academy of Sciences.
# Wick, R. R., et al. (2019). Genome biology.
nanopore_segments = [
    "AGA", "GAG", "CTC", "TCT"
]


if __name__ == "__main__":
    screens = {"1": LocalScreen(max_homopolymer_runs=2, max_gc_bias=0.0, ignore_motifs=cut_segments),
               "2": LocalScreen(max_homopolymer_runs=1),
               "3": LocalScreen(max_homopolymer_runs=2, max_gc_bias=0.1, ignore_motifs=nanopore_segments),
               "4": LocalScreen(max_homopolymer_runs=2, max_gc_bias=0.1),
               "5": LocalScreen(max_homopolymer_runs=4, max_gc_bias=0.1),
               "6": LocalScreen(max_homopolymer_runs=3)}

    for screen_index, screen in screens.items():
        print("Calculate graph " + str(screen_index) + ".")
        vertices = find_vertices(length=10, screen=screen, save_path="../outputs/screen" + screen_index)
        connect_graph(length=10, vertices=vertices, threshold=1, save_path="../outputs/screen" + screen_index)
        print()
