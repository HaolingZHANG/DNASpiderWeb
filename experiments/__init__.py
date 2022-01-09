from os import path, makedirs

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


def obtain_filters():
    return {
        "01": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=2,
                             gc_range=balanced_gc_range,
                             undesired_motifs=cut_segments),
        "02": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=1),  # Goldman method
        "03": LocalBioFilter(observed_length=10,
                             gc_range=low_gc_bias),
        "04": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=2,
                             gc_range=accepted_gc_bias,
                             undesired_motifs=nanopore_segments),
        "05": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=2,
                             gc_range=accepted_gc_bias),
        "06": LocalBioFilter(observed_length=10,
                             gc_range=high_gc_bias),
        "07": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=3,
                             gc_range=accepted_gc_bias),
        "08": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=4,
                             gc_range=accepted_gc_bias),  # HEDGES Code
        "09": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=3),  # Church method
        "10": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=4),
        "11": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=5),
        "12": LocalBioFilter(observed_length=10,
                             max_homopolymer_runs=6)
    }


def create_folders():
    if not path.exists("./results/"):
        makedirs("./results/")
    if not path.exists("./results/data/"):
        makedirs("./results/data/")
    if not path.exists("./results/figures/"):
        makedirs("./results/figures/")
    if not path.exists("./results/temp/"):
        makedirs("./results/temp/")


constraint_names = ["in vivo\n(site)", "goldman\norigin", "in vivo\n(10-30)",
                    "nanopore\nA/G and C/T", "fountain\n(-2)", "in vivo\n(50-70)",
                    "fountain\n(-1)", "fountain\norigin", "church\norigin",
                    "church\n(+1)", "church\n(+2)", "church\n(+3)"]


colors = {
    "trad1": "#81B8DF",
    "trad2": "#B1CCDF",
    "algo1": "#FE817D",
    "algo2": "#FCBBAE",
    "yyco1": "#00C382",
    "yyco2": "#ADEDD9",
    "foun1": "#AD94FF",
    "foun2": "#D3C2FE",
    "hedc1": "#FFC000",
    "hedc2": "#FFDB75",
    "diffs": "#F1F1F1",
    "goodG": "#BAE49A",
    "normG": "#FFDF5F",
    "hardG": "#FB875E"
}
