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
|cite| Zhi Ping et al. (2021) Nature Computational Science
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


local_bio_filters = {
    "01": LocalBioFilter(observed_length=10, max_homopolymer_runs=2, gc_range=balanced_gc_range,
                         undesired_motifs=cut_segments),
    "02": LocalBioFilter(observed_length=10, max_homopolymer_runs=1),
    "03": LocalBioFilter(observed_length=10, gc_range=low_gc_bias),
    "04": LocalBioFilter(observed_length=10, max_homopolymer_runs=2, gc_range=accepted_gc_bias,
                         undesired_motifs=nanopore_segments),
    "05": LocalBioFilter(observed_length=10, max_homopolymer_runs=2, gc_range=accepted_gc_bias),
    "06": LocalBioFilter(observed_length=10, gc_range=high_gc_bias),
    "07": LocalBioFilter(observed_length=10, max_homopolymer_runs=3, gc_range=accepted_gc_bias),
    "08": LocalBioFilter(observed_length=10, max_homopolymer_runs=4, gc_range=accepted_gc_bias),
    "09": LocalBioFilter(observed_length=10, max_homopolymer_runs=3),
    "10": LocalBioFilter(observed_length=10, max_homopolymer_runs=4),
    "11": LocalBioFilter(observed_length=10, max_homopolymer_runs=5),
    "12": LocalBioFilter(observed_length=10, max_homopolymer_runs=6)
}
