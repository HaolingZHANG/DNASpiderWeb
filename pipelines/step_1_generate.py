__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import math
import numpy
from matplotlib import pyplot
from dsw.generator import find_vertices
from dsw.generator import connect_default_graph, connect_fixed_graph, connect_variable_graph
from dsw.biofilter import LocalBioFilter
from pipelines import colors

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


def display():
    pyplot.figure(figsize=(10, 2.4))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.25)

    # calculate vertex number.
    default, variable, fixed = [], [], []
    for filter_index in ["1", "2", "3", "4", "5", "6"]:
        default_vertices = numpy.load(file="../entities/Default" + str(filter_index) + "[vertices].npy")
        default.append(math.log(numpy.sum(default_vertices), 4))
        variable_vertices = numpy.load(file="../entities/VLC" + str(filter_index) + "[vertices].npy")
        variable.append(math.log(numpy.sum(variable_vertices), 4))
        fixed_vertices = numpy.load(file="../entities/FLC" + str(filter_index) + "[vertices].npy")
        fixed.append(math.log(numpy.sum(fixed_vertices), 4))
    default = numpy.array(default)
    variable = numpy.array(variable)
    fixed = numpy.array(fixed)

    pyplot.subplot(1, 2, 1)
    pyplot.title("vertex number", fontsize=12)
    x = numpy.array(range(6)) + 1
    pyplot.bar(x - 0.2, default - 5, width=0.2, color=colors["default"], edgecolor="black", linewidth=0.5, zorder=1)
    pyplot.bar(x, variable - 5, width=0.2, color=colors["variable"], edgecolor="black", linewidth=0.5, zorder=1)
    pyplot.bar(x + 0.2, fixed - 5, width=0.2, color=colors["fixed"], edgecolor="black", linewidth=0.5, zorder=1)
    for value in [0, 1, 2, 3, 4, 5]:
        pyplot.hlines(value, 0.5, 6.5, color="silver", linewidth=0.5, linestyle="--", zorder=0)
    pyplot.xlabel("biochemical requirement index", fontsize=10)
    pyplot.ylabel("log4(vertex number)", fontsize=10)
    pyplot.xticks(x, x, fontsize=8)
    pyplot.yticks([0, 1, 2, 3, 4, 5], [5, 6, 7, 8, 9, 10], fontsize=8)
    pyplot.xlim(0.5, 6.5)
    pyplot.ylim(-0.3, 5.3)

    # calculate out-degree.
    default, variable, fixed = [], [], []
    for filter_index in ["1", "2", "3", "4", "5", "6"]:
        default_vertices = numpy.load(file="../entities/Default" + str(filter_index) + "[vertices].npy")
        default_transform = numpy.load(file="../entities/Default" + str(filter_index) + "[graph].npy")
        default_transform[default_transform > 0] = 0
        average = (4 ** 11 + numpy.sum(default_transform)) / numpy.sum(default_vertices)
        minimum, maximum = 4, 0
        for vertex in default_transform:
            count = len(numpy.where(vertex >= 0)[0])
            if count > 0:
                minimum = min(minimum, count)
                maximum = max(maximum, count)
        default.append([minimum, average, maximum])
        variable_vertices = numpy.load(file="../entities/VLC" + str(filter_index) + "[vertices].npy")
        variable_transform = numpy.load(file="../entities/VLC" + str(filter_index) + "[graph].npy")
        variable_transform[variable_transform > 0] = 0
        average = (4 ** 11 + numpy.sum(variable_transform)) / numpy.sum(variable_vertices)
        minimum, maximum = 4, 1
        for vertex in variable_transform:
            count = len(numpy.where(vertex >= 0)[0])
            if count > 0:
                minimum = min(minimum, count)
                maximum = max(maximum, count)
        variable.append([minimum, average, maximum])
        fixed_vertices = numpy.load(file="../entities/FLC" + str(filter_index) + "[vertices].npy")
        fixed_transform = numpy.load(file="../entities/FLC" + str(filter_index) + "[graph].npy")
        fixed_transform[fixed_transform > 0] = 0
        average = (4 ** 11 + numpy.sum(fixed_transform)) / numpy.sum(fixed_vertices)
        minimum, maximum = 4, 1
        for vertex in fixed_transform:
            count = len(numpy.where(vertex >= 0)[0])
            if count > 0:
                minimum = min(minimum, count)
                maximum = max(maximum, count)
        fixed.append([minimum, average, maximum])
    default = numpy.array(default).T
    variable = numpy.array(variable).T
    fixed = numpy.array(fixed).T

    pyplot.subplot(1, 2, 2)
    pyplot.title("composition of out-degree", fontsize=12)
    x = numpy.array(range(6)) + 1
    pyplot.bar(x - 0.2, default[1], width=0.2, color=colors["default"], edgecolor="black", linewidth=0.5, zorder=1)
    for position, minimum, maximum in zip(x, default[0], default[2]):
        pyplot.vlines(position - 0.2, minimum, maximum, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(minimum, position - 0.25, position - 0.15, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(maximum, position - 0.25, position - 0.15, color="black", linewidth=1.5, zorder=2)
    pyplot.bar(x, variable[1], width=0.2, color=colors["variable"], edgecolor="black", linewidth=0.5, zorder=1)
    for position, minimum, maximum in zip(x, variable[0], variable[2]):
        pyplot.vlines(position, minimum, maximum, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(minimum, position - 0.05, position + 0.05, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(maximum, position - 0.05, position + 0.05, color="black", linewidth=1.5, zorder=2)
    pyplot.bar(x + 0.2, fixed[1], width=0.2, color=colors["fixed"], edgecolor="black", linewidth=0.5, zorder=1)
    for position, minimum, maximum in zip(x, fixed[0], fixed[2]):
        pyplot.vlines(position + 0.2, minimum, maximum, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(minimum, position + 0.15, position + 0.25, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(maximum, position + 0.15, position + 0.25, color="black", linewidth=1.5, zorder=2)
    for value in [0, 1, 2, 3, 4]:
        pyplot.hlines(value, 0.5, 6.5, color="silver", linewidth=0.5, linestyle="--", zorder=0)
    pyplot.xlabel("biochemical requirement index", fontsize=10)
    pyplot.ylabel("out-degree", fontsize=10)
    pyplot.xticks(x, x, fontsize=8)
    pyplot.yticks([0, 1, 2, 3, 4], [0, 1, 2, 3, 4], fontsize=8)
    pyplot.xlim(0.5, 6.5)
    pyplot.ylim(-0.25, 4.25)
    pyplot.savefig("../outputs/generation statistics.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


if __name__ == "__main__":
    bio_filters = {"1": LocalBioFilter(max_homopolymer_runs=2, max_gc_bias=0.0, ignore_motifs=cut_segments),
                   "2": LocalBioFilter(max_homopolymer_runs=1),  # Goldman method
                   "3": LocalBioFilter(max_homopolymer_runs=2, max_gc_bias=0.1, ignore_motifs=nanopore_segments),
                   "4": LocalBioFilter(max_homopolymer_runs=2, max_gc_bias=0.1),
                   "5": LocalBioFilter(max_homopolymer_runs=4, max_gc_bias=0.1),
                   "6": LocalBioFilter(max_homopolymer_runs=3)}  # Church method

    for index, bio_filter in bio_filters.items():
        print("Calculate graph " + str(index) + ".")
        vertices = find_vertices(length=10, bio_filter=bio_filter, save_path="../entities/Default" + index)
        connect_default_graph(length=10, vertices=vertices, save_path="../entities/Default" + index)
        connect_fixed_graph(length=10, vertices=vertices, save_path="../entities/FLC" + index)
        connect_variable_graph(length=10, vertices=vertices, threshold=1, save_path="../entities/VLC" + index)
        print()

    display()
