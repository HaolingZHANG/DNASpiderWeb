__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import os
import random
import numpy
from collections import Counter
from matplotlib import pyplot
from dsw import b_system
from dsw.coder import encode
from dsw.repairer import repair_oligo
from dsw.monitor import Monitor
from pipelines import colors


def examples():
    print("Case estimate the fixed-length graph coder created by filter 1.")
    vertices = numpy.load(file="../entities/FLC1[vertices].npy")
    graph = numpy.load(file="../entities/FLC1[graph].npy")
    start_index = numpy.where(vertices == 1)[0][0]
    bits = b_system * 20
    wrong_location = 12
    right_oligo = encode(bits=bits, graph=graph, start_index=start_index)
    print("Right: " + right_oligo)
    print()
    for nucleotide in list(filter(lambda n: n != right_oligo[wrong_location], ["A", "C", "G", "T"])):
        wrong_oligo = list(right_oligo)
        wrong_oligo[wrong_location] = nucleotide
        wrong_oligo = "".join(wrong_oligo)
        print("Wrong: " + wrong_oligo)
        repaired_oligos = repair_oligo(graph=graph, oligo=wrong_oligo, start_index=start_index, need_logs=True)
        print("Repaired:")
        for repaired_oligo in repaired_oligos:
            print("       " + repaired_oligo + " | repair successfully? " + str(repaired_oligo == right_oligo))
        print()

    print("Case estimate the variable-length graph coder created by filter 1.")
    vertices = numpy.load(file="../entities/VLC1[vertices].npy")
    graph = numpy.load(file="../entities/VLC1[graph].npy")
    start_index = numpy.where(vertices == 1)[0][0]
    bits = b_system * 20
    wrong_location = 12
    right_oligo = encode(bits=bits, graph=graph, start_index=start_index)
    print("Right: " + right_oligo)
    print()
    for nucleotide in list(filter(lambda n: n != right_oligo[wrong_location], ["A", "C", "G", "T"])):
        wrong_oligo = list(right_oligo)
        wrong_oligo[wrong_location] = nucleotide
        wrong_oligo = "".join(wrong_oligo)
        print("Wrong: " + wrong_oligo)
        repaired_oligos = repair_oligo(graph=graph, oligo=wrong_oligo, start_index=start_index, need_logs=True)
        print("Repaired:")
        for repaired_oligo in repaired_oligos:
            print("       " + repaired_oligo + " | repair successfully? " + str(repaired_oligo == right_oligo))
        print()


def calculate_totally(task_seed, graph_type, filter_index, wrong_location, bit_length, times):
    """
    Calculate self-repair ability of the specific graph coder.

    :param task_seed: seed of the task (for generating tested bit matrix).
    :type task_seed: int

    :param graph_type: type of the graph (fixed-length code or variable-length code).
    :type graph_type: str

    :param filter_index: index of the filter.
    :type filter_index: int

    :param wrong_location: location in substitution error in the oligo.
    :type wrong_location: int

    :param bit_length: length of bit array.
    :type bit_length: int

    :param times: generating times of bit array.
    :type times: int
    """
    random.seed(task_seed)
    matrix = [[random.choice(b_system) for _ in range(bit_length)] for _ in range(times)]
    if not os.path.exists("../outputs/" + graph_type + "-" + str(filter_index) + " repair.npy"):
        monitor = Monitor()
        print("Totally calculate " + graph_type + " created by filter " + str(filter_index) + ".")
        vertices = numpy.load(file="../entities/" + graph_type + str(filter_index) + "[vertices].npy")
        graph = numpy.load(file="../entities/" + graph_type + str(filter_index) + "[graph].npy")
        start_indices = numpy.where(vertices == 1)[0]
        results = numpy.zeros(shape=(len(start_indices), times * 3))
        for current, start_index in enumerate(start_indices):
            position = 0
            for array_index, bits in enumerate(matrix):
                right_oligo = encode(bits=bits, graph=graph, start_index=start_index)
                for nucleotide in list(filter(lambda n: n != right_oligo[wrong_location], ["A", "C", "G", "T"])):
                    wrong_oligo = list(right_oligo)
                    wrong_oligo[wrong_location] = nucleotide
                    wrong_oligo = "".join(wrong_oligo)
                    recover_oligos = repair_oligo(graph=graph, oligo=wrong_oligo, start_index=start_index,
                                                  need_logs=False)
                    found = right_oligo in recover_oligos
                    if found:
                        results[current][position] = len(recover_oligos)
                    else:
                        results[current][position] = 0
                    position += 1
            monitor.output(current + 1, len(start_indices))
        numpy.save(file="../outputs/" + graph_type + "-" + str(filter_index) + " repair.npy", arr=results)


def analyze():
    results = []
    for screen_index in [1, 2, 3, 4, 5, 6]:
        repaired_results = numpy.load(file="../outputs/FLC-" + str(screen_index) + " repair.npy")
        repaired_results = repaired_results.reshape(-1)
        counter = Counter(repaired_results.tolist())
        result = [0 for _ in range(21)]
        for repaired_oligo_number in range(21):
            if repaired_oligo_number in counter:
                result[repaired_oligo_number] = counter[repaired_oligo_number]
        results.append(numpy.array(result))

    pyplot.figure(figsize=(10, 4.8))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.25, hspace=0.5)
    for index, result in enumerate(results):
        pyplot.subplot(2, 3, index + 1)
        pyplot.title("VLC-" + str(index + 1), fontsize=12)
        pyplot.bar(range(21), result / numpy.sum(result), width=1, color=colors["fixed"], edgecolor="black",
                   linewidth=0.5, zorder=1)
        for value in [0.2, 0.4, 0.6, 0.8]:
            pyplot.hlines(value, -0.5, 20.5, color="silver", linewidth=0.5, linestyle="--", zorder=0)
        pyplot.xlabel("repaired oligo number", fontsize=10)
        pyplot.ylabel("percentage (%)", fontsize=10)
        pyplot.xticks(range(21), range(21), fontsize=8)
        pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0", "20", "40", "60", "80", "100"], fontsize=8)
        pyplot.xlim(-0.5, 20.5)
        pyplot.ylim(0, 1)
    pyplot.savefig("../outputs/self-repair FLC statistics.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()

    results = []
    for screen_index in [1, 2, 3, 4, 5, 6]:
        repaired_results = numpy.load(file="../outputs/VLC-" + str(screen_index) + " repair.npy")
        repaired_results = repaired_results.reshape(-1)
        counter = Counter(repaired_results.tolist())
        result = [0 for _ in range(21)]
        for repaired_oligo_number in range(21):
            if repaired_oligo_number in counter:
                result[repaired_oligo_number] = counter[repaired_oligo_number]
        results.append(numpy.array(result))

    pyplot.figure(figsize=(10, 4.8))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.25, hspace=0.5)
    for index, result in enumerate(results):
        pyplot.subplot(2, 3, index + 1)
        pyplot.title("VLC-" + str(index + 1), fontsize=12)
        pyplot.bar(range(21), result / numpy.sum(result), width=1, color=colors["variable"], edgecolor="black",
                   linewidth=0.5, zorder=1)
        for value in [0.2, 0.4, 0.6, 0.8]:
            pyplot.hlines(value, -0.5, 20.5, color="silver", linewidth=0.5, linestyle="--", zorder=0)
        pyplot.xlabel("repaired oligo number", fontsize=10)
        pyplot.ylabel("percentage (%)", fontsize=10)
        pyplot.xticks(range(21), range(21), fontsize=8)
        pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0", "20", "40", "60", "80", "100"], fontsize=8)
        pyplot.xlim(-0.5, 20.5)
        pyplot.ylim(0, 1)
    pyplot.savefig("../outputs/self-repair VLC statistics.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


if __name__ == "__main__":
    examples()
    for t in ["FLC", "VLC"]:
        for fi in [1, 2, 3, 4, 5, 6]:
            calculate_totally(task_seed=2021, graph_type=t, filter_index=fi,
                              wrong_location=12, bit_length=50, times=100)
    analyze()
