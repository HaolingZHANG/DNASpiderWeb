__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import os
import random
import numpy
from matplotlib import pyplot
from dsw.coder import encode, decode, obtain_fixed_length
from dsw.evaluator import approximate_upper_bound
from dsw.monitor import Monitor
from pipelines import colors


def transcode(name, graph, start_indices, task_seed, bit_length, test_time, fixed_length=None):
    """
    Run the encoding and decoding pipeline.

    :param name: name of this graph coder.
    :type name: str

    :param graph: graph of DNA Spider-Web.
    :type graph: numpy.ndarray

    :param start_indices: virtual vertices to start decoding (consistent with the encoding process).
    :type start_indices: list

    :param task_seed: task seed for random generating the bit matrix.
    :type task_seed: int

    :param bit_length: length of bit array.
    :type bit_length: int

    :param test_time: transcoding test iteration.
    :type test_time: int

    :param fixed_length: length of oligo (it is required for fixed-length coder).
    :type fixed_length: int

    :raise: ValueError, if encode bit matrix is not equal to decoded bit matrix.

    :returns: decoded bit matrix and nucleotide number during encoding process.
    :rtype: tuple
    """
    random.seed(task_seed)
    encoded_matrix = [[random.randint(0, 1) for _ in range(bit_length)] for _ in range(test_time)]

    print("Evaluate " + name + ".")
    monitor, nucleotide_numbers = Monitor(), []
    for current, start_index in enumerate(start_indices):
        oligos, nucleotide_number = [], 0
        for bits in encoded_matrix:
            oligo = encode(bits=bits, graph=graph, start_index=start_index, length=fixed_length)
            nucleotide_number += len(oligo)
            oligos.append(oligo)

        nucleotide_numbers.append(nucleotide_number)

        decoded_matrix = []
        for oligo in oligos:
            decoded_matrix.append(decode(oligo=oligo, bit_length=bit_length, graph=graph, start_index=start_index))

        if encoded_matrix != decoded_matrix:
            raise ValueError("Encoded matrix is not equal to decoded matrix!")

        monitor.output(current + 1, len(start_indices))

    return nucleotide_numbers


def display():
    pyplot.figure(figsize=(10, 4.8))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.25, hspace=0.5)
    for filter_index in [1, 2, 3, 4, 5, 6]:
        pyplot.subplot(2, 3, filter_index)
        pyplot.title("biochemical filter " + str(filter_index), fontsize=12)
        d_bound = approximate_upper_bound(numpy.load(file="../entities/Default" + str(filter_index) + "[graph].npy"))
        v_bound = approximate_upper_bound(numpy.load(file="../entities/VLC" + str(filter_index) + "[graph].npy"))
        f_bound = approximate_upper_bound(numpy.load(file="../entities/FLC" + str(filter_index) + "[graph].npy"))

        v_trans = numpy.load("../outputs/VLC-" + str(filter_index) + " transcode.npy")
        f_trans = numpy.load("../outputs/FLC-" + str(filter_index) + " transcode.npy")

        pyplot.text(x=0, y=d_bound + 0.01, s=str(round(d_bound, 3)), fontsize=8,
                    horizontalalignment="center", verticalalignment="bottom")
        pyplot.hlines(d_bound, -0.25, 0.25, color="black", linewidth=1.5)
        pyplot.text(x=1, y=v_bound + 0.01, s=str(round(v_bound, 3)), fontsize=8,
                    horizontalalignment="center", verticalalignment="bottom")
        pyplot.hlines(v_bound, 0.75, 1.25, color=colors["variable"], linewidth=1.5)
        pyplot.text(x=2, y=f_bound + 0.01, s=str(round(f_bound, 3)), fontsize=8,
                    horizontalalignment="center", verticalalignment="bottom")
        pyplot.hlines(f_bound, 1.75, 2.25, color=colors["fixed"], linewidth=1.5)

        violin = pyplot.violinplot(dataset=v_trans, positions=[1], showextrema=False)
        for patch in violin["bodies"]:
            patch.set_facecolor(colors["variable"])
            patch.set_edgecolor("black")
            patch.set_linewidth(0.5)
            patch.set_alpha(1)
        min_value, median, max_value = numpy.percentile(v_trans, [0, 50, 100])
        pyplot.scatter([1], median, color=colors["variable"], s=5, linewidths=1, edgecolors="black", zorder=4)

        violin = pyplot.violinplot(dataset=f_trans, positions=[2], showextrema=False)
        for patch in violin["bodies"]:
            patch.set_facecolor(colors["fixed"])
            patch.set_edgecolor("black")
            patch.set_linewidth(0.5)
            patch.set_alpha(1)
        min_value, median, max_value = numpy.percentile(f_trans, [0, 50, 100])
        pyplot.scatter([2], median, color=colors["fixed"], s=5, linewidths=1, edgecolors="black", zorder=4)

        for value in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            pyplot.hlines(value, -0.5, 2.5, color="silver", linewidth=0.5, linestyle="--", zorder=0)
        pyplot.xlabel("upper bound & actual performances", fontsize=10)
        pyplot.ylabel("capacity", fontsize=10)
        pyplot.xticks(range(3), ["benchmark", "VLC", "FLC"], fontsize=8)
        pyplot.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0], [0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize=8)
        pyplot.xlim(-0.5, 2.5)
        pyplot.ylim(0.425, 1.075)
    pyplot.savefig("../outputs/transcode statistics.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


if __name__ == "__main__":
    l, t = 100, 100
    for index in [1, 2, 3, 4, 5, 6]:
        if not os.path.exists("../outputs/FLC-" + str(index) + " transcode.npy"):
            g = numpy.load(file="../entities/FLC" + str(index) + "[graph].npy")
            v = numpy.where(numpy.load(file="../entities/FLC" + str(index) + "[vertices].npy") == 1)[0]
            numbers = transcode(name="FLC-" + str(index), graph=g, start_indices=v,
                                task_seed=2021, bit_length=l, test_time=t,
                                fixed_length=obtain_fixed_length(graph=g, name="FLC-" + str(index), bit_length=l))
            numbers = ((l * t) / numpy.array(numbers)) / 2
            numpy.save(file="../outputs/FLC-" + str(index) + " transcode.npy", arr=numbers)

        if not os.path.exists("../outputs/VLC-" + str(index) + " transcode.npy"):
            g = numpy.load(file="../entities/VLC" + str(index) + "[graph].npy")
            v = numpy.where(numpy.load(file="../entities/VLC" + str(index) + "[vertices].npy") == 1)[0]
            numbers = transcode(name="VLC-" + str(index), graph=g, start_indices=v,
                                task_seed=2021, bit_length=l, test_time=t)
            numbers = ((l * t) / numpy.array(numbers)) / 2
            numpy.save(file="../outputs/VLC-" + str(index) + " transcode.npy", arr=numbers)
    display()
