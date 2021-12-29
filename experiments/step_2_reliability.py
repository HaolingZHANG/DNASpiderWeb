__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from itertools import product
from math import log
# noinspection PyPackageRequirements
from matplotlib import pyplot
from numpy import array, zeros, zeros_like, random, linspace, where
from numpy import sum, min, max, log10, linalg, percentile, real, polyfit, std, sqrt, mean, abs
from os import path
from pickle import load, dump

from dsw import Monitor, approximate_capacity, dna_to_number
from dsw import accessor_to_adjacency_matrix, adjacency_matrix_to_accessor, get_complete_accessor

from experiments import colors, create_folders


def display_cases():

    def create_matrix(dna_strings):
        matrix = zeros(shape=(16, 16), dtype=int)
        vertices = [dna_to_number(dna_string=vertex, is_string=False) for vertex in dna_strings]
        for former, latter in product(range(len(dna_strings)), repeat=2):
            if dna_strings[former][1] == dna_strings[latter][0]:
                matrix[vertices[former], vertices[latter]] = 1
        return matrix

    cases = [create_matrix(dna_strings=["AC",
                                        "CG",
                                        "GT",
                                        "TA"]),  # a cycle.
             create_matrix(dna_strings=["AC", "AG",
                                        "CA", "CG",
                                        "GA", "GT",
                                        "TC", "TG"]),  # GC content is 50%.
             create_matrix(dna_strings=["AC", "AG", "AT",
                                        "CA", "CG", "CT",
                                        "GA", "GC", "GT",
                                        "TA", "TC", "TG"]),  # no homopolymer.
             create_matrix(dna_strings=["AA", "AC", "AG", "AT",
                                        "CA", "CC", "CG", "CT",
                                        "GA", "GC", "GG", "GT",
                                        "TA", "TC", "TG", "TT"])]  # complete graph.

    figure = pyplot.figure(figsize=(10, 5), tight_layout=True)
    pyplot.subplots_adjust(wspace=0.15, hspace=0.3)
    for index, case in enumerate(cases):
        pyplot.subplot(2, len(cases), index + 1)
        for x in range(16):
            for y in range(16):
                if case[x, y] == 1:
                    pyplot.fill_between(x=[x, x + 1], y1=[y, y], y2=[y + 1, y + 1], color=colors["trad1"])
        pyplot.xticks(array(list(range(16))) + 0.5,
                      ["0", "", "", "3", "", "", "6", "", "", "9", "", "", "12", "", "", "15"], fontsize=8)
        if index > 0:
            pyplot.yticks(array(list(range(16))) + 0.5,
                          ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""])
        else:
            pyplot.yticks(array(list(range(16))) + 0.5,
                          ["0", "", "", "3", "", "", "6", "", "", "9", "", "", "12", "", "", "15"], fontsize=8)
            pyplot.ylabel("index of latter DNA string", fontsize=8)
        for value in range(16):
            pyplot.vlines(value, 0, 16, color="black", linewidth=0.5)
            pyplot.hlines(value, 0, 16, color="black", linewidth=0.5)
        pyplot.xlabel("index of former DNA string", fontsize=8)
        pyplot.xlim(0, 16)
        pyplot.ylim(0, 16)
        axes = pyplot.subplot(2, len(cases), index + len(cases) + 1)
        _, process = approximate_capacity(accessor=adjacency_matrix_to_accessor(case), repeats=1, need_process=True)
        if len(process) == 1:
            pyplot.scatter([0], process, color=colors["trad1"])
        pyplot.plot(range(len(process)), process, color=colors["trad1"], marker="o", linewidth=2)
        pyplot.text(len(process) - 0.8, process[-1], "%.3f" % process[-1], va="center", ha="left")
        pyplot.xlim(-0.1, 2.1)
        pyplot.ylim(-0.1, 2.1)
        if index > 0:
            pyplot.yticks([0, 0.5, 1, 1.5, 2], ["", "", "", "", ""], fontsize=8)
        else:
            pyplot.ylabel("capacity", fontsize=8)
        pyplot.xticks([0, 1, 2], [1, 2, 3], fontsize=8)
        pyplot.xlabel("iteration", fontsize=8)
        # noinspection PyUnresolvedReferences
        axes.spines["right"].set_visible(False)
        # noinspection PyUnresolvedReferences
        axes.spines["top"].set_visible(False)

    figure.align_labels()

    pyplot.savefig("./results/figures/[2-1] reliability regular.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()


def evaluate_length_2(test_times):
    if not path.exists("./results/data/step_2_reliability_detail.pkl"):
        random.seed(2021)
        feed_matrix = accessor_to_adjacency_matrix(get_complete_accessor(observed_length=2))
        data, errors, monitor = [], [], Monitor()
        for time in range(test_times):
            matrix = feed_matrix.copy()
            for position in array(list(where(feed_matrix == 1)[:2])).T:
                if random.random(1)[0] <= 0.5:
                    matrix[position[0], position[1]] = 0
            monitor.output(current_state=time + 1, total_state=test_times)
            calculated_eigenvalue = real(max(linalg.eig(matrix)[0]))
            calculated_value = log(calculated_eigenvalue, 2) if calculated_eigenvalue > 0.0 else 0.0
            approximate_value = approximate_capacity(accessor=adjacency_matrix_to_accessor(matrix), repeats=10)
            data.append(approximate_value)
            errors.append(abs(approximate_value - calculated_value))
        with open("./results/data/step_2_reliability_detail.pkl", "wb") as file:
            dump((data, errors), file)
    else:
        with open("./results/data/step_2_reliability_detail.pkl", "rb") as file:
            data, errors = load(file)

    errors = log10(array(errors))

    violin = pyplot.violinplot(dataset=errors, positions=[0], showextrema=False, vert=False, widths=1)
    paint_data = violin["bodies"][0].get_paths()[0].vertices
    paint_data = paint_data[len(paint_data) // 2 + 1: -1]
    pyplot.close()

    x, y = paint_data[:, 0], paint_data[:, 1]
    y = y / sum(y) * 100.0

    [v_25, v_50, v_75] = percentile(a=errors, q=[25, 50, 75])
    lower, upper = v_25 - 1.5 * (v_75 - v_25), 1.5 * (v_75 - v_25) + v_75

    figure = pyplot.figure(figsize=(10, 5), tight_layout=True)
    pyplot.subplots_adjust(wspace=0.3, hspace=0.3)
    pyplot.subplot(2, 1, 1)

    pyplot.plot(x, y, color=colors["trad1"], linewidth=1.5, zorder=4)
    pyplot.fill_between(x, y, zeros_like(y), color=colors["trad2"], alpha=0.75, zorder=3)

    pyplot.vlines(x=v_50, ymin=0, ymax=6.5, color="black", linewidth=1, zorder=2)
    pyplot.vlines(x=v_50, ymin=6.9, ymax=8.5, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=8.5, xmin=v_50 - 0.3, xmax=v_50 + 0.3, color="black", linewidth=1, zorder=2)
    pyplot.text(x=v_50, y=8.6, s="median", va="bottom", ha="center", fontsize=8, zorder=2)

    pyplot.vlines(x=v_25, ymin=0, ymax=7, color="black", linewidth=1, zorder=2)
    pyplot.vlines(x=v_75, ymin=0, ymax=7, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=6.7, xmin=v_25, xmax=v_75, color="black", linewidth=1, zorder=2)
    pyplot.text(x=v_75 + 0.1, y=6.7, s="interquartile range", va="center", ha="left", fontsize=8, zorder=2)
    pyplot.vlines(x=lower, ymin=0, ymax=4, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=3.7, xmin=lower, xmax=min(errors) + 2, color="black", linewidth=1, zorder=2)
    pyplot.text(x=min(errors) + 1.9, y=3.7, s="outlier range", va="center", ha="right", fontsize=8, zorder=2)
    pyplot.vlines(x=upper, ymin=0, ymax=4, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=3.7, xmin=upper, xmax=max(errors) - 2, color="black", linewidth=1, zorder=2)
    pyplot.text(x=max(errors) - 1.9, y=3.7, s="outlier range", va="center", ha="left", fontsize=8, zorder=2)

    pyplot.xlim(-14.1, -3.9)
    pyplot.ylim(-0.5, 10.5)
    pyplot.ylabel("frequency", fontsize=8)
    pyplot.xticks([-14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4],
                  ["", "", "", "", "", "", "", "", "", "", ""], fontsize=8)
    pyplot.yticks([0, 2, 4, 6, 8, 10], ["0%", "2%", "4%", "6%", "8%", "10%"], fontsize=8)

    pyplot.subplot(2, 1, 2)
    pyplot.scatter(errors, data, color=colors["trad1"], edgecolor="black", zorder=3)

    temp_errors, temp_data = [], []
    for value_1, value_2 in sorted(zip(errors, data), key=lambda v: v[0]):
        temp_errors.append(value_1)
        temp_data.append(value_2)
    errors, data = array(temp_errors), array(temp_data)

    a, b = polyfit(errors, data, deg=1)
    estimated = a * errors + b
    bounds = abs(std(errors) * sqrt(1.0 / len(errors)
                                    + (errors - mean(errors)) ** 2 / sum((errors - mean(errors)) ** 2)))
    pyplot.fill_between(errors, estimated - bounds, estimated + bounds, color=colors["trad2"], alpha=0.75, zorder=1)
    x = linspace(min(errors), max(errors), 100)
    y = a * x + b
    pyplot.plot(x, y, color="black", linestyle="--", linewidth=1, zorder=2)

    pyplot.xlim(-14.1, -3.9)
    pyplot.ylim(-0.1, 2.1)
    pyplot.xticks([-14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4],
                  ["$10^{-14}$", "$10^{-13}$", "$10^{-12}$", "$10^{-11}$", "$10^{-10}$",
                   "$10^{-09}$", "$10^{-08}$", "$10^{-07}$", "$10^{-06}$", "$10^{-05}$", "$10^{-04}$"],
                  fontsize=8)
    pyplot.yticks([0, 0.5, 1.0, 1.5, 2.0], ["0.0", "0.5", "1.0", "1.5", "2.0"], fontsize=8)
    pyplot.xlabel("relative error with NumPy \"linalg.eig\" function", fontsize=8)
    pyplot.ylabel("capacity using SPIDER-WEB", fontsize=8)

    figure.text(0.018, 0.98, "A", va="center", ha="center")
    figure.text(0.018, 0.52, "B", va="center", ha="center")

    figure.align_labels()
    pyplot.savefig("./results/figures/[2-2] reliability detailed.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()


def evaluate_growing(terminal_length, test_times):
    if not path.exists("./results/data/step_2_reliability_extend.pkl"):
        random.seed(2021)
        whole_data = {}
        for length in range(2, terminal_length + 1):
            print("calculate length " + str(length))
            feed_matrix = accessor_to_adjacency_matrix(get_complete_accessor(observed_length=length))
            values, errors, monitor = [], [], Monitor()
            for time in range(test_times):
                matrix = feed_matrix.copy()
                for position in array(list(where(feed_matrix == 1)[:2])).T:
                    if random.random(1)[0] <= 0.5:
                        matrix[position[0], position[1]] = 0
                monitor.output(current_state=time + 1, total_state=test_times)
                calculated_eigenvalue = real(max(linalg.eig(matrix)[0]))
                calculated_value = log(calculated_eigenvalue, 2) if calculated_eigenvalue > 0.0 else 0.0
                approximate_value = approximate_capacity(accessor=adjacency_matrix_to_accessor(matrix), repeats=10)
                values.append(approximate_value)
                errors.append(abs(approximate_value - calculated_value))
            whole_data[length] = (values, errors)

        with open("./results/data/step_2_reliability_extend.pkl", "wb") as file:
            dump(whole_data, file)
    else:
        with open("./results/data/step_2_reliability_extend.pkl", "rb") as file:
            whole_data = load(file)

    pyplot.figure(figsize=(10, 5), tight_layout=True)

    pyplot.fill_between([1.5, 2.5], [-15, -15], [-3, -3], color=colors["diffs"], zorder=0)
    pyplot.fill_between([3.5, 4.5], [-15, -15], [-3, -3], color=colors["diffs"], zorder=0)

    used = []
    for length, data in whole_data.items():
        used.append(log10(data[1]))

    violin = pyplot.violinplot(dataset=used, positions=[0.9, 1.9, 2.9, 3.9, 4.9], showextrema=False, widths=0.35)
    for patch in violin["bodies"]:
        patch.set_edgecolor(colors["trad1"])
        patch.set_facecolor(colors["trad2"])
        patch.set_linewidth(1)
        patch.set_alpha(1)

    for index, data in enumerate(used):
        [v_25, v_50, v_75] = percentile(a=data, q=[25, 50, 75])
        pyplot.scatter(x=[index + 0.9], y=[v_50],
                       color="white", edgecolor=colors["trad1"], linewidth=1, s=25, zorder=4)
        lower, upper = v_25 - 1.5 * (v_75 - v_25), 1.5 * (v_75 - v_25) + v_75
        pyplot.vlines(x=index + 1.2, ymin=lower, ymax=upper,
                      color=colors["trad1"], linewidth=1, zorder=2)
        pyplot.hlines(y=lower, xmin=index + 1.2 - 0.04, xmax=index + 1.2 + 0.04,
                      color=colors["trad1"], linewidth=1, zorder=2)
        pyplot.hlines(y=upper, xmin=index + 1.2 - 0.04, xmax=index + 1.2 + 0.04,
                      color=colors["trad1"], linewidth=1, zorder=2)
        pyplot.fill_between(x=[index + 1.2 - 0.04, index + 1.2 + 0.04], y1=[v_25, v_25], y2=[v_75, v_75],
                            color=colors["trad1"], zorder=3)
        for value in data:
            if value < lower or value > upper:
                pyplot.scatter(x=[index + 1.2], y=[value],
                               color="white", edgecolor=colors["trad1"], linewidth=1, s=20, zorder=4)

    pyplot.xlim(0.5, 5.5)
    pyplot.ylim(-15, -3)
    pyplot.xticks([1, 2, 3, 4, 5],
                  [r"$4^2 \times 4^2$", r"$4^3 \times 4^3$", r"$4^4 \times 4^4$",
                   r"$4^5 \times 4^5$", r"$4^6 \times 4^6$"],
                  fontsize=8)
    pyplot.yticks([-15, -13, -11, -9, -7, -5, -3],
                  ["$10^{-15}$", "$10^{-13}$", "$10^{-11}$", "$10^{-09}$", "$10^{-07}$", "$10^{-05}$", "$10^{-03}$"],
                  fontsize=8)
    pyplot.xlabel("size of adjacency matrix", fontsize=8)
    pyplot.ylabel("relative error with NumPy \"linalg.eig\" function", fontsize=8)
    pyplot.savefig("./results/figures/[2-3] reliability extended.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()
    display_cases()
    evaluate_length_2(test_times=100)
    evaluate_growing(terminal_length=6, test_times=100)
