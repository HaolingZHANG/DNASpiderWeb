__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from itertools import product
from math import log
from matplotlib import pyplot
from numpy import array, zeros, ones, zeros_like, percentile, sum, min, max, log10, linalg, real, random, where
from os import path
from pickle import load, dump

from dsw import Monitor, calculate_capacity
from dsw.spiderweb import to_adjacency_matrix, obtain_latters, to_graph, obtain_number


colors = {
    "trad1": "#81B8DF",
    "trad2": "#B1CCDF",
}


# noinspection PyUnresolvedReferences
def display_cases():
    def create_matrix(oligos):
        matrix = zeros(shape=(16, 16), dtype=int)
        vertices = [obtain_number(oligo=vertex) for vertex in oligos]

        for former, latter in product(range(len(oligos)), repeat=2):
            if oligos[former][1] == oligos[latter][0]:
                matrix[vertices[former], vertices[latter]] = 1
        return matrix

    cases = [create_matrix(oligos=["AC",
                                   "CG",
                                   "GT",
                                   "TA"]),  # a cycle.
             create_matrix(oligos=["AC", "AG",
                                   "CA", "CG",
                                   "GA", "GT",
                                   "TC", "TG"]),  # GC content is 50%.
             create_matrix(oligos=["AC", "AG", "AT",
                                   "CA", "CG", "CT",
                                   "GA", "GC", "GT",
                                   "TA", "TC", "TG"]),  # no homopolymer.
             create_matrix(oligos=["AA", "AC", "AG", "AT",
                                   "CA", "CC", "CG", "CT",
                                   "GA", "GC", "GG", "GT",
                                   "TA", "TC", "TG", "TT"])]  # complete graph.

    pyplot.figure(figsize=(10, 5))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.15, hspace=0.3)
    for index, case in enumerate(cases):
        pyplot.subplot(2, len(cases), index + 1)
        for x in range(16):
            for y in range(16):
                if case[x, y] == 1:
                    pyplot.fill_between(x=[x, x + 1], y1=[y, y], y2=[y + 1, y + 1],
                                        color=colors["trad1"])
        pyplot.xticks(array(list(range(16))) + 0.5,
                      ["0", "", "", "3", "", "", "6", "", "", "9", "", "", "12", "", "", "15"])
        if index > 0:
            pyplot.yticks(array(list(range(16))) + 0.5,
                          ["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""])
        else:
            pyplot.yticks(array(list(range(16))) + 0.5,
                          ["0", "", "", "3", "", "", "6", "", "", "9", "", "", "12", "", "", "15"])
            pyplot.ylabel("index of latter DNA string")
        for value in range(16):
            pyplot.vlines(value, 0, 16, color="black", linewidth=0.5)
            pyplot.hlines(value, 0, 16, color="black", linewidth=0.5)
        pyplot.xlabel("index of former DNA string")
        pyplot.xlim(0, 16)
        pyplot.ylim(0, 16)
        axes = pyplot.subplot(2, len(cases), index + len(cases) + 1)
        _, process = calculate_capacity(graph=to_graph(case), replay=1, need_process=True)
        if len(process) == 1:
            pyplot.scatter([0], process, color=colors["trad1"])
        pyplot.plot(range(len(process)), process, color=colors["trad1"], marker="o", linewidth=2)
        pyplot.text(len(process) - 0.8, process[-1], "%.3f" % process[-1], va="center", ha="left")
        pyplot.xlim(-0.1, 2.1)
        pyplot.ylim(-0.1, 2.1)
        if index > 0:
            pyplot.yticks([0, 0.5, 1, 1.5, 2], ["", "", "", "", ""])
        else:
            pyplot.ylabel("information density")
        pyplot.xticks([0, 1, 2], [1, 2, 3])
        pyplot.xlabel("iteration")
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)

    pyplot.savefig("../results/matrices.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


def display_length_2_detailed(data, errors):
    errors = log10(errors)

    violin = pyplot.violinplot(dataset=errors, positions=[0], showextrema=False, vert=False, widths=1)
    paint_data = violin["bodies"][0].get_paths()[0].vertices
    paint_data = paint_data[len(paint_data) // 2 + 1: -1]
    pyplot.close()

    x, y = paint_data[:, 0], paint_data[:, 1]
    y = y / sum(y) * 100.0

    pyplot.figure(figsize=(10, 4))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplot(2, 1, 2)

    pyplot.plot(x, y, color=colors["trad1"], linewidth=1.5, zorder=4)
    pyplot.fill_between(x, y, zeros_like(y), color=colors["trad2"], alpha=0.75, zorder=3)

    [v_25, v_50, v_75] = percentile(a=errors, q=[25, 50, 75])

    pyplot.vlines(x=v_50, ymin=0, ymax=6.5, color="black", linewidth=1, zorder=2)
    pyplot.vlines(x=v_50, ymin=6.9, ymax=8.5, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=8.5, xmin=v_50 - 0.3, xmax=v_50 + 0.3, color="black", linewidth=1, zorder=2)
    pyplot.text(x=v_50, y=8.6, s="median", va="bottom", ha="center", zorder=2)

    pyplot.vlines(x=v_25, ymin=0, ymax=7, color="black", linewidth=1, zorder=2)
    pyplot.vlines(x=v_75, ymin=0, ymax=7, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=6.7, xmin=v_25, xmax=v_75, color="black", linewidth=1, zorder=2)
    pyplot.text(x=v_75 + 0.1, y=6.7, s="interquartile range", va="center", ha="left", zorder=2)
    lower, upper = v_25 - 1.5 * (v_75 - v_25), 1.5 * (v_75 - v_25) + v_75
    pyplot.vlines(x=lower, ymin=0, ymax=4, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=3.7, xmin=lower, xmax=min(errors), color="black", linewidth=1, zorder=2)
    pyplot.text(x=min(errors) - 0.1, y=3.7, s="outlier range", va="center", ha="right", zorder=2)
    pyplot.vlines(x=upper, ymin=0, ymax=4, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=3.7, xmin=upper, xmax=max(errors), color="black", linewidth=1, zorder=2)
    pyplot.text(x=max(errors) + 0.1, y=3.7, s="outlier range", va="center", ha="left", zorder=2)

    pyplot.xlim(-15.1, -2.9)
    pyplot.ylim(-0.5, 10.5)
    pyplot.xlabel("systematic error")
    pyplot.ylabel("frequency")
    pyplot.xticks([-15, -13, -11, -9, -7, -5, -3],
                  ["$10^{-15}$", "$10^{-13}$", "$10^{-11}$", "$10^{-09}$", "$10^{-07}$", "$10^{-05}$", "$10^{-03}$"])
    pyplot.yticks([0, 2, 4, 6, 8, 10], ["0%", "2%", "4%", "6%", "8%", "10%"])

    pyplot.subplot(2, 1, 1)
    pyplot.scatter(errors, data, color=colors["trad2"], edgecolor=colors["trad1"], zorder=2)
    pyplot.vlines(x=v_25, ymin=0, ymax=2, color="black", linestyle="--", linewidth=1, zorder=1)
    pyplot.vlines(x=v_50, ymin=0, ymax=2, color="black", linestyle="--", linewidth=1, zorder=1)
    pyplot.vlines(x=v_75, ymin=0, ymax=2, color="black", linestyle="--", linewidth=1, zorder=1)
    pyplot.vlines(x=lower, ymin=0, ymax=2, color="black", linestyle="--", linewidth=1, zorder=1)
    pyplot.vlines(x=upper, ymin=0, ymax=2, color="black", linestyle="--", linewidth=1, zorder=1)

    pyplot.xlim(-15.1, -2.9)
    pyplot.ylim(-0.1, 2.1)
    pyplot.xticks([-15, -13, -11, -9, -7, -5, -3], ["", "", "", "", "", "", ""])
    pyplot.yticks([0, 0.5, 1.0, 1.5, 2.0], ["0.0  ", "0.5  ", "1.0  ", "1.5  ", "2.0  "])
    pyplot.ylabel("information density")

    pyplot.savefig("../results/errors.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


def display_growing(whole_data):
    pyplot.figure(figsize=(10, 3))
    pyplot.rc("font", family="Times New Roman")

    for value in [1.5, 2.5, 3.5, 4.5]:
        pyplot.vlines(x=value, ymin=-15, ymax=-3, color="black", linestyle="--", linewidth=0.5)

    used = []
    for length, data in whole_data.items():
        used.append(log10(data[1]))
    violin = pyplot.violinplot(dataset=used, positions=[0.9, 1.9, 2.9, 3.9, 4.9], showextrema=False, widths=0.35)
    for patch in violin["bodies"]:
        patch.set_edgecolor(colors["trad1"])
        patch.set_facecolor(colors["trad2"])
        patch.set_linewidth(1.5)
        patch.set_alpha(1)

    for index, data in enumerate(used):
        [v_25, v_50, v_75] = percentile(a=data, q=[25, 50, 75])
        pyplot.scatter(x=[index + 0.9], y=[v_50],
                       color="white", edgecolor=colors["trad1"], linewidth=1.5, s=25, zorder=4)
        lower, upper = v_25 - 1.5 * (v_75 - v_25), 1.5 * (v_75 - v_25) + v_75
        pyplot.vlines(x=index + 1.2, ymin=lower, ymax=upper,
                      color=colors["trad1"], linewidth=1.5, zorder=2)
        pyplot.hlines(y=lower, xmin=index + 1.2 - 0.04, xmax=index + 1.2 + 0.04,
                      color=colors["trad1"], linewidth=1.5, zorder=2)
        pyplot.hlines(y=upper, xmin=index + 1.2 - 0.04, xmax=index + 1.2 + 0.04,
                      color=colors["trad1"], linewidth=1.5, zorder=2)
        pyplot.fill_between(x=[index + 1.2 - 0.04, index + 1.2 + 0.04], y1=[v_25, v_25], y2=[v_75, v_75],
                            color=colors["trad1"], zorder=3)
        for value in data:
            if value < lower or value > upper:
                pyplot.scatter(x=[index + 1.2], y=[value],
                               color="white", edgecolor=colors["trad1"], linewidth=1.5, s=20, zorder=4)

    pyplot.xlim(0.5, 5.5)
    pyplot.ylim(-15, -3)
    pyplot.xticks([1, 2, 3, 4, 5],
                  [r"$4^2 \times 4^2$", r"$4^3 \times 4^3$", r"$4^4 \times 4^4$",
                   r"$4^5 \times 4^5$", r"$4^6 \times 4^6$"])
    pyplot.yticks([-15, -13, -11, -9, -7, -5, -3],
                  ["$10^{-15}$", "$10^{-13}$", "$10^{-11}$", "$10^{-09}$", "$10^{-07}$", "$10^{-05}$", "$10^{-03}$"])
    pyplot.xlabel("size of adjacency matrix")
    pyplot.ylabel("systematic error")

    pyplot.savefig("../results/extend_errors.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


def get_feed_matrix(length):
    accessor = -ones(shape=(int(4 ** length), 4), dtype=int)
    for vertex_index in range(int(4 ** length)):
        for position, latter_vertex_index in enumerate(obtain_latters(current=vertex_index, length=length)):
            accessor[vertex_index][position] = latter_vertex_index
    return to_adjacency_matrix(graph=accessor)


# noinspection PyTypeChecker
def evaluate_length_2(test_times):
    random.seed(2021)
    feed_matrix = get_feed_matrix(length=2)
    values, errors, monitor = [], [], Monitor()
    for time in range(test_times):
        matrix = feed_matrix.copy()
        for position in array(list(where(feed_matrix == 1)[:2])).T:
            if random.random(1)[0] <= 0.5:
                matrix[position[0], position[1]] = 0
        monitor.output(current_state=time + 1, total_state=test_times)
        calculated_eigenvalue = real(max(linalg.eig(matrix)[0]))
        calculated_value = log(calculated_eigenvalue, 2) if calculated_eigenvalue > 0.0 else 0.0
        approximate_value = calculate_capacity(graph=to_graph(matrix), replay=10)
        values.append(approximate_value)
        errors.append(abs(approximate_value - calculated_value))

    return values, errors


# noinspection PyTypeChecker
def evaluate_growing(terminal_length, test_times):
    if not path.exists("../results/Practicability.pkl"):
        random.seed(2021)
        whole_data = {}
        for length in range(2, terminal_length + 1):
            print("calculate length " + str(length))
            feed_matrix = get_feed_matrix(length=length)
            values, errors, monitor = [], [], Monitor()
            for time in range(test_times):
                matrix = feed_matrix.copy()
                for position in array(list(where(feed_matrix == 1)[:2])).T:
                    if random.random(1)[0] <= 0.5:
                        matrix[position[0], position[1]] = 0
                monitor.output(current_state=time + 1, total_state=test_times)
                calculated_eigenvalue = real(max(linalg.eig(matrix)[0]))
                calculated_value = log(calculated_eigenvalue, 2) if calculated_eigenvalue > 0.0 else 0.0
                approximate_value = calculate_capacity(graph=to_graph(matrix), replay=10)
                values.append(approximate_value)
                errors.append(abs(approximate_value - calculated_value))
            whole_data[length] = (values, errors)

        with open("../results/Practicability.pkl", "wb") as file:
            dump(whole_data, file)

        return whole_data
    else:
        with open("../results/Practicability.pkl", "rb") as file:
            return load(file)


if __name__ == "__main__":
    display_cases()
    results, differs = evaluate_length_2(test_times=100)
    display_length_2_detailed(data=results, errors=differs)
    data_group = evaluate_growing(terminal_length=6, test_times=100)
    display_growing(whole_data=data_group)
