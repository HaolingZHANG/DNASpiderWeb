__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import math
import numpy
from matplotlib import pyplot
from dsw.evaluator import approximate_upper_bound


colors = {
    "default": "#05A5F3",
    "variable": "#F25413",
    "fixed": "#85B803"
}


def draw_graph_information():
    pyplot.figure(figsize=(7, 3.2))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.25)

    # calculate vertex number.
    default, fixed, variable = [], [], []
    for filter_index in ["1", "2", "3", "4", "5", "6"]:
        default_vertices = numpy.load(file="../entities/BM" + str(filter_index) + "[v].npy")
        default.append(math.log(numpy.sum(default_vertices), 4))
        fixed_vertices = numpy.load(file="../entities/FLC" + str(filter_index) + "[v].npy")
        fixed.append(math.log(numpy.sum(fixed_vertices), 4))
        variable_vertices = numpy.load(file="../entities/VLC" + str(filter_index) + "[v].npy")
        variable.append(math.log(numpy.sum(variable_vertices), 4))
    default = numpy.array(default)
    fixed = numpy.array(fixed)
    variable = numpy.array(variable)

    pyplot.subplot(1, 2, 1)
    pyplot.title("vertex number", fontsize=12)
    x = numpy.array(range(6)) + 1
    pyplot.bar(x - 0.2, default - 5, label="BM",
               width=0.2, color=colors["default"], edgecolor="black", linewidth=0.5, zorder=1)
    pyplot.bar(x, fixed - 5, label="FLC",
               width=0.2, color=colors["fixed"], edgecolor="black", linewidth=0.5, zorder=1)
    pyplot.bar(x + 0.2, variable - 5, label="VLC",
               width=0.2, color=colors["variable"], edgecolor="black", linewidth=0.5, zorder=1)
    for value in [0, 1, 2, 3, 4, 5]:
        pyplot.hlines(value, 0.5, 6.5, color="silver", linewidth=0.5, linestyle="--", zorder=0)
    pyplot.xlabel("biochemical requirement index", fontsize=10)
    pyplot.ylabel("log4(vertex number)", fontsize=10)
    pyplot.xticks(x, x, fontsize=8)
    pyplot.yticks([0, 1, 2, 3, 4, 5], [5, 6, 7, 8, 9, 10], fontsize=8)
    pyplot.xlim(0.5, 6.5)
    pyplot.ylim(-0.3, 5.3)

    # calculate out-degree.
    default, fixed, variable = [], [], []
    for filter_index in ["1", "2", "3", "4", "5", "6"]:
        default_vertices = numpy.load(file="../entities/BM" + filter_index + "[v].npy")
        default_transform = numpy.load(file="../entities/BM" + filter_index + "[g].npy")
        default_transform[default_transform >= 0] = 0
        average = (4 ** 11 + numpy.sum(default_transform)) / numpy.sum(default_vertices)
        minimum, maximum = 4, 0
        for vertex in default_transform:
            count = len(numpy.where(vertex >= 0)[0])
            if count > 0:
                minimum = min(minimum, count)
                maximum = max(maximum, count)
        default.append([minimum, average, maximum])
        fixed_vertices = numpy.load(file="../entities/FLC" + filter_index + "[v].npy")
        fixed_transform = numpy.load(file="../entities/FLC" + filter_index + "[g].npy")
        count = len(numpy.where(fixed_transform[numpy.where(fixed_vertices != 0)[0][0]] >= 0)[0])
        stride = math.log(len(fixed_transform[0]), 4)
        if stride > 1:
            average = math.log(count, stride)
        else:
            average = count
        minimum, maximum = average, average
        fixed.append([minimum, average, maximum])
        variable_vertices = numpy.load(file="../entities/VLC" + filter_index + "[v].npy")
        variable_transform = numpy.load(file="../entities/VLC" + filter_index + "[g].npy")
        variable_transform[variable_transform >= 0] = 0
        average = (4 ** 11 + numpy.sum(variable_transform)) / numpy.sum(variable_vertices)
        minimum, maximum = 4, 1
        for vertex in variable_transform:
            count = len(numpy.where(vertex >= 0)[0])
            if count > 0:
                minimum = min(minimum, count)
                maximum = max(maximum, count)
        variable.append([minimum, average, maximum])
    default = numpy.array(default).T
    fixed = numpy.array(fixed).T
    variable = numpy.array(variable).T

    pyplot.subplot(1, 2, 2)
    pyplot.title("composition of out-degree", fontsize=12)
    x = numpy.array(range(6)) + 1
    pyplot.bar(x - 0.2, default[1], width=0.2, color=colors["default"], edgecolor="black", linewidth=0.5, zorder=1)
    for position, minimum, maximum in zip(x, default[0], default[2]):
        pyplot.vlines(position - 0.2, minimum, maximum, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(minimum, position - 0.25, position - 0.15, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(maximum, position - 0.25, position - 0.15, color="black", linewidth=1.5, zorder=2)
    pyplot.bar(x, fixed[1], width=0.2, color=colors["fixed"], edgecolor="black", linewidth=0.5, zorder=1)
    for position, minimum, maximum in zip(x, fixed[0], fixed[2]):
        pyplot.vlines(position, minimum, maximum, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(minimum, position - 0.05, position + 0.05, color="black", linewidth=1.5, zorder=2)
        pyplot.hlines(maximum, position - 0.05, position + 0.05, color="black", linewidth=1.5, zorder=2)
    pyplot.bar(x + 0.2, variable[1], width=0.2, color=colors["variable"], edgecolor="black", linewidth=0.5, zorder=1)
    for position, minimum, maximum in zip(x, variable[0], variable[2]):
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
    pyplot.savefig("../outputs/graph information.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


def draw_fixed_graph_evolution(maximum_stride):
    paths = []
    for index in [1, 2, 3, 4, 5, 6]:
        paths.append(numpy.load(file="../entities/FLC" + str(index) + "[p].npy"))

    pyplot.figure(figsize=(10, 5))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.25, hspace=0.5)

    for path_index, path in enumerate(paths):
        pyplot.subplot(2, 3, path_index + 1)
        pyplot.title("biochemical filter " + str(path_index + 1))
        for value in [1, 1.2, 1.4, 1.6, 1.8, 2.0]:
            pyplot.hlines(value, 0.5, 5.5, color="silver", linewidth=0.5, linestyle="--", zorder=0)
        b_bound = approximate_upper_bound(numpy.load(file="../entities/BM" + str(path_index + 1) + "[g].npy"),
                                          replay=10) * 2.0
        pyplot.hlines(b_bound, 0.5, 5.5, color="black", linewidth=1.5, zorder=0)
        pyplot.text(x=3, y=b_bound + 0.03, s="%.3f" % b_bound, fontsize=9,
                    horizontalalignment="center", verticalalignment="bottom")
        path = numpy.array(path).T
        strides = path[0]
        information_densities = [math.log2(out_degree) / stride for out_degree, stride in zip(path[1], path[0])]
        pyplot.plot(strides, information_densities, marker="o", color=colors["fixed"], ms=5, zorder=2)
        if path_index < 5:
            for stride, information_density in enumerate(information_densities):
                pyplot.vlines(stride + 1, max(information_densities) + 0.1, information_density,
                              color="black", linestyle="dotted", linewidth=1, zorder=1)
                pyplot.text(x=stride + 1, y=max(information_densities) + 0.1, s="%.3f" % information_density,
                            bbox=dict(boxstyle="round", ec="black", fc="white"),
                            ha="center", va="bottom", fontsize=8, zorder=1)
        else:
            for stride, information_density in enumerate(information_densities):
                pyplot.vlines(stride + 1, min(information_densities) - 0.1, information_density,
                              color="black", linestyle="dotted", linewidth=1, zorder=1)
                pyplot.text(x=stride + 1, y=min(information_densities) - 0.1, s="%.3f" % information_density,
                            bbox=dict(boxstyle="round", ec="black", fc="white"),
                            ha="center", va="top", fontsize=8, zorder=1)

        for stride, out_degree in enumerate(path[1]):
            if stride == maximum_stride - 1:
                break
            upper_bound = out_degree + 1
            lower_bound = out_degree
            upper_bound = upper_bound ** ((stride + 2) / (stride + 1))
            if upper_bound == math.floor(upper_bound):
                upper_bound = math.floor(upper_bound) - 1
            else:
                upper_bound = math.floor(upper_bound)

            lower_bound = lower_bound ** ((stride + 2) / (stride + 1))
            if lower_bound == math.ceil(lower_bound):
                lower_bound = math.ceil(lower_bound) + 1
            else:
                lower_bound = math.ceil(lower_bound)
            pyplot.fill_between(x=[stride + 1, stride + 2],
                                y1=[math.log2(out_degree) / (stride + 1), math.log2(upper_bound) / (stride + 2)],
                                y2=[math.log2(out_degree) / (stride + 1), math.log2(lower_bound) / (stride + 2)],
                                color="silver", zorder=0)
            if stride + 2 > len(path[1]):
                x = [stride + 2]
                y = [((math.log2(upper_bound) / (stride + 2)) + (math.log2(lower_bound) / (stride + 2))) / 2]
                pyplot.scatter(x, y, color="red", marker="x", s=15)
        pyplot.xlabel("stride per step", fontsize=10)
        pyplot.ylabel("information density (bit/nt)", fontsize=10)
        pyplot.xticks([1, 2, 3, 4, 5], [1, 2, 3, 4, 5], fontsize=8)
        pyplot.yticks([1, 1.2, 1.4, 1.6, 1.8, 2.0], [1.0, 1.2, 1.4, 1.6, 1.8, 2.0], fontsize=8)
        pyplot.xlim(0.5, 5.5)
        pyplot.ylim(0.85, 2.15)
    pyplot.savefig("../outputs/statistics fixed coders.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


def draw_variable_graph_evolution():
    pyplot.figure(figsize=(10, 5))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.25, hspace=0.5)
    for filter_index in [1, 2, 3, 4, 5, 6]:
        pyplot.subplot(2, 3, filter_index)
        pyplot.title("biochemical filter " + str(filter_index), fontsize=12)
        b_bound = approximate_upper_bound(numpy.load(file="../entities/BM" + str(filter_index) + "[g].npy"),
                                          replay=10)
        v_bound = approximate_upper_bound(numpy.load(file="../entities/VLC" + str(filter_index) + "[g].npy"),
                                          replay=10)
        v_trans = numpy.load("../outputs/VLC-" + str(filter_index) + " transcode.npy")
        pyplot.text(x=0, y=b_bound + 0.01, s="%.3f" % b_bound, fontsize=9,
                    horizontalalignment="center", verticalalignment="bottom")
        pyplot.hlines(b_bound, -0.25, 0.25, color="black", linewidth=1.5)
        pyplot.text(x=1, y=v_bound + 0.01, s="%.3f" % v_bound, fontsize=9,
                    horizontalalignment="center", verticalalignment="bottom")
        pyplot.hlines(v_bound, 0.75, 1.25, color=colors["variable"], linewidth=1.5)

        violin = pyplot.violinplot(dataset=v_trans, positions=[1], showextrema=False)
        for patch in violin["bodies"]:
            patch.set_facecolor(colors["variable"])
            patch.set_edgecolor("black")
            patch.set_linewidth(0.5)
            patch.set_alpha(1)
        min_value, median, max_value = numpy.percentile(v_trans, [0, 50, 100])
        pyplot.scatter([1], median, color=colors["variable"], s=5, linewidths=1, edgecolors="black", zorder=4)

        for value in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            pyplot.hlines(value, -0.5, 1.5, color="silver", linewidth=0.5, linestyle="--", zorder=0)
        pyplot.xlabel("actual transcoding performances", fontsize=10)
        pyplot.ylabel("spectral radius (capacity)", fontsize=10)
        pyplot.xticks(range(2), ["benchmark", "VLC"], fontsize=8)
        pyplot.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0], [0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize=8)
        pyplot.xlim(-0.5, 1.5)
        pyplot.ylim(0.425, 1.075)
    pyplot.savefig("../outputs/statistics variable coders.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()
