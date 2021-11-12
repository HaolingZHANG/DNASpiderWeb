__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from math import log2, floor, ceil
from matplotlib import pyplot
from numpy import array, load, max, min, median

from dsw import calculate_capacity


colors = {
    "trad1": "#81B8DF",
    "trad2": "#B1CCDF",
    "algo1": "#FE817D",
    "algo2": "#FCBBAE"
}


def calculate_upper_bounds():
    bounds = []
    for bound_index in range(1, 8 + 1):
        bound_graph = load(file="../entities/BM" + str(bound_index) + "[g].npy")
        bounds.append(calculate_capacity(graph=bound_graph, replay=10))
    return bounds


def evaluate_fixed_graph(bounds, maximum_stride):
    paths = []
    for graph_index in range(1, 9):
        paths.append(load(file="../entities/FLC" + str(graph_index) + "[p].npy"))
    figure = pyplot.figure(figsize=(10, 4))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.15, hspace=0.3)
    for path_index, path in enumerate(paths):
        pyplot.subplot(2, 4, path_index + 1)
        pyplot.title("group " + str(path_index + 1), fontsize=11)
        bound = bounds[path_index]
        pyplot.hlines(bound, 0.5, 5.5, color=colors["trad1"], linewidth=2, zorder=0)
        pyplot.text(x=3, y=bound + 0.03, s="%.3f" % bound, fontsize=9,
                    horizontalalignment="center", verticalalignment="bottom")
        path = array(path).T
        strides = path[0]
        information_densities = [log2(out_degree) / stride for out_degree, stride in zip(path[1], path[0])]
        pyplot.plot(strides, information_densities, marker="o", color=colors["algo1"], ms=5, zorder=2)
        for stride, out_degree in enumerate(path[1]):
            if stride == maximum_stride - 1:
                break
            upper_bound = out_degree + 1
            lower_bound = out_degree
            upper_bound = upper_bound ** ((stride + 2) / (stride + 1))
            if upper_bound == floor(upper_bound):
                upper_bound = floor(upper_bound) - 1
            else:
                upper_bound = floor(upper_bound)

            lower_bound = lower_bound ** ((stride + 2) / (stride + 1))
            if lower_bound == ceil(lower_bound):
                lower_bound = ceil(lower_bound) + 1
            else:
                lower_bound = ceil(lower_bound)
            pyplot.fill_between(x=[stride + 1, stride + 2],
                                y1=[log2(out_degree) / (stride + 1), log2(upper_bound) / (stride + 2)],
                                y2=[log2(out_degree) / (stride + 1), log2(lower_bound) / (stride + 2)],
                                color=colors["algo2"], zorder=0)
            if stride + 2 > len(path[1]):
                x = [stride + 2]
                y = [((log2(upper_bound) / (stride + 2)) + (log2(lower_bound) / (stride + 2))) / 2]
                pyplot.scatter(x, y, color=colors["algo1"], marker="x", s=15)
        if path_index >= 4:
            pyplot.xticks([1, 2, 3, 4, 5], [1, 2, 3, 4, 5], fontsize=8)
        else:
            pyplot.xticks([1, 2, 3, 4, 5], ["", "", "", "", ""], fontsize=8)

        if path_index % 4 == 0:
            pyplot.yticks([1, 1.2, 1.4, 1.6, 1.8, 2.0], [1.0, 1.2, 1.4, 1.6, 1.8, 2.0], fontsize=8)
        else:
            pyplot.yticks([1, 1.2, 1.4, 1.6, 1.8, 2.0], ["", "", "", "", "", ""], fontsize=8)
        pyplot.xlim(0.5, 5.5)
        pyplot.ylim(0.85, 2.15)

    figure.text(0.475, 0.03, "stride per step", va="center", fontsize=10)
    figure.text(0.085, 0.5, "information density", va="center", rotation="vertical", fontsize=10)
    pyplot.savefig("../results/fixed_results.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


def evaluate_variable_graph(bounds):
    figure = pyplot.figure(figsize=(10, 4))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplots_adjust(wspace=0.3)
    for filter_index in range(1, 9):
        pyplot.subplot(1, 8, filter_index)
        pyplot.title("group " + str(filter_index), fontsize=11)
        bound = bounds[filter_index - 1]
        transcode_results = load("../results/VLC" + str(filter_index) + "[t].npy")
        pyplot.text(x=0, y=bound + 0.01, s="%.3f" % bound, fontsize=9,
                    horizontalalignment="center", verticalalignment="bottom")
        pyplot.hlines(bound, -0.5, 0.5, color=colors["trad1"], linewidth=2)
        if max(transcode_results) - min(transcode_results) < 0.01:
            transcode_result = median(transcode_results)
            if transcode_result > bound:
                transcode_result = bound
            pyplot.hlines(transcode_result, -0.25, 0.25, linewidths=2, edgecolors=colors["algo1"], zorder=3)
            pyplot.scatter([0], transcode_result, color="white", edgecolor=colors["algo1"], linewidth=2, s=15, zorder=4)
        else:
            violin = pyplot.violinplot(dataset=transcode_results, positions=[0], bw_method=0.5, showextrema=False)
            for patch in violin["bodies"]:
                patch.set_edgecolor(colors["algo1"])
                patch.set_facecolor(colors["algo2"])
                patch.set_linewidth(1.5)
                patch.set_alpha(1)
            pyplot.scatter([0], median(transcode_results), color="white", edgecolor=colors["algo1"],
                           linewidth=1.5, s=15, zorder=4)
        pyplot.xticks([])
        if filter_index > 1:
            pyplot.yticks([1, 1.2, 1.4, 1.6, 1.8, 2.0], ["", "", "", "", "", ""], fontsize=8)
        else:
            pyplot.yticks([1, 1.2, 1.4, 1.6, 1.8, 2.0], [1.0, 1.2, 1.4, 1.6, 1.8, 2.0], fontsize=8)
        pyplot.xlim(-0.5, 0.5)
        pyplot.ylim(0.92, 2.08)

    figure.text(0.46, 0.075, "actual performance", va="center", fontsize=10)
    figure.text(0.085, 0.5, "information density", va="center", rotation="vertical", fontsize=10)
    pyplot.savefig("../results/variable_results.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


if __name__ == "__main__":
    b_bounds = calculate_upper_bounds()  # [1.000, 1.585, 1.630, 1.670, 1.776, 1.796, 1.815, 1.982]
    for index, b_bound in enumerate(b_bounds):
        print("biochemical constraint group " + str(index + 1) + ": %.3f" % b_bound)

    evaluate_fixed_graph(bounds=b_bounds, maximum_stride=5)
    evaluate_variable_graph(bounds=b_bounds)
