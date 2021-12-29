__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


# noinspection PyPackageRequirements
from matplotlib import pyplot
from numpy import load, array, zeros, max, median, mean, deg2rad

from experiments import colors


def load_data():
    algorithm_data = {
        "SPIDER-WEB": [],
        "HEDGES": [],
        "DNA Fountain": [],
        "Yin-Yang Code": [],
    }

    results = load(file="./results/data/step_3_stability_evaluation.npy")

    # evaluate code rate. (constraint 9, 10, 11, 12)
    values = [[[] for _ in range(4)] for _ in range(4)]
    for sample in results:
        if sample[1] >= 9:
            values[int(sample[1]) - 9][int(sample[0])].append((8 * 1024 * 64) / sample[3])
    rates = zeros(shape=(4, 4))
    for filter_index in range(4):
        for algorithm_index in range(4):
            if len(values[filter_index][algorithm_index]) > 0:
                rates[filter_index, algorithm_index] = median(array(values[filter_index][algorithm_index]))
    for index in range(len(rates)):
        rates[index] = rates[index] / max(rates[index])

    rates = mean(rates, axis=0)
    algorithm_data["SPIDER-WEB"].append(rates[0])
    # considering the section (Imposing DNA Output Sequence Constraints) of HEDGES,
    # transcoding efficiency of SPIDER-WEB and HEDGES should be equivalent.
    algorithm_data["HEDGES"].append(rates[0])
    algorithm_data["DNA Fountain"].append(rates[1])
    algorithm_data["Yin-Yang Code"].append(rates[2])

    # evaluate stable.
    counts = array([0, 0, 0, 0])
    for sample in results:
        counts[int(sample[0])] += sample[-1] > 0
    counts = counts / max(counts)
    counts[3] = 1
    algorithm_data["SPIDER-WEB"].append(counts[0])
    algorithm_data["HEDGES"].append(counts[1])
    algorithm_data["DNA Fountain"].append(counts[2])
    algorithm_data["Yin-Yang Code"].append(counts[3])

    # evaluate repairability.
    # unknown for DNA Fountain and Yin-Yang Code.
    # Theoretically, both SPIDER-WEB and HEDGES are tree-based searches, so the heap size is equivalent.
    algorithm_data["SPIDER-WEB"].append(1.0)
    algorithm_data["HEDGES"].append(1.0)
    algorithm_data["DNA Fountain"].append(0.0)
    algorithm_data["Yin-Yang Code"].append(0.0)

    return algorithm_data


# noinspection PyUnresolvedReferences
def draw_data(algorithm_data):
    labels = ["code rate", "stability", "repairability"]
    angles = array([30, 150, 270])
    positions = deg2rad(angles)
    x = positions.tolist() + [positions[0]]

    figure = pyplot.figure(figsize=(10, 2.5), tight_layout=True)
    pyplot.subplots_adjust(wspace=0.4, hspace=0.4)
    ax = pyplot.subplot(1, 4, 1, polar=True)
    pyplot.title("SPIDER-WEB", fontsize=8)
    y = algorithm_data["SPIDER-WEB"] + [algorithm_data["SPIDER-WEB"][0]]
    ax.plot(x, y, color=colors["algo1"], linewidth=1.5)
    ax.fill(positions, algorithm_data["SPIDER-WEB"], facecolor=colors["algo2"])
    ax.set_thetagrids(angles, labels, fontsize=8, color="white")
    ax.set_rmax(1)
    ax.set_rticks([])
    ax.set_rlabel_position(0)
    ax.tick_params(grid_color="black")
    ax.grid(True)

    ax = pyplot.subplot(1, 4, 2, polar=True)
    pyplot.title("HEDGES", fontsize=8)
    y = algorithm_data["HEDGES"] + [algorithm_data["HEDGES"][0]]
    ax.plot(x, y, color=colors["hedc1"], linewidth=1.5)
    ax.fill(positions, algorithm_data["HEDGES"], facecolor=colors["hedc2"])
    ax.set_thetagrids(angles, labels, fontsize=8, color="white")
    ax.set_rmax(1)
    ax.set_rticks([])
    ax.set_rlabel_position(0)
    ax.tick_params(grid_color="black")
    ax.grid(True)

    ax = pyplot.subplot(1, 4, 3, polar=True)
    pyplot.title("DNA Fountain", fontsize=8)
    y = algorithm_data["DNA Fountain"] + [algorithm_data["DNA Fountain"][0]]
    ax.plot(x, y, color=colors["foun1"], linewidth=1.5)
    ax.fill(positions, algorithm_data["DNA Fountain"], facecolor=colors["foun2"])
    ax.set_thetagrids(angles, labels, fontsize=8, color="white")
    ax.set_rmax(1)
    ax.set_rticks([])
    ax.set_rlabel_position(0)
    ax.tick_params(grid_color="black")
    ax.grid(True)

    ax = pyplot.subplot(1, 4, 4, polar=True)
    pyplot.title("Yin-Yang Code", fontsize=8)
    y = algorithm_data["Yin-Yang Code"] + [algorithm_data["Yin-Yang Code"][0]]
    ax.plot(x, y, color=colors["yyco1"], linewidth=1.5)
    ax.fill(positions, algorithm_data["Yin-Yang Code"], facecolor=colors["yyco2"])
    ax.set_thetagrids(angles, labels, fontsize=8, color="white")
    ax.set_rmax(1)
    ax.set_rticks([])
    ax.set_rlabel_position(0)
    ax.tick_params(grid_color="black")
    ax.grid(True)

    figure.text(0.095, 0.095, "repairability", fontsize=8)
    figure.text(0.341, 0.095, "repairability", fontsize=8)
    figure.text(0.587, 0.095, "repairability", fontsize=8)
    figure.text(0.834, 0.095, "repairability", fontsize=8)

    x_bias, y_bias = 0, -0.102
    figure.text(0.193 + x_bias, 0.698 + y_bias, "code rate", fontsize=8, rotation=-60)
    figure.text(0.439 + x_bias, 0.698 + y_bias, "code rate", fontsize=8, rotation=-60)
    figure.text(0.686 + x_bias, 0.698 + y_bias, "code rate", fontsize=8, rotation=-60)
    figure.text(0.932 + x_bias, 0.698 + y_bias, "code rate", fontsize=8, rotation=-60)

    x_bias, y_bias = 0.014, -0.089
    figure.text(0.015 + x_bias, 0.698 + y_bias, "stability", fontsize=8, rotation=60)
    figure.text(0.261 + x_bias, 0.698 + y_bias, "stability", fontsize=8, rotation=60)
    figure.text(0.507 + x_bias, 0.698 + y_bias, "stability", fontsize=8, rotation=60)
    figure.text(0.754 + x_bias, 0.698 + y_bias, "stability", fontsize=8, rotation=60)

    pyplot.savefig("./results/figures/[6-0] result overviews.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    data = load_data()
    draw_data(algorithm_data=data)
