# noinspection PyPackageRequirements
from matplotlib import pyplot
from numpy import load, array, zeros, median, mean, max, std, deg2rad

from experiments import colors


def load_data():
    algorithm_data = {
        "SPIDER-WEB": [],
        "HEDGES": [],
        "DNA Fountain": [],
        "Yin-Yang Code": [],
    }

    results = load(file="./results/data/step_3_stability_evaluation.npy")

    # evaluate code rate.
    values = [[[] for _ in range(12)] for _ in range(4)]
    for sample in results:
        if sample[2] > 0:
            values[int(sample[0])][int(sample[1]) - 1].append((8 * 1024 * 64) / sample[2])

    code_rates = zeros(shape=(4, 12))
    for algorithm_index in range(4):
        for filter_index in range(12):
            if len(values[algorithm_index][filter_index]) > 0:
                code_rates[algorithm_index, filter_index] = median(array(values[algorithm_index][filter_index]))

    for filter_index in range(12):
        code_rates[:, filter_index] /= max(code_rates[:, filter_index])

    scores = array([mean(code_rates[0][code_rates[0] > 0]), mean(code_rates[1][code_rates[1] > 0]),
                    mean(code_rates[2][code_rates[2] > 0]), mean(code_rates[3][code_rates[3] > 0])])
    scores /= max(scores)
    algorithm_data["SPIDER-WEB"].append(scores[0])
    algorithm_data["HEDGES"].append(scores[1])
    algorithm_data["DNA Fountain"].append(scores[2])
    algorithm_data["Yin-Yang Code"].append(scores[3])

    # evaluate the codability.
    counts, total_counts = array([0, 0, 0, 0]), array([0, 0, 0, 0])
    for sample in results:
        counts[int(sample[0])] += sample[-1] > 0
        total_counts[int(sample[0])] += 1

    codabilities = counts / total_counts
    algorithm_data["SPIDER-WEB"].append(codabilities[0])
    algorithm_data["HEDGES"].append(codabilities[1])
    algorithm_data["DNA Fountain"].append(codabilities[2])
    algorithm_data["Yin-Yang Code"].append(codabilities[3])

    # evaluate the parameter sensibility.
    sensibilities = zeros(shape=(4, 12))
    for algorithm_index in range(4):
        for filter_index in range(12):
            if len(values[algorithm_index][filter_index]) > 0:
                sensibilities[algorithm_index, filter_index] = std(array(values[algorithm_index][filter_index]))

    scores = array([mean(sensibilities[0][sensibilities[0] > 0]), mean(sensibilities[1][sensibilities[1] > 0]),
                    mean(sensibilities[2][sensibilities[2] > 0]), mean(sensibilities[3][sensibilities[3] > 0])])
    scores = 1.0 / scores
    scores /= max(scores)

    algorithm_data["SPIDER-WEB"].append(scores[0])
    algorithm_data["HEDGES"].append(scores[1])
    algorithm_data["DNA Fountain"].append(scores[2])
    algorithm_data["Yin-Yang Code"].append(scores[3])

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
    labels = ["[C-R]", "[C-A]", "[P-S]", "[R-A]"]
    angles = array([45, 135, 225, 315])
    positions = deg2rad(angles)
    x = positions.tolist() + [positions[0]]

    pyplot.figure(figsize=(10, 2.5), tight_layout=True)

    pyplot.subplots_adjust(wspace=0.4, hspace=0.4)
    ax = pyplot.subplot(1, 4, 1, polar=True)
    pyplot.title("SPIDER-WEB", fontsize=8)
    y = algorithm_data["SPIDER-WEB"] + [algorithm_data["SPIDER-WEB"][0]]
    ax.plot(x, y, color=colors["algo1"], linewidth=1.5)
    ax.fill(positions, algorithm_data["SPIDER-WEB"], facecolor=colors["algo2"])
    ax.set_thetagrids(angles, labels, fontsize=8)
    ax.set_rmax(1)
    ax.set_rticks([])
    ax.set_rlabel_position(0)
    ax.tick_params(grid_color="silver")
    ax.grid(True)

    ax = pyplot.subplot(1, 4, 2, polar=True)
    pyplot.title("HEDGES", fontsize=8)
    y = algorithm_data["HEDGES"] + [algorithm_data["HEDGES"][0]]
    ax.plot(x, y, color=colors["hedc1"], linewidth=1.5)
    ax.fill(positions, algorithm_data["HEDGES"], facecolor=colors["hedc2"])
    ax.set_thetagrids(angles, labels, fontsize=8)
    ax.set_rmax(1)
    ax.set_rticks([])
    ax.set_rlabel_position(0)
    ax.tick_params(grid_color="silver")
    ax.grid(True)

    ax = pyplot.subplot(1, 4, 3, polar=True)
    pyplot.title("DNA Fountain", fontsize=8)
    y = algorithm_data["DNA Fountain"] + [algorithm_data["DNA Fountain"][0]]
    ax.plot(x, y, color=colors["foun1"], linewidth=1.5)
    ax.fill(positions, algorithm_data["DNA Fountain"], facecolor=colors["foun2"])
    ax.set_thetagrids(angles, labels, fontsize=8)
    ax.set_rmax(1)
    ax.set_rticks([])
    ax.set_rlabel_position(0)
    ax.tick_params(grid_color="silver")
    ax.grid(True)

    ax = pyplot.subplot(1, 4, 4, polar=True)
    pyplot.title("Yin-Yang Code", fontsize=8)
    y = algorithm_data["Yin-Yang Code"] + [algorithm_data["Yin-Yang Code"][0]]
    ax.plot(x, y, color=colors["yyco1"], linewidth=1.5)
    ax.fill(positions, algorithm_data["Yin-Yang Code"], facecolor=colors["yyco2"])
    ax.set_thetagrids(angles, labels, fontsize=8)
    ax.set_rmax(1)
    ax.set_rticks([])
    ax.set_rlabel_position(0)
    ax.tick_params(grid_color="silver")
    ax.grid(True)

    pyplot.savefig("./results/figures/[6-0] result overviews.svg", format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    data = load_data()
    draw_data(algorithm_data=data)
