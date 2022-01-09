# noinspection PyPackageRequirements
from matplotlib import pyplot
from numpy import random, load, save, array, zeros, linspace, log, min, median, max, sum, where
from os import path
from pickle import load as pload
from scipy.optimize import curve_fit

from dsw import Monitor, approximate_capacity

from experiments import create_folders


def reconstruction(random_seed, sequencing_length, repeats, threshold=10000):
    nucleotides = ["A", "C", "G", "T"]
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    if not path.exists("./results/data/step_5_privacy_reconstruction.npy"):
        reconstructions = []
        for index, filter_index in enumerate(filter_indices):
            if not path.exists("./results/temp/reconstruction" + filter_index + ".npy"):
                accessor = load(file="./results/data/a" + filter_index + "[g].npy")
                vertices = where(load(file="./results/data/a" + filter_index + "[v].npy") == 1)[0]
                records = zeros(shape=(repeats, threshold + 1))
                random.seed(random_seed)
                for repeat in range(repeats):
                    print("Do the reconstruction test for the constraint set " + filter_index
                          + " with sequencing random " + str(sequencing_length) + "bp as repeat " + str(repeat + 1))

                    monitor = Monitor()
                    visited, number, total, used = zeros(shape=(4 ** 10,), dtype=bool), 0, len(vertices), []
                    while sum(visited) < total:
                        vertex_index = random.choice(vertices)
                        dna_string = ""
                        number += 1
                        for _ in range(sequencing_length):
                            used_index = random.choice(where(accessor[vertex_index] >= 0)[0])
                            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
                            visited[vertex_index] = True
                            dna_string += nucleotide
                        used.append(dna_string)

                        records[repeat, number] = sum(visited)
                        monitor.output(current_state=sum(visited), total_state=total,
                                       extra={"number of DNA string": number})

                        if number == threshold:
                            monitor.output(current_state=total, total_state=total,
                                           extra={"number of DNA string": "> " + str(threshold)})
                            break

                save(file="./results/temp/reconstruction" + filter_index + ".npy", arr=records)
                random.seed(None)
            else:
                records = load(file="./results/temp/reconstruction" + filter_index + ".npy")

            reconstructions.append(records.tolist())

        reconstructions = array(reconstructions)
        save(file="./results/data/step_5_privacy_reconstruction.npy", arr=reconstructions)


def shuffle_evaluation():
    if not path.exists("./results/data/step_5_privacy_shuffle_evaluation.npy"):
        follow_ups = zeros(shape=(12,))
        for index, filter_index in enumerate(["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]):
            accessor = load(file="./results/data/a" + filter_index + "[g].npy")
            capacity = approximate_capacity(accessor=accessor, repeats=10)
            average = 2.0 ** capacity
            minimum, maximum = 4, -1

            for vertex in accessor:
                count = len(vertex[vertex >= 0])
                if count > 0:
                    minimum = min(minimum, count)
                    maximum = max(maximum, count)
            follow_ups[index] = [minimum, average, maximum]

        save(file="./results/data/step_5_privacy_shuffle_evaluation.npy", arr=follow_ups)


def draw(bit_length):
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))
    reference_strength = log(2.0 ** bit_length) / log(4)

    figure = pyplot.figure(figsize=(10, 5), tight_layout=True)

    pyplot.subplot(1, 2, 1)
    reconstructions = load(file="./results/data/step_5_privacy_reconstruction.npy")
    with open("./results/data/step_1_compatibility_capacities.pkl", "rb") as file:
        capacities, _ = pload(file)

    display_data, percentages = {}, linspace(0, 1, 101)
    for index, data in enumerate(reconstructions):
        medians = [0]
        for iteration_data in data.T[1:]:
            used_data = iteration_data[iteration_data > 0]
            if len(used_data) == 0:
                break
            medians.append(median(used_data))

        number = len(where(load(file="./results/data/a" + filter_indices[index] + "[v].npy") == 1)[0])
        medians = array(medians) / number

        def estimated_equation(x, a):
            return a ** x

        numbers = array(list(range(len(medians))))
        parameter = curve_fit(estimated_equation, xdata=medians[1:], ydata=numbers[1:], p0=(2,))[0][0]
        used_capacities = parameter ** percentages * 100 * capacities[index] / 8.0
        display_data[index] = used_capacities
    for index, values in display_data.items():
        value, size = values[-1], 0
        sizes = ["B", "KB", "MB", "GB", "TB"]
        while value > 1024:
            size += 1
            value /= 1024.0

        info = "%.2f" % value
        if len(info) < 6:
            info = "  " * (6 - len(info)) + info
        threshold = log(values[-1]) / log(1024.0)
        pyplot.scatter([0.992], [threshold], label="[" + filter_indices[index] + "] " + info + sizes[size],
                       color=gradient_colors[index], edgecolor="black", s=20, zorder=4)
        pyplot.hlines(threshold, 0, 1, color="silver", linestyle="--", linewidth=0.75, zorder=1)

        pyplot.plot(percentages, log(values) / log(1024.0), color=gradient_colors[index],
                    alpha=0.5, linewidth=2, zorder=2)

    pyplot.legend(loc="upper left", ncol=2, fontsize=8)
    pyplot.xlabel("percentage of graph reconstruction", fontsize=8)
    pyplot.xticks([0, 0.25, 0.5, 0.75, 1.0], ["0%", "25%", "50%", "75%", "100%"], fontsize=8)
    pyplot.xlim(0, 1)
    pyplot.ylabel("transmitted data capacity", fontsize=8)
    pyplot.yticks([0, 1, 2, 3, 4], ["B", "KB", "MB", "GB", "TB"], fontsize=8)
    pyplot.ylim(0, 4)

    pyplot.subplot(1, 2, 2)

    follow_ups = load(file="./results/data/step_5_privacy_shuffle_evaluation.npy")
    pyplot.hlines(reference_strength, 0, 320,
                  color="silver", linewidth=0.75, linestyle="--", zorder=1)
    pyplot.text(4, 128 + 2, "reference", va="bottom", ha="left", fontsize=8)
    for index, filter_index in enumerate(filter_indices):
        lengths = linspace(0, 320, 321)
        medians = log(follow_ups[index][1] ** lengths) / log(4)
        pyplot.plot(lengths, medians, color=gradient_colors[index], alpha=0.5, linewidth=2, zorder=2)

        reference_length = log(2.0 ** 256) / log(follow_ups[index][1])
        pyplot.vlines(reference_length, 0, reference_strength,
                      color="silver", linewidth=0.75, linestyle="--", zorder=1)
        pyplot.scatter(reference_length, reference_strength,
                       color=gradient_colors[index], edgecolor="black", s=20, zorder=4,
                       label="[" + filter_index + "] %.2f" % reference_length)

    pyplot.legend(loc="upper left", ncol=2, fontsize=8)
    pyplot.xlabel("DNA string length", fontsize=8)
    pyplot.xlim(0, 320)
    pyplot.xticks([0, 64, 128, 192, 256, 320], [0, 64, 128, 192, 256, 320], fontsize=8)
    pyplot.ylabel("key strength", fontsize=8)
    pyplot.ylim(0, 256)
    pyplot.yticks([0, 64, 128, 192, 256], [0, 128, 256, 384, 512], fontsize=8)

    figure.text(0.019, 0.98, "A", va="center", ha="center")
    figure.text(0.511, 0.98, "B", va="center", ha="center")
    pyplot.savefig("./results/figures/[5-1] privacy results.svg", format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()

    # why do we need to pay attention to privacy?
    reconstruction(random_seed=2021, sequencing_length=100, repeats=100)

    shuffle_evaluation()

    # bit length or key strength is 256
    # |cite| https://web.archive.org/web/20200404041712/https://www.keylength.com/en/4/
    draw(bit_length=256)
