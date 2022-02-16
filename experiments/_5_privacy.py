from collections import Counter
# noinspection PyPackageRequirements
from matplotlib import pyplot
from numpy import random, load, save, array, zeros, linspace, log, log2, median, sum, min, max, argmax, where
from os import path
from pickle import load as pload
from scipy.optimize import curve_fit
from warnings import filterwarnings

from dsw import Monitor

from experiments import create_folders

filterwarnings("ignore")


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


def draw_total(bit_length):
    def estimated_equation(x, a):
        return a ** x

    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))
    figure = pyplot.figure(figsize=(10, 5), tight_layout=True)

    pyplot.subplot(1, 2, 1)
    reconstructions = load(file="./results/data/step_5_privacy_reconstruction.npy")

    with open("./results/data/step_1_compatibility_capacities.pkl", "rb") as file:
        capacities, _ = pload(file)

    collected_data = zeros(shape=(12, 100))
    for index, data in enumerate(reconstructions):
        number = sum(load(file="./results/data/a" + filter_indices[index] + "[v].npy"))
        for iteration_index, iteration_data in enumerate(data):
            if number > max(iteration_data):
                sequences = linspace(0, 10000, 10001)
                parameter = curve_fit(estimated_equation, xdata=iteration_data, ydata=sequences)[0][0]
                log_sequence_number = number * log(parameter)
            else:
                log_sequence_number = log(argmax(iteration_data))
            value = (log_sequence_number + log(100) + log(capacities[index]) - log(8)) / log(1024.0)
            collected_data[index, iteration_index] = value

    for index in range(12):
        if max(collected_data[index]) - min(collected_data[index]) > 0.05:
            pyplot.hlines(max(collected_data[index]), index - 0.25, index + 0.25, linewidth=1.0)
            pyplot.hlines(min(collected_data[index]), index - 0.25, index + 0.25, linewidth=1.0)
            pyplot.vlines(index, min(collected_data[index]), max(collected_data[index]), linewidth=1.0)
        else:
            pyplot.hlines(median(collected_data[index]), index - 0.25, index + 0.25, linewidth=1.5)

        pyplot.scatter([index], [median(collected_data[index])],
                       color=gradient_colors[index], edgecolor="black", linewidth=1, s=20, zorder=2)

    pyplot.xlabel("constraint set", fontsize=8)
    pyplot.xticks(range(12), filter_indices, fontsize=8)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylabel("transmitted data capacity", fontsize=8)
    pyplot.yticks([0, 1, 2, 3], ["B", "KB", "MB", "GB"], fontsize=8)
    pyplot.ylim(0, 3)

    pyplot.subplot(1, 2, 2)
    follow_ups = zeros(shape=(12, 3))  # [2, 3, 4]
    for index, filter_index in enumerate(filter_indices):
        accessor = load(file="./results/data/a" + filter_index + "[g].npy")
        out_degrees, counts = array(list(Counter(where(accessor >= 0)[0]).values())), zeros(shape=(3,), dtype=float)
        for out_degree in [2, 3, 4]:
            counts[out_degree - 2] = len(where(out_degrees == out_degree)[0])
        counts /= sum(counts)
        follow_ups[index] = counts

    pyplot.hlines(bit_length, 0, 320, color="silver", linewidth=0.75, linestyle="--", zorder=1)
    pyplot.text(4, 256 + 4, "reference", va="bottom", ha="left", fontsize=8)
    for index, filter_index in enumerate(filter_indices):
        lengths = linspace(0, 320, 321)
        rate = follow_ups[index][0] * 2.0 + follow_ups[index][1] * 6.0 + follow_ups[index][2] * 24.0
        values = lengths * log2(rate)
        pyplot.plot(lengths, values, color=gradient_colors[index], alpha=0.5, linewidth=2, zorder=2)

        reference_length = bit_length / log2(rate)
        pyplot.vlines(reference_length, 0, bit_length,
                      color="silver", linewidth=0.75, linestyle="--", zorder=1)
        if index > 0:
            pyplot.scatter(reference_length, bit_length,
                           color=gradient_colors[index], edgecolor="black", s=20, zorder=4,
                           label=" %.2f" % reference_length)
        else:
            pyplot.scatter(reference_length, bit_length,
                           color=gradient_colors[index], edgecolor="black", s=20, zorder=4,
                           label="%.2f" % reference_length)

    pyplot.legend(loc="upper right", ncol=2, fontsize=8)
    pyplot.xlabel("DNA string length", fontsize=8)
    pyplot.xlim(0, 320)
    pyplot.xticks([0, 64, 128, 192, 256, 320], ["0nt", "64nt", "128nt", "192nt", "256nt", "320nt"], fontsize=8)
    pyplot.ylabel("equivalent key strength", fontsize=8)
    pyplot.ylim(0, 384)
    pyplot.yticks([0, 128, 256, 384], [0, 128, 256, 384], fontsize=8)

    figure.align_labels()
    figure.text(0.018, 0.98, "A", va="center", ha="center")
    figure.text(0.500, 0.98, "B", va="center", ha="center")

    pyplot.savefig("./results/figures/[5-1] privacy evaluation.svg",
                   format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


def draw_reason():
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))

    pyplot.figure(figsize=(10, 8))
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
        used_capacities = parameter ** percentages * 100 * capacities[index] / 8
        display_data[index] = used_capacities

    for index, locations in display_data.items():
        pyplot.plot(percentages, locations, color=gradient_colors[index], linewidth=2, zorder=2)

    pyplot.xlabel("graph reconstruction percentage", fontsize=8)
    pyplot.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                  ["0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"], fontsize=8)
    pyplot.xlim(0, 1)
    pyplot.ylabel("transmitted data capacity", fontsize=8)
    pyplot.yticks([0.0 * 1e8, 0.2 * 1e8, 0.4 * 1e8, 0.6 * 1e8, 0.8 * 1e8, 1.0 * 1e8,
                   1.2 * 1e8, 1.4 * 1e8, 1.6 * 1e8, 1.8 * 1e8, 2.0 * 1e8],
                  ["0.0E+8 bytes", "0.2E+8 bytes", "0.4E+8 bytes", "0.6E+8 bytes", "0.8E+8 bytes", "1.0E+8 bytes",
                   "1.2E+8 bytes", "1.4E+8 bytes", "1.6E+8 bytes", "1.8E+8 bytes", "2.0E+8 bytes"],
                  fontsize=8)
    pyplot.ylim(0, 2 * 1e8)

    pyplot.savefig("./results/figures/[5-2] privacy reason.svg",
                   format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()

    # why do we need to pay attention to privacy?
    reconstruction(random_seed=2021, sequencing_length=100, repeats=100)

    # bit length or key strength is 256.
    # British Standards Institution (2020). Cryptographic Mechanisms: Recommendations and Key Lengths.
    draw_total(bit_length=256)
    draw_reason()
