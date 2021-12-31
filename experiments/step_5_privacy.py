__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


# noinspection PyPackageRequirements
from matplotlib import pyplot
from numpy import random, load, save, array, zeros, zeros_like, linspace, log, min, median, max, sum, where
from os import path, remove
from pickle import load as pload
from pickle import dump as psave
from scipy.special import comb
from scipy.optimize import curve_fit

from dsw import Monitor, approximate_capacity

from experiments import create_folders


def reconstruction(random_seed, sequencing_length, repeats, threshold=10000):
    nucleotides = ["A", "C", "G", "T"]
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    if not path.exists("./results/data/step_5_security_reconstruction.npy"):
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
        save(file="./results/data/step_5_security_reconstruction.npy", arr=reconstructions)

        for filter_index in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
            remove(path="./results/temp/reconstruction" + filter_index + ".npy")


def shuffle_evaluation():
    if not path.exists("./results/data/step_5_security_shuffle_evaluation.npy"):
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

        save(file="./results/data/step_5_security_shuffle_evaluation.npy", arr=follow_ups)


def remove_evaluation(maximum_count, maximum_remove):
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    if not path.exists("./results/data/step_5_security_remove_evaluation.pkl"):
        remove_data = {}
        for filter_index in filter_indices:
            if not path.exists("./results/temp/remove" + filter_index + ".pkl"):
                accessor = load(file="./results/data/a" + filter_index + "[g].npy")
                vertices = where(load(file="./results/data/a" + filter_index + "[v].npy") == 1)[0]
                outs = zeros_like(vertices)
                for index, vertex in enumerate(vertices):
                    outs[index] = len(where(accessor[vertex] >= 0)[0])

                print("Do the remove test for the constraint set " + filter_index)
                usage_3, usage_4 = len(outs[outs == 3]), len(outs[outs == 4]) * 2
                print("There are " + str(usage_3) + " vertices with out-degree is 3 and "
                      + str(usage_4) + " vertices with out-degree is 4.")
                records = [[] for _ in range(maximum_remove + 1)]
                for count_3_1 in range(maximum_count + 1):
                    assign = [0, 0, 0]  # 3/1, 4/1, 4/2
                    if count_3_1 <= usage_3:
                        assign[0] = count_3_1
                        for count_4_1 in range(maximum_count + 1):
                            if count_4_1 <= usage_4:
                                assign[1] = count_4_1
                                for count_4_2 in range(maximum_count + 1):
                                    if count_4_1 + count_4_2 * 2 <= usage_4:
                                        assign[2] = count_4_2

                                        if count_3_1 + count_4_1 + 2 * count_4_2 <= maximum_remove:
                                            strength = 0.0

                                            # delete directed edges from vertex with 3 out-degrees
                                            if assign[0] > 0:
                                                strength += log(comb(N=usage_3, k=assign[0], exact=False)) / log(2.0)
                                                strength += assign[0] * (log(3.0) / log(2.0))

                                            # delete directed edges from vertex with 4 out-degrees
                                            if assign[1] > 0 and assign[2] > 0:
                                                strength += log(comb(N=usage_4, k=assign[1] + assign[2],
                                                                     exact=False)) / log(2.0)
                                                strength += log(comb(N=assign[1] + assign[2], k=assign[1],
                                                                     exact=False)) / log(2.0)
                                                strength += assign[1] * (log(4.0) / log(2.0))
                                                strength += assign[2] * (log(6.0) / log(2.0))
                                            elif assign[1] > 0 and assign[2] == 0:
                                                strength += log(comb(N=usage_4, k=assign[1],
                                                                     exact=False)) / log(2.0)
                                                strength += assign[1] * (log(4.0) / log(2.0))
                                            elif assign[1] == 0 and assign[2] > 0:
                                                strength += log(comb(N=usage_4, k=assign[2],
                                                                     exact=False)) / log(2.0)
                                                strength += assign[2] * (log(6.0) / log(2.0))

                                            remove_number = count_3_1 + count_4_1 + 2 * count_4_2
                                            percentage = comb(N=usage_3, k=assign[0], exact=False)
                                            percentage /= comb(N=usage_3 + usage_4, k=remove_number, exact=False)
                                            percentage *= comb(N=usage_4, k=assign[1] + assign[2], exact=False)

                                            records[int(remove_number)].append((strength, percentage))

                with open("./results/temp/remove" + filter_index + ".pkl", "wb") as file:
                    psave(records, file)
            else:
                with open("./results/temp/remove" + filter_index + ".pkl", "rb") as file:
                    records = pload(file)

            remove_data[int(filter_index)] = records

        with open("./results/data/step_5_security_remove_evaluation.pkl", "wb") as file:
            psave(remove_data, file)

        for filter_index in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
            remove(path="./results/temp/remove" + filter_index + ".pkl")


def draw(bit_length):
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))

    reconstructions = load(file="./results/data/step_5_security_reconstruction.npy")
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

    pyplot.figure(figsize=(10, 5))
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
        pyplot.scatter([0.995], [threshold], label="[" + filter_indices[index] + "] " + info + sizes[size],
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
    pyplot.savefig("./results/figures/[5-1] privacy reconstruction.svg", format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()

    reference_strength = log(2.0 ** bit_length) / log(4)

    figure = pyplot.figure(figsize=(10, 5))
    pyplot.subplots_adjust(wspace=0.1)

    pyplot.subplot(1, 2, 1)
    follow_ups = load(file="./results/data/step_5_security_shuffle_evaluation.npy")
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
    pyplot.xlabel("DNA length", fontsize=8)
    pyplot.xlim(0, 320)
    pyplot.xticks([0, 64, 128, 192, 256, 320], [0, 64, 128, 192, 256, 320], fontsize=8)
    pyplot.ylabel("key strength", fontsize=8)
    pyplot.ylim(0, 256)
    pyplot.yticks([0, 64, 128, 192, 256], [0, 128, 256, 384, 512], fontsize=8)

    pyplot.subplot(1, 2, 2)
    with open("./results/data/step_5_security_remove_evaluation.pkl", "rb") as file:
        remove_data = pload(file)

    pyplot.text(0.5, 256 + 4, "reference", va="bottom", ha="left", fontsize=8)
    pyplot.hlines(bit_length, 0, 48,
                  color="silver", linewidth=0.75, linestyle="--", zorder=1)
    for index, records in remove_data.items():
        strengths = []
        for record in records:
            if len(record) > 0:
                average, total_count = 0, 0
                for strength, percentage in record:
                    average += strength * percentage
                    total_count += percentage
                average /= total_count
                strengths.append(average)
            else:
                strengths.append(0)

        if index > 1:
            pyplot.plot(list(range(33)), strengths[:33],
                        color=gradient_colors[index - 1], linewidth=2, alpha=0.5, zorder=2)
        else:
            pyplot.plot(list(range(33)), zeros(shape=(33,)) + 2,
                        color=gradient_colors[index - 1], linewidth=2, alpha=0.5, zorder=2)

        reference_length = 0
        for location, (former, latter) in enumerate(zip(strengths[:-1], strengths[1:])):
            if former <= bit_length <= latter:
                reference_length = location
                reference_length += (bit_length - former) / (latter - former)
                break

        if reference_length > 0:
            pyplot.scatter(reference_length, bit_length, color=gradient_colors[index - 1], edgecolor="black", s=20,
                           zorder=4, label="[" + filter_indices[index - 1] + "] %.2f" % reference_length)
            pyplot.vlines(reference_length, 0, bit_length,
                          color="silver", linewidth=0.75, linestyle="--", zorder=1)
        else:
            pyplot.scatter(33, bit_length, color=gradient_colors[index - 1], edgecolor="black", s=20, zorder=4,
                           label="[" + filter_indices[index - 1] + "] N/A")

    pyplot.legend(loc="upper left", ncol=2, fontsize=8)
    pyplot.xlabel("removed edge number", fontsize=8)
    pyplot.xticks([0, 8, 16, 24, 32], [0, 8, 16, 24, 32], fontsize=8)
    pyplot.xlim(0, 32)
    pyplot.yticks([0, 128, 256, 384, 512], ["", "", "", "", ""], fontsize=8)
    pyplot.ylim(0, 512)

    figure.text(0.080, 0.9, "A", va="center", ha="center")
    figure.text(0.515, 0.9, "B", va="center", ha="center")
    pyplot.savefig("./results/figures/[5-2] privacy optional strategy.svg", format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()

    # why do we need to pay attention to privacy?
    reconstruction(random_seed=2021, sequencing_length=100, repeats=100)

    shuffle_evaluation()
    remove_evaluation(maximum_count=64, maximum_remove=64)  # 64 is enough

    # bit length or key strength is 256
    # |cite| https://web.archive.org/web/20200404041712/https://www.keylength.com/en/4/
    draw(bit_length=256)
