from itertools import permutations
# noinspection PyPackageRequirements
from matplotlib import pyplot, rcParams
from numpy import array, random, load, save, max, min, sum, std, median, any, where
from os import path
from pickle import load as pload
from pickle import dump as psave

from dsw import Monitor, encode, decode, bit_to_number, calculus_division
from dsw import find_vertices, connect_valid_graph, connect_coding_graph, approximate_capacity

from experiments import colors, obtain_filters, create_folders


def generate_by_filters():
    filters = obtain_filters()
    for index, bio_filter in filters.items():
        if not path.exists("./results/data/v" + index + "[v].npy"):
            print("Calculate graph " + str(index) + ".")
            vertices = find_vertices(observed_length=10, bio_filter=bio_filter, verbose=True)
            save(file="./results/data/v" + index + "[v].npy", arr=vertices)
            print()
        else:
            vertices = load(file="./results/data/v" + str(index) + "[v].npy")

        if not path.exists("./results/data/v" + index + "[g].npy"):
            accessor = connect_valid_graph(observed_length=10, vertices=vertices, verbose=True)
            save(file="./results/data/v" + index + "[g].npy", arr=accessor)
            print()

        if not (path.exists("./results/data/a" + index + "[v].npy")
                and path.exists("./results/data/a" + index + "[g].npy")):
            vertices, accessor = connect_coding_graph(observed_length=10, vertices=vertices, threshold=2, verbose=True)
            save(file="./results/data/a" + index + "[v].npy", arr=vertices)
            save(file="./results/data/a" + index + "[g].npy", arr=accessor)
            print()


def evaluate_performances(task_seed, bit_length, repeats):
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    if not path.exists("./results/data/step_1_compatibility_capacities.pkl"):
        capacities, data = [], []
        for filter_index in filter_indices:
            accessor = load(file="./results/data/v" + filter_index + "[g].npy")
            capacity, processes = approximate_capacity(accessor=accessor, repeats=10, need_process=True)
            capacities.append(capacity)
            data.append((capacity, processes))
        with open("./results/data/step_1_compatibility_capacities.pkl", "wb") as file:
            psave((capacities, data), file)

    if not path.exists("./results/data/step_1_compatibility_code_rates.pkl"):
        transcode_results, monitor = [], Monitor()
        for filter_index in filter_indices:
            print("Evaluate coding algorithm " + filter_index + ".")
            accessor = load(file="./results/data/a" + filter_index + "[g].npy")
            vertices = where(load(file="./results/data/a" + filter_index + "[v].npy") == 1)[0]

            random.seed(task_seed)

            encoded_matrix, nucleotide_numbers = random.randint(0, 2, size=(repeats, bit_length)), []
            for current, start_index in enumerate(vertices):
                dna_strings, nucleotide_number = [], 0
                for binary_message in encoded_matrix:
                    dna_string = encode(binary_message=binary_message, accessor=accessor, start_index=start_index)
                    nucleotide_number += len(dna_string)
                    dna_strings.append(dna_string)

                nucleotide_numbers.append(nucleotide_number)

                # temp the correctness of transcoding.
                decoded_matrix = []
                for dna_string in dna_strings:
                    binary_message = decode(dna_string=dna_string, bit_length=bit_length,
                                            accessor=accessor, start_index=start_index)
                    decoded_matrix.append(binary_message)
                decoded_matrix = array(decoded_matrix)

                if any(encoded_matrix != decoded_matrix):
                    raise ValueError("Encoded matrix is not equal to decoded matrix!")

                monitor.output(current + 1, len(vertices))

            transcode_results.append((bit_length * repeats) / array(nucleotide_numbers))

        with open("./results/data/step_1_compatibility_code_rates.pkl", "wb") as file:
            psave(transcode_results, file)


def evaluate_additional_performance(task_seed, bit_length, repeats):
    bio_filters = obtain_filters()
    if not path.exists("./results/data/step_1_compatibility_out_degree_1.pkl"):
        transcode_results, monitor = {}, Monitor()
        for filter_index, bio_filter in bio_filters.items():
            accessor = load(file="./results/data/v" + filter_index + "[g].npy")
            out_degrees = sum(where(accessor >= 0, 1, 0), axis=1)
            if len(out_degrees[out_degrees == 1]) == 0:
                continue

            print("Evaluate coding algorithm " + filter_index + " [out-degree 1].")

            random.seed(task_seed)

            encoded_matrix, nucleotide_numbers = random.randint(0, 2, size=(repeats, bit_length)), []
            current = 0
            for mapping in permutations(["A", "C", "G", "T"]):
                for binary_message in encoded_matrix:
                    available, nucleotide_number = mapping, 0
                    quotient, dna_string, monitor = bit_to_number(binary_message), "", Monitor()
                    while quotient != "0":
                        if len(available) > 1:  # current vertex contains information.
                            quotient, remainder = calculus_division(number=quotient, base=str(len(available)))
                            dna_string += available[int(remainder)]

                        elif len(available) == 1:  # current vertex does not contain information.
                            dna_string += available[0]

                        available = []
                        for potential_nucleotide in mapping:
                            if bio_filter.valid(dna_string + potential_nucleotide, only_last=True):
                                available.append(potential_nucleotide)

                        if len(available) == 0:
                            dna_string = ""
                            break

                    if len(dna_string) > 0:
                        nucleotide_numbers.append(len(dna_string))
                    else:
                        nucleotide_numbers.append(-1)

                    monitor.output(current + 1, len(encoded_matrix) * 24)
                    current += 1

            transcode_results[filter_index] = bit_length / array(nucleotide_numbers)

        with open("./results/data/step_1_compatibility_out_degree_1.pkl", "wb") as file:
            psave(transcode_results, file)


def display_difference():
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    rates = []
    for filter_index in filter_indices:
        accessor = load(file="./results/data/v" + filter_index + "[g].npy")
        out_degrees = sum(where(accessor >= 0, 1, 0), axis=1)
        rates.append(len(out_degrees[out_degrees == 1]) / len(out_degrees[out_degrees >= 1]))

    figure = pyplot.figure(figsize=(10, 8), tight_layout=True)
    rcParams["font.family"] = "Linux Libertine"
    pyplot.subplot(2, 1, 1)
    visited_indices = []
    for index in range(12):
        if rates[index] > 0:
            pyplot.bar([len(visited_indices)], [rates[index]], color="silver", edgecolor="black", width=0.6)
            pyplot.text(len(visited_indices), rates[index] + 0.001, "%.1f" % (rates[index] * 100) + "%", fontsize=12,
                        va="bottom", ha="center")
            visited_indices.append(filter_indices[index])

    pyplot.xlabel("selected constraint set", fontsize=14)
    pyplot.xticks(range(len(visited_indices)), visited_indices, fontsize=14)
    pyplot.xlim(-0.35, 3.35)
    pyplot.ylabel("probability of out-degree 1", fontsize=14)
    pyplot.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10], ["0%", "2%", "4%", "6%", "8%", "10%"], fontsize=14)
    pyplot.ylim(-0.005, 0.105)

    with open("./results/data/step_1_compatibility_code_rates.pkl", "rb") as file:
        proposed_results = pload(file)

    with open("./results/data/step_1_compatibility_out_degree_1.pkl", "rb") as file:
        out_1_results = pload(file)

    pyplot.subplot(2, 1, 2)
    visited_indices = []
    for index, (filter_index, out_values) in enumerate(out_1_results.items()):
        visited_indices.append(filter_index)
        value_1, value_2 = std(out_values), std(proposed_results[int(filter_index) - 1])
        if index > 0:
            pyplot.bar([index - 0.15], [value_1], color=colors["trad1"], edgecolor="black", width=0.3)
            pyplot.bar([index + 0.15], [value_2], color=colors["algo1"], edgecolor="black", width=0.3)
        else:
            pyplot.bar([index - 0.15], [value_1],
                       color=colors["trad1"], edgecolor="black", width=0.3, label="retain vertex of out-degree 1")
            pyplot.bar([index + 0.15], [value_2],
                       color=colors["algo1"], edgecolor="black", width=0.3, label="remove vertex of out-degree 1")
        pyplot.text(index - 0.15, value_1 + 0.001, "%.3f" % value_1, fontsize=12, va="bottom", ha="center")
        pyplot.text(index + 0.15, value_2 + 0.001, "%.3f" % value_2, fontsize=12, va="bottom", ha="center")

    pyplot.legend(loc="upper right", ncol=2, fontsize=12)
    pyplot.xlabel("selected constraint set", fontsize=14)
    pyplot.xticks(range(len(visited_indices)), visited_indices, fontsize=14)
    pyplot.xlim(-0.35, 3.35)
    pyplot.ylabel("standard deviation of code rate", fontsize=14)
    pyplot.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10], ["0.00", "0.02", "0.04", "0.06", "0.08", "0.10"], fontsize=14)
    pyplot.ylim(-0.005, 0.105)

    figure.align_labels()
    figure.text(0.024, 1.00, "A", va="center", ha="center", fontsize=14)
    figure.text(0.024, 0.51, "B", va="center", ha="center", fontsize=14)

    pyplot.savefig("./results/figures/Figure S1.pdf",
                   format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.show()


def display_performances():
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    with open("./results/data/step_1_compatibility_capacities.pkl", "rb") as file:
        capacities, data = pload(file)
    with open("./results/data/step_1_compatibility_code_rates.pkl", "rb") as file:
        transcode_results = pload(file)

    pyplot.figure(figsize=(10, 8))
    rcParams["font.family"] = "Linux Libertine"
    for index, filter_index in enumerate(filter_indices):
        capacity, code_rates = capacities[index], transcode_results[index]
        pyplot.text(x=index, y=capacity + 0.01, s="%.4f" % capacity, ha="center", va="bottom", fontsize=12)
        pyplot.hlines(capacity, index - 0.4, index + 0.4, color=colors["trad1"], linewidth=1)
        if max(code_rates) - min(code_rates) < 0.01 or max(code_rates) > capacity:
            code_rate = median(code_rates)
            if code_rate > capacity:
                code_rate = capacity
            pyplot.hlines(code_rate, index - 0.25, index + 0.25, linewidths=1, edgecolors=colors["algo1"], zorder=3)
            pyplot.scatter([index], code_rate, color="white", edgecolor=colors["algo1"], linewidth=1, s=15, zorder=4)
        else:
            violin = pyplot.violinplot(dataset=code_rates, positions=[index], bw_method=0.5, showextrema=False)
            for patch in violin["bodies"]:
                patch.set_edgecolor(colors["algo1"])
                patch.set_facecolor(colors["algo2"])
                patch.set_linewidth(1)
                patch.set_alpha(1)
            pyplot.scatter([index], median(code_rates), color="white", edgecolor=colors["algo1"],
                           linewidth=1, s=15, zorder=4)
        if index % 2 != 0:
            pyplot.fill_between([index - 0.5, index + 0.5], [0.92, 0.92], [2.08, 2.08],
                                color=colors["diffs"], zorder=0)

    pyplot.xlabel("constraint set", fontsize=14)
    pyplot.xticks(range(12), filter_indices, fontsize=14)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylabel("code rate", fontsize=14)
    pyplot.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0],
                  ["1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0"], fontsize=14)
    pyplot.ylim(0.95, 2.05)
    pyplot.savefig("./results/figures/Figure S2.pdf",
                   format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()
    # generate algorithms from 12 sets of local biochemical constraint sets.
    generate_by_filters()
    # evaluate the performances of the generated algorithm.
    evaluate_performances(task_seed=2021, bit_length=100, repeats=100)
    evaluate_additional_performance(task_seed=2021, bit_length=100, repeats=100)
    # display the result of performances.
    display_difference()
    display_performances()
