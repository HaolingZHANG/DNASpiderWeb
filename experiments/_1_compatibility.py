# noinspection PyPackageRequirements
from matplotlib import pyplot
from numpy import array, random, load, save, max, min, median, any, where
from os import path
from pickle import load as pload
from pickle import dump as psave

from dsw import Monitor, encode, decode
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


def display_performances():
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    with open("./results/data/step_1_compatibility_capacities.pkl", "rb") as file:
        capacities, data = pload(file)
    with open("./results/data/step_1_compatibility_code_rates.pkl", "rb") as file:
        transcode_results = pload(file)

    pyplot.figure(figsize=(10, 8))
    for index, filter_index in enumerate(filter_indices):
        capacity, code_rates = capacities[index], transcode_results[index]
        pyplot.text(x=index, y=capacity + 0.02, s="%.4f" % capacity, ha="center", va="bottom", fontsize=8)
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

    pyplot.xlabel("constraint set", fontsize=8)
    pyplot.xticks(range(12), filter_indices, fontsize=8)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylabel("code rate", fontsize=8)
    pyplot.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0],
                  ["1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0"], fontsize=8)
    pyplot.ylim(0.95, 2.05)
    pyplot.savefig("./results/figures/[1-1] compatibility code rates.svg",
                   format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()
    # generate algorithms from 12 sets of local biochemical constraint sets.
    generate_by_filters()
    # evaluate the performances of the generated algorithm.
    evaluate_performances(task_seed=2021, bit_length=100, repeats=100)
    # display the result of performances.
    display_performances()
