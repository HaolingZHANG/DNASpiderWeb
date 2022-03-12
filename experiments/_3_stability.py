from itertools import permutations
# noinspection PyPackageRequirements
from matplotlib import pyplot, patches, rcParams
from numpy import random, array, load, save, sum, median, min, max, where
from os import path

from experiments import colors, create_folders, obtain_filters
from experiments.code_encode import generated, hedges, fountain, yinyang, generate, insert_index


def trans_fountain(dataset, shape, index_length, param_number):
    print("Calculate DNA Fountain.")
    header_size = index_length // 8
    parameter_set = random.uniform(low=array([0.1, 0.01]), high=array([1, 0.1]), size=(param_number, 2))
    records, current, total = [], 0, 12 * len(dataset) * param_number

    for filter_index, bio_filter in obtain_filters().items():
        for data_index, data in enumerate(dataset):
            for delta, c in parameter_set:
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                nucleotide_number = fountain(binary_messages=data.reshape(shape).tolist(), bio_filter=bio_filter,
                                             header_size=header_size, c=c, delta=delta)
                records.append([0, int(filter_index), nucleotide_number])
                current += 1

    return records


def trans_yinyang(dataset, shape, index_length, random_seed, param_number):
    print("Calculate Yin-Yang Code.")
    random.seed(random_seed)
    total_indices = array(list(range(6144)))
    random.shuffle(total_indices)
    chosen_indices = total_indices[:param_number]
    records, current, total = [], 0, 12 * len(dataset) * len(chosen_indices)

    for filter_index, bio_filter in obtain_filters().items():
        for data_index, data in enumerate(dataset):
            matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
            for rule_index in chosen_indices:
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                random.seed(random_seed)
                nucleotide_number = yinyang(binary_messages=matrix.tolist(), bio_filter=bio_filter,
                                            rule_index=rule_index)
                records.append([1, int(filter_index), nucleotide_number])
                current += 1

    return records


def trans_hedges(dataset, shape, index_length):
    print("Calculate HEDGES Code.")
    records, current, total = [], 0, 12 * len(dataset) * 24

    for filter_index, bio_filter in obtain_filters().items():
        for data_index, data in enumerate(dataset):
            matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
            for mapping in permutations(["A", "C", "G", "T"]):
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                nucleotide_number = hedges(binary_messages=matrix, mapping=mapping, bio_filter=bio_filter)
                records.append([2, int(filter_index), nucleotide_number])
                current += 1

    return records


def trans_generated(dataset, shape, index_length, random_seed, param_number):
    print("Calculate generated algorithms.")
    records, current, total = [], 0, 12 * len(dataset) * param_number  # number of vertex.

    for filter_index in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
        accessor = load(file="./results/data/a" + filter_index + "[g].npy")
        vertices = where(load(file="./results/data/a" + filter_index + "[v].npy") == 1)[0]
        random.seed(random_seed)
        random.shuffle(vertices)
        chosen_indices = vertices[:param_number]
        for data_index, data in enumerate(dataset):
            matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
            for start_index in chosen_indices:
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                nucleotide_number = generated(binary_messages=matrix, accessor=accessor, start_index=start_index)
                records.append([3, int(filter_index), nucleotide_number])
                current += 1

    return records


def evaluate(dataset):
    if not path.exists("./results/data/step_3_stability_evaluation.npy"):
        records = []
        records += trans_fountain(dataset=dataset, shape=(2048, 256), index_length=32,
                                  param_number=100)
        records += trans_yinyang(dataset=dataset, shape=(4096, 128), index_length=16,
                                 random_seed=2021, param_number=100)
        records += trans_hedges(dataset=dataset, shape=(4096, 128), index_length=16)
        records += trans_generated(dataset=dataset, shape=(4096, 128), index_length=16,
                                   random_seed=2021, param_number=100)
        save(file="./results/data/step_3_stability_evaluation.npy", arr=array(records))


def draw():
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    display_data = [[[] for _ in range(12)] for _ in range(4)]
    for sample in load(file="./results/data/step_3_stability_evaluation.npy"):
        if sample[2] > 0:
            display_data[int(sample[0])][int(sample[1]) - 1].append((8 * 1024 * 64) / sample[2])
        else:
            display_data[int(sample[0])][int(sample[1]) - 1].append(-1)

    pyplot.figure(figsize=(10, 5), tight_layout=True)
    rcParams["font.family"] = "Linux Libertine"
    bias = [-0.375, -0.125, 0.125, 0.375]
    used_colors = [[colors["foun1"], colors["foun2"]], [colors["yyco1"], colors["yyco2"]],
                   [colors["hedc1"], colors["hedc2"]], [colors["algo1"], colors["algo2"]]]
    labels = ["DNA Fountain", "Yin-Yang Code", "HEDGES", "SPIDER-WEB"]
    for bias_index, method_data in enumerate(display_data):
        for location, data in enumerate(method_data):
            if sum(data) == -len(data):
                pyplot.scatter([location + bias[bias_index]], [0], color=used_colors[bias_index][0], marker="x", s=15)
            else:
                if max(data) - min(data) < 0.02:
                    result = median(data)
                    pyplot.hlines(result, location + bias[bias_index] - 0.08, location + bias[bias_index] + 0.08,
                                  linewidths=1, edgecolors=used_colors[bias_index][0], zorder=3)
                    pyplot.scatter([location + bias[bias_index]], result,
                                   color="white", edgecolor=used_colors[bias_index][0], linewidth=1, s=5, zorder=4)
                else:
                    violin = pyplot.violinplot(dataset=data, positions=[location + bias[bias_index]],
                                               bw_method=0.5, showextrema=False, widths=0.16)
                    for patch in violin["bodies"]:
                        patch.set_edgecolor(used_colors[bias_index][0])
                        patch.set_facecolor(used_colors[bias_index][1])
                        patch.set_linewidth(1)
                        patch.set_alpha(1)
                    pyplot.scatter([location + bias[bias_index]], median(data),
                                   color="white", edgecolor=used_colors[bias_index][0], linewidth=1, s=5, zorder=4)

            if location % 2 != 0:
                pyplot.fill_between([location - 0.5, location + 0.5], [-0.08, -0.08], [2.08, 2.08],
                                    color=colors["diffs"], zorder=0)

    legends = [patches.Patch(facecolor=used_colors[index][1], edgecolor=used_colors[index][0],
                             linewidth=1, label=labels[index]) for index in list(range(4))]
    pyplot.legend(handles=legends, loc="upper right", ncol=4, fontsize=12)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylim(-0.08, 2.08)
    pyplot.xticks(range(12), filter_indices, fontsize=14)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
                  ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0"], fontsize=14)
    pyplot.xlabel("constraint set", fontsize=14)
    pyplot.ylabel("code rate", fontsize=14)
    pyplot.savefig("./results/figures/Figure M4.pdf",
                   format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()

    whole_data = generate(random_seed=2021, total_length=8 * 1024 * 64, times=100)
    evaluate(dataset=whole_data)

    draw()
