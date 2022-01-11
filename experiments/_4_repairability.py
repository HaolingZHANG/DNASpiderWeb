# noinspection PyPackageRequirements
from matplotlib import pyplot, patches
from numpy import load, save, ones, array, linspace, concatenate, median, mean, sum, random, where
from os import path
from pickle import load as pload
from pickle import dump as psave

from experiments import colors, create_folders
from experiments.code_repair import show_single_examples, show_multiple_examples
from experiments.code_repair import evaluate_single_error, evaluate_repair_multiple_errors


def single_evaluation(task_seed, repeats, vertex_number):
    if not path.exists("./results/data/step_4_repairability_single_error.npy"):
        filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
        whole_data, observed_length = [], 10
        random.seed(task_seed)
        for filter_index in filter_indices:
            if not path.exists("./results/temp/single" + filter_index + ".npy"):
                print("Analyze the generated algorithm under constraint set " + filter_index)
                accessor = load(file="./results/data/a" + filter_index + "[g].npy")
                vertices = where(load(file="./results/data/a" + filter_index + "[v].npy") == 1)[0]
                random.shuffle(vertices)
                info = evaluate_single_error(accessor=accessor, start_indices=vertices[:vertex_number],
                                             dna_length=30, vt_length=10, repeats=repeats)
                records = concatenate((ones(shape=(1, len(info))).T * int(filter_index), info), axis=1).astype(int)
                save(file="./results/temp/single" + filter_index + ".npy", arr=records)
            else:
                records = load(file="./results/temp/single" + filter_index + ".npy")

            whole_data += records.tolist()

        save(file="./results/data/step_4_repairability_single_error.npy", arr=array(whole_data))


def multiple_evaluation(task_seed, repeats):
    accessor, vertices = load("./results/data/a01[g].npy"), where(load("./results/data/a01[v].npy") == 1)[0]
    dna_lengths = linspace(start=100, stop=800, num=8, dtype=int)
    for error_time in range(1, 11):
        for dna_length in dna_lengths:
            save_path = "./results/temp/multiple" + str(dna_length).zfill(4) + "." + str(error_time).zfill(2) + ".pkl"
            if not path.exists(save_path):
                records = evaluate_repair_multiple_errors(random_seed=task_seed, accessor=accessor, vertices=vertices,
                                                          observed_length=10, vt_length=10,
                                                          repeats=repeats, dna_length=dna_length,
                                                          error_times=error_time, check_iterations=error_time + 1)
                with open(save_path, "wb") as file:
                    psave(records, file)


def draw_normal_situation(vertex_number):
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

    data = load(file="./results/data/step_4_repairability_single_error.npy")
    results = [[[[], [], []] for _ in range(len(filter_indices))] for _ in range(4)]
    for index in range(0, len(data) - vertex_number * 8 + 1, vertex_number * 8):
        sub_data = data[index: index + vertex_number * 8]
        filter_index = sub_data[0, 0]
        for mark, time in zip([0, 1, 2], [3, 4, 1]):  # mark 0 = substitution, 1 = insertion, 2 = deletion
            locations = where(sub_data[:, 2] == mark)[0]
            detect_rate_1 = sum(sub_data[locations][:, 3]) / (vertex_number * time)
            strategy_numbers = sub_data[locations][:, 4]
            strategy_numbers[strategy_numbers == 0] = 1
            detect_rate_2 = sum(sub_data[locations][:, 5]) / (vertex_number * time)
            repair_rate = sum(sub_data[locations][:, 6]) / (vertex_number * time)
            results[0][filter_index - 1][mark].append(detect_rate_1)
            results[1][filter_index - 1][mark].append(detect_rate_2)
            if len(strategy_numbers) > 0:
                results[2][filter_index - 1][mark] += strategy_numbers.tolist()
            results[3][filter_index - 1][mark].append(repair_rate)

    bias = [-0.3, 0, 0.3]
    used_colors = [[colors["algo1"], colors["algo2"]],
                   [colors["trad1"], colors["trad2"]],
                   [colors["foun1"], colors["foun2"]]]
    labels = ["substitution", "insertion", "deletion"]

    figure = pyplot.figure(figsize=(10, 5), tight_layout=True)
    pyplot.subplot(2, 1, 1)
    for filter_index in range(len(filter_indices)):
        for index, values in enumerate(results[0][filter_index]):
            if sum(values) == 0:
                pyplot.scatter([filter_index + bias[index]], [0],
                               color=used_colors[index][0], marker="x", s=15)
            elif sum(values) == len(values):
                pyplot.scatter([filter_index + bias[index]], [1],
                               color=used_colors[index][0], marker="^", s=15)
            elif max(values) - min(values) < 0.02:
                value = median(values)
                pyplot.hlines(value, filter_index + bias[index] - 0.1, filter_index + bias[index] + 0.1,
                              linewidths=1, edgecolors=used_colors[index][0], zorder=3)
                pyplot.scatter([filter_index + bias[index]], value, color="white", edgecolor=used_colors[index][0],
                               linewidth=1, s=8, zorder=4)
            else:
                violin = pyplot.violinplot(dataset=values, positions=[filter_index + bias[index]],
                                           widths=0.2, bw_method=0.5, showextrema=False)
                for patch in violin["bodies"]:
                    patch.set_edgecolor(used_colors[index][0])
                    patch.set_facecolor(used_colors[index][1])
                    patch.set_linewidth(1)
                    patch.set_alpha(1)
                pyplot.scatter([filter_index + bias[index]], median(values),
                               color="white", edgecolor=used_colors[index][0],
                               linewidth=1, s=8, zorder=4)

        if filter_index % 2 != 0:
            pyplot.fill_between([filter_index - 0.5, filter_index + 0.5], [-0.05, -0.05], [1.05, 1.05],
                                color=colors["diffs"], zorder=0)

    handles = [patches.Patch(facecolor=used_colors[index][1], edgecolor=used_colors[index][0],
                             linewidth=1, label=labels[index]) for index in [0, 1, 2]]
    pyplot.legend(handles=handles, loc="upper right", fontsize=8)
    pyplot.xlim(-0.5, 11.5)
    pyplot.xticks(range(12), ["" for _ in range(12)], fontsize=8)
    pyplot.ylabel("detection rate", fontsize=8)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=8)

    pyplot.subplot(2, 1, 2)
    for filter_index in range(len(filter_indices)):
        for index, values in enumerate(results[2][filter_index]):
            if sum(values) == 0:
                pyplot.scatter([filter_index + bias[index]], [0],
                               color=used_colors[index][0], marker="x", s=15)
            elif sum(values) == len(values):
                pyplot.scatter([filter_index + bias[index]], [1],
                               color=used_colors[index][0], marker="^", s=15)
            elif max(values) - min(values) < 0.02:
                value = median(values)
                pyplot.hlines(value, filter_index + bias[index] - 0.1, filter_index + bias[index] + 0.1,
                              linewidths=1, edgecolors=used_colors[index][0], zorder=3)
                pyplot.scatter([filter_index + bias[index]], value, color="white", edgecolor=used_colors[index][0],
                               linewidth=1, s=8, zorder=4)
            else:
                violin = pyplot.violinplot(dataset=values, positions=[filter_index + bias[index]],
                                           widths=0.2, bw_method=0.5, showextrema=False)
                for patch in violin["bodies"]:
                    patch.set_edgecolor(used_colors[index][0])
                    patch.set_facecolor(used_colors[index][1])
                    patch.set_linewidth(1)
                    patch.set_alpha(1)
                pyplot.scatter([filter_index + bias[index]], median(values),
                               color="white", edgecolor=used_colors[index][0],
                               linewidth=1, s=8, zorder=4)

        if filter_index % 2 != 0:
            pyplot.fill_between([filter_index - 0.5, filter_index + 0.5], [-4, -4], [84, 84],
                                color=colors["diffs"], zorder=0)

    handles = [patches.Patch(facecolor=used_colors[index][1], edgecolor=used_colors[index][0],
                             linewidth=1, label=labels[index]) for index in [0, 1, 2]]
    pyplot.legend(handles=handles, loc="upper right", fontsize=8)
    pyplot.xlabel("constraint set\n", fontsize=8)
    pyplot.xlim(-0.5, 11.5)
    pyplot.xticks(range(12), filter_indices, fontsize=8)
    pyplot.ylabel("number of repair strategy", fontsize=8)
    pyplot.ylim(-3, 63)
    pyplot.yticks([0, 12, 24, 36, 48, 60], [0, 12, 24, 36, 48, 60], fontsize=8)
    figure.align_labels()
    figure.text(0.02, 0.99, "A", va="center", ha="center")
    figure.text(0.02, 0.53, "B", va="center", ha="center")
    pyplot.savefig("./results/figures/[4-3] repairability normal results.svg",
                   format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


def draw_single_situations(vertex_number):
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

    data = load(file="./results/data/step_4_repairability_single_error.npy")
    results = [[[[], [], []] for _ in range(len(filter_indices))] for _ in range(4)]
    for index in range(0, len(data) - vertex_number * 8 + 1, vertex_number * 8):
        sub_data = data[index: index + vertex_number * 8]
        filter_index = sub_data[0, 0]
        for mark, time in zip([0, 1, 2], [3, 4, 1]):  # mark 0 = substitution, 1 = insertion, 2 = deletion
            locations = where(sub_data[:, 2] == mark)[0]
            detect_rate_1 = sum(sub_data[locations][:, 3]) / (vertex_number * time)
            strategy_numbers = sub_data[locations][:, 4]
            strategy_numbers = strategy_numbers[strategy_numbers > 0]
            detect_rate_2 = sum(sub_data[locations][:, 5]) / (vertex_number * time)
            repair_rate = sum(sub_data[locations][:, 6]) / (vertex_number * time)
            results[0][filter_index - 1][mark].append(detect_rate_1)
            results[1][filter_index - 1][mark].append(detect_rate_2)
            if len(strategy_numbers) > 0:
                results[2][filter_index - 1][mark].append(mean(strategy_numbers))
            results[3][filter_index - 1][mark].append(repair_rate)

    bias = [-0.3, 0, 0.3]
    used_colors = [[colors["algo1"], colors["algo2"]],
                   [colors["trad1"], colors["trad2"]],
                   [colors["foun1"], colors["foun2"]]]
    labels = ["substitution", "insertion", "deletion"]

    pyplot.figure(figsize=(10, 5))
    for filter_index in range(len(filter_indices)):
        for index, values in enumerate(results[3][filter_index]):
            if sum(values) == 0:
                pyplot.scatter([filter_index + bias[index]], [0],
                               color=used_colors[index][0], marker="x", s=15)
            elif sum(values) == len(values):
                pyplot.scatter([filter_index + bias[index]], [1],
                               color=used_colors[index][0], marker="^", s=15)
            elif max(values) - min(values) < 0.02:
                value = median(values)
                pyplot.hlines(value, filter_index + bias[index] - 0.1, filter_index + bias[index] + 0.1,
                              linewidths=1, edgecolors=used_colors[index][0], zorder=3)
                pyplot.scatter([filter_index + bias[index]], value, color="white", edgecolor=used_colors[index][0],
                               linewidth=1, s=8, zorder=4)
            else:
                violin = pyplot.violinplot(dataset=values, positions=[filter_index + bias[index]],
                                           widths=0.2, bw_method=0.5, showextrema=False)
                for patch in violin["bodies"]:
                    patch.set_edgecolor(used_colors[index][0])
                    patch.set_facecolor(used_colors[index][1])
                    patch.set_linewidth(1)
                    patch.set_alpha(1)
                pyplot.scatter([filter_index + bias[index]], median(values),
                               color="white", edgecolor=used_colors[index][0],
                               linewidth=1, s=8, zorder=4)

        if filter_index % 2 != 0:
            pyplot.fill_between([filter_index - 0.5, filter_index + 0.5], [-0.05, -0.05], [1.05, 1.05],
                                color=colors["diffs"], zorder=0)
    handles = [patches.Patch(facecolor=used_colors[index][1], edgecolor=used_colors[index][0],
                             linewidth=1, label=labels[index]) for index in [0, 1, 2]]
    pyplot.legend(handles=handles, loc="upper right", fontsize=8)
    pyplot.xlabel("constraint set\n", fontsize=8)
    pyplot.xlim(-0.5, 11.5)
    pyplot.xticks(range(12), filter_indices, fontsize=8)
    pyplot.ylabel("correction rate", fontsize=8)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=8)

    pyplot.savefig("./results/figures/[4-2] repairability vt-based correction.svg",
                   format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    # create_folders()
    #
    # show_single_examples()
    # show_multiple_examples()
    #
    # single_evaluation(task_seed=2021, repeats=100, vertex_number=100)
    # multiple_evaluation(task_seed=2021, repeats=10000)

    draw_normal_situation(vertex_number=100)
    # draw_single_situations(vertex_number=100)
