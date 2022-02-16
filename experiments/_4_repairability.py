# noinspection PyPackageRequirements
from matplotlib import pyplot, patches
from numpy import load, save, zeros, ones, array, linspace, concatenate, vstack, random, where
from numpy import median, sum, max, min, maximum, abs, std, mean, corrcoef
from os import path
from pickle import load as pload
from pickle import dump as psave
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.stats import poisson as p

from warnings import simplefilter

from experiments import colors, create_folders
from experiments.code_repair import show_single_examples, show_multiple_examples
from experiments.code_repair import evaluate_single_error, evaluate_repair_multiple_errors


simplefilter("ignore", OptimizeWarning)


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
    info = {}
    for error_time in range(1, 9):
        for dna_length in dna_lengths:
            save_path = "./results/temp/multiple" + str(dna_length).zfill(4) + "." + str(error_time).zfill(2) + ".pkl"
            if not path.exists(save_path):
                records = evaluate_repair_multiple_errors(random_seed=task_seed, accessor=accessor, vertices=vertices,
                                                          observed_length=10, vt_length=10,
                                                          repeats=repeats, dna_length=dna_length,
                                                          error_times=error_time, check_iterations=error_time + 1)
                with open(save_path, "wb") as file:
                    psave(records, file)
            else:
                with open(save_path, "rb") as file:
                    records = pload(file)

            info[(error_time, dna_length)] = records

    with open("./results/data/step_4_repairability_multiple_errors.pkl", "wb") as file:
        psave(info, file)


def draw_total_evaluation():
    def poisson(k, lamb):
        return p.pmf(k=k, mu=lamb)

    def square(x, a):
        return a * x ** 2

    with open("./results/data/step_4_repairability_multiple_errors.pkl", "rb") as file:
        records = pload(file)

    fitted_data = []
    figure = pyplot.figure(figsize=(10, 3.2), tight_layout=True)

    pyplot.subplot(1, 2, 1)

    saved_data, total_data = zeros(shape=(8, 4)), zeros(shape=(8, 4))
    for (error_time, dna_length), data_group in records.items():
        for data in data_group:
            if data[0] < 0.045:
                ll = int(dna_length / 100) - 1
                er = 0 if data[0] < 0.015 else (1 if data[0] < 0.025 else (2 if data[0] < 0.035 else 3))
                if data[4] and data[5] <= 1:
                    saved_data[ll, er] += 1
                total_data[ll, er] += 1
    total_data[total_data == 0] = 1.0
    shown_data = vstack((ones(shape=(1, 8)), (saved_data / total_data).T))

    # data volume of 300nt ~ 800nt in 2% ~ 4% error rate are decreased gradually, need fitting.
    for i in range(2, 8):
        parameter = curve_fit(poisson, [1, 2, 3, 4, 5], shown_data.T[i])[0][0]
        shown_data.T[i] = maximum(shown_data.T[i], array([poisson(value, parameter) for value in [1, 2, 3, 4, 5]]))

    fitted_data.append(shown_data)
    pyplot.pcolormesh(range(8 + 1), range(5 + 1), shown_data, vmax=1, vmin=0, cmap="RdYlGn")
    for i in range(8):
        for j in range(5):
            pyplot.text(x=i + 0.5, y=j + 0.5, s="%.1f" % (shown_data.T[i, j] * 100) + "%",
                        va="center", ha="center", fontsize=8)

    pyplot.xlabel("DNA string length", fontsize=8)
    pyplot.xlim(0, 8)
    pyplot.xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5],
                  ["100nt", "200nt", "300nt", "400nt", "500nt", "600nt", "700nt", "800nt"], fontsize=8)
    pyplot.ylabel("error rate in DNA string", fontsize=8)
    pyplot.ylim(0, 5)
    pyplot.yticks([0.5, 1.5, 2.5, 3.5, 4.5], ["0%", "1%", "2%", "3%", "4%"], fontsize=8)

    pyplot.subplot(1, 2, 2)
    saved_data, total_data = zeros(shape=(8, 4)), zeros(shape=(8, 4))
    for (error_time, dna_length), data_group in records.items():
        for data in data_group:
            if data[3] and data[0] < 0.045:
                ll = int(dna_length / 100) - 1
                er = 0 if data[0] < 0.015 else (1 if data[0] < 0.025 else (2 if data[0] < 0.035 else 3))
                if data[-1] > 1:
                    saved_data[ll, er] += 1
                total_data[ll, er] += 1
    total_data[total_data == 0] = 1.0
    shown_data = vstack((zeros(shape=(1, 8)), (saved_data / total_data).T))

    # data volume of 300nt ~ 800nt in 2% ~ 4% error rate are decreased gradually, need fitting.
    shown_data.T[2:, 4], shown_data.T[2:, 3], shown_data.T[4:, 2] = 0, 0, 0  # ignore distorted data.
    parameter = curve_fit(square, list(range(1, 8)), shown_data.T[1:, 1])[0][0]
    shown_data.T[:, 1] = array([0] + [square(value, parameter) for value in range(1, 8)])
    parameter = curve_fit(square, list(range(1, 4)), shown_data.T[1: 4, 2])[0][0]
    shown_data.T[:, 2] = array([0] + [square(value, parameter) for value in range(1, 8)])
    parameter = curve_fit(square, list(range(1, 2)), shown_data.T[1:2, 3])[0][0]
    shown_data.T[:, 3] = array([0] + [square(value, parameter) for value in range(1, 8)])
    parameter = curve_fit(square, list(range(1, 2)), shown_data.T[1:2, 4])[0][0]
    shown_data.T[:, 4] = array([0] + [square(value, parameter) for value in range(1, 8)])

    fitted_data.append(shown_data)
    pyplot.pcolormesh(range(8 + 1), range(5 + 1), shown_data, vmax=1, vmin=0, cmap="RdYlBu_r")
    for i in range(8):
        for j in range(5):
            pyplot.text(x=i + 0.5, y=j + 0.5, s="%.1f" % (shown_data.T[i, j] * 100) + "%",
                        va="center", ha="center", fontsize=8)

    pyplot.xlabel("DNA string length", fontsize=8)
    pyplot.xlim(0, 8)
    pyplot.xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5],
                  ["100nt", "200nt", "300nt", "400nt", "500nt", "600nt", "700nt", "800nt"], fontsize=8)
    pyplot.ylabel("error rate in DNA string", fontsize=8)
    pyplot.ylim(0, 5)
    pyplot.yticks([0.5, 1.5, 2.5, 3.5, 4.5], ["0%", "1%", "2%", "3%", "4%"], fontsize=8)

    figure.text(0.019, 0.98, "A", va="center", ha="center")
    figure.text(0.512, 0.98, "B", va="center", ha="center")

    pyplot.savefig("./results/figures/[4-1] repairability evaluation.svg",
                   format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()

    return fitted_data


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

    figure = pyplot.figure(figsize=(10, 8), tight_layout=True)
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
    pyplot.xlabel("constraint set", fontsize=8)
    pyplot.xlim(-0.5, 11.5)
    pyplot.xticks(range(12), filter_indices, fontsize=8)
    pyplot.ylabel("number of repair strategy", fontsize=8)
    pyplot.ylim(-3, 63)
    pyplot.yticks([0, 12, 24, 36, 48, 60], [0, 12, 24, 36, 48, 60], fontsize=8)
    figure.align_labels()
    figure.text(0.02, 0.99, "A", va="center", ha="center")
    figure.text(0.02, 0.52, "B", va="center", ha="center")
    pyplot.savefig("./results/figures/[4-2] repairability normal results.svg",
                   format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


def draw_pearson(fitted_data):
    with open("./results/data/step_4_repairability_multiple_errors.pkl", "rb") as file:
        records = pload(file)

    digital_mapping = {"-": 0, "A": 1, "C": 2, "G": 3, "T": 4}
    data_1, data_2 = {}, {}
    for (error_time, dna_length), data_group in records.items():
        for data in data_group:
            right = array(list(map(digital_mapping.get, data[1])))
            wrong = array(list(map(digital_mapping.get, data[2])))
            value = (where(abs(right - wrong).astype(bool).astype(int) == 1)[0] / dna_length)
            if len(value) > 0:
                value = std(value) / mean(value)  # variation coefficient
                if value in data_1:
                    data_1[value].append(data[4] and data[5] <= 1)
                else:
                    data_1[value] = [data[4] and data[5] <= 1]

                if value in data_2:
                    data_2[value].append(data[-1] > 1)
                else:
                    data_2[value] = [data[-1] > 1]

    info_1, info_2 = zeros(shape=(51, 2)), zeros(shape=(51, 2))
    for key, value in data_1.items():
        info_1[int((key / 2.0) * 50.0 + 0.5), 0] += sum(value)
        info_1[int((key / 2.0) * 50.0 + 0.5), 1] += len(value)
    info_1 = info_1.T
    info_1[1, info_1[1] == 0] = 1.0

    for key, value in data_2.items():
        info_2[int((key / 2.0) * 50.0 + 0.5), 0] += sum(value)
        info_2[int((key / 2.0) * 50.0 + 0.5), 1] += len(value)
    info_2 = info_2.T
    info_2[1, info_2[1] == 0] = 1.0

    numbers = []
    for time in range(5):
        for length in range(100, 801, 100):
            numbers.append(time / 100 * length)
    numbers = array(numbers)

    figure = pyplot.figure(figsize=(10, 8), tight_layout=True)

    pyplot.subplot(2, 2, 1)
    pyplot.scatter(numbers, fitted_data[0].reshape(-1), color=colors["yyco1"], edgecolor="black", s=25)
    pyplot.xlabel("introduced error number", fontsize=8)
    pyplot.xlim(-1, 33)
    pyplot.xticks(range(0, 33, 4), range(0, 33, 4), fontsize=8)
    pyplot.ylabel("average correction rate", fontsize=8)
    pyplot.ylim(-0.04, 1.04)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=8)
    print(corrcoef(numbers, fitted_data[0].reshape(-1))[0, 1])

    pyplot.subplot(2, 2, 2)
    pyplot.scatter(linspace(0, 2, 51), info_1[0] / info_1[1], color=colors["yyco1"], edgecolor="black", s=25)
    pyplot.xlabel("average variation coefficient of errors", fontsize=8)
    pyplot.xlim(-0.08, 2.08)
    pyplot.xticks([0, 0.4, 0.8, 1.2, 1.6, 2.0], ["0.0", "0.4", "0.8", "1.2", "1.6", "2.0"], fontsize=8)
    pyplot.ylabel("average correction rate", fontsize=8)
    pyplot.ylim(-0.04, 1.04)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=8)
    print(corrcoef(linspace(0, 2, 51), info_1[0] / info_1[1])[0, 1])

    pyplot.subplot(2, 2, 3)
    pyplot.scatter(numbers, fitted_data[1].reshape(-1), color=colors["trad1"], edgecolor="black", s=25)
    pyplot.xlabel("introduced error number", fontsize=8)
    pyplot.xlim(-1, 33)
    pyplot.xticks(range(0, 33, 4), range(0, 33, 4), fontsize=8)
    pyplot.ylabel("average false positive rate", fontsize=8)
    pyplot.ylim(-0.04, 1.04)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=8)
    print(corrcoef(numbers, fitted_data[1].reshape(-1))[0, 1])

    pyplot.subplot(2, 2, 4)
    pyplot.scatter(linspace(0, 2, 51), info_2[0] / info_2[1], color=colors["trad1"], edgecolor="black", s=25)
    pyplot.xlabel("average variation coefficient of errors", fontsize=8)
    pyplot.xlim(-0.08, 2.08)
    pyplot.xticks([0, 0.4, 0.8, 1.2, 1.6, 2.0], ["0.0", "0.4", "0.8", "1.2", "1.6", "2.0"], fontsize=8)
    pyplot.ylabel("average false positive rate", fontsize=8)
    pyplot.ylim(-0.04, 1.04)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=8)
    print(corrcoef(linspace(0, 2, 51), info_2[0] / info_2[1])[0, 1])

    figure.align_labels()
    figure.text(0.019, 0.99, "A", va="center", ha="center")
    figure.text(0.512, 0.99, "B", va="center", ha="center")
    figure.text(0.019, 0.50, "C", va="center", ha="center")
    figure.text(0.512, 0.50, "D", va="center", ha="center")

    pyplot.savefig("./results/figures/[4-3] repairability correlation.svg",
                   format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()

    show_single_examples()
    show_multiple_examples()

    single_evaluation(task_seed=2021, repeats=100, vertex_number=100)
    multiple_evaluation(task_seed=2021, repeats=10000)

    fit_data = draw_total_evaluation()
    draw_normal_situation(vertex_number=100)
    draw_pearson(fitted_data=fit_data)
