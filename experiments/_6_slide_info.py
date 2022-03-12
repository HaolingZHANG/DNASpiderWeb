from collections import Counter
from itertools import product
# noinspection PyPackageRequirements
from matplotlib import pyplot, patches
from numpy import load, array, zeros, zeros_like, linspace, percentile, corrcoef, polyfit, where
from numpy import median, mean, std, sqrt, min, max, sum, argmax, log, log2, log10, clip, abs, inf
from pickle import load as pload
from scipy.optimize import curve_fit
from warnings import filterwarnings

from dsw import dna_to_number, approximate_capacity, adjacency_matrix_to_accessor

from experiments import colors

filterwarnings("ignore")


def stable():
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    display_data = [[[] for _ in range(12)] for _ in range(4)]
    for sample in load(file="./results/data/step_3_stability_evaluation.npy"):
        if sample[2] > 0:
            display_data[int(sample[0])][int(sample[1]) - 1].append((8 * 1024 * 64) / sample[2])
        else:
            display_data[int(sample[0])][int(sample[1]) - 1].append(-1)

    pyplot.figure(figsize=(10, 5))
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
                             linewidth=1, label=labels[index]) for index in list(range(4))[::-1]]
    pyplot.legend(handles=legends, loc="upper left", fontsize=10)
    pyplot.xlim(-0.5, 5.5)
    pyplot.ylim(-0.08, 2.08)
    pyplot.xticks(range(6), filter_indices[:6], fontsize=12)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
                  ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0"], fontsize=12)
    pyplot.xlabel("constraint set", fontsize=12)
    pyplot.ylabel("code rate", fontsize=12)
    pyplot.savefig("./results/figures/1.1.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
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
                             linewidth=1, label=labels[index]) for index in list(range(4))[::-1]]
    pyplot.legend(handles=legends, loc="lower right", fontsize=10)
    pyplot.xlim(5.5, 11.5)
    pyplot.ylim(-0.08, 2.08)
    pyplot.xticks(range(6, 12), filter_indices[6:], fontsize=12)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
                  ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0"], fontsize=12)
    pyplot.xlabel("constraint set", fontsize=12)
    pyplot.ylabel("code rate", fontsize=12)
    pyplot.savefig("./results/figures/1.2.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    with open("./results/data/step_1_compatibility_capacities.pkl", "rb") as file:
        capacities, data = pload(file)
    with open("./results/data/step_1_compatibility_code_rates.pkl", "rb") as file:
        transcode_results = pload(file)

    pyplot.figure(figsize=(10, 5))
    for index, filter_index in enumerate(filter_indices):
        capacity, code_rates = capacities[index], transcode_results[index]
        pyplot.text(x=index, y=capacity + 0.02, s="%.4f" % capacity, ha="center", va="bottom", fontsize=12)
        pyplot.hlines(capacity, index - 0.4, index + 0.4, color=colors["trad1"], linewidth=2)
        if max(code_rates) - min(code_rates) < 0.01 or max(code_rates) > capacity:
            code_rate = median(code_rates)
            if code_rate > capacity:
                code_rate = capacity
            pyplot.hlines(code_rate, index - 0.25, index + 0.25, linewidths=2, edgecolors=colors["algo1"], zorder=3)
            pyplot.scatter([index], code_rate, color="white", edgecolor=colors["algo1"], linewidth=2, s=15, zorder=4)
        else:
            violin = pyplot.violinplot(dataset=code_rates, positions=[index], bw_method=0.5, showextrema=False)
            for patch in violin["bodies"]:
                patch.set_edgecolor(colors["algo1"])
                patch.set_facecolor(colors["algo2"])
                patch.set_linewidth(2)
                patch.set_alpha(1)
            pyplot.scatter([index], median(code_rates), color="white", edgecolor=colors["algo1"],
                           linewidth=2, s=15, zorder=4)
        if index % 2 != 0:
            pyplot.fill_between([index - 0.5, index + 0.5], [0.9, 0.9], [2.1, 2.1],
                                color=colors["diffs"], zorder=0)

    pyplot.xlabel("constraint set", fontsize=12)
    pyplot.xticks(range(12), filter_indices, fontsize=12)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylabel("code rate", fontsize=12)
    pyplot.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0],
                  ["1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0"], fontsize=12)
    pyplot.ylim(0.9, 2.1)
    pyplot.savefig("./results/figures/1.3.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()


def repair():
    with open("./results/data/step_4_repairability_multiple_errors.pkl", "rb") as file:
        task_1, task_2, task_3 = pload(file)

    pyplot.figure(figsize=(10, 5))
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))
    valid_numbers, observed_rates = [], []
    for filter_index in range(1, 13):
        valid_numbers.append(sum(load("./results/data/a" + str(filter_index).zfill(2) + "[v].npy")))
        rates = [sum(task_1[filter_index, error_time][:, -1]) / 2000.0 for error_time in range(1, 5)]
        observed_rates.append(rates[1])
        pyplot.plot(range(4), rates, color=gradient_colors[filter_index - 1],
                    linewidth=2, marker="o", label=str(filter_index).zfill(2))
    pyplot.legend(loc="upper right", ncol=3, fontsize=12)
    pyplot.xlabel("introduced error rate", fontsize=12)
    pyplot.xlim(-0.06, 3.06)
    pyplot.xticks([0, 1, 2, 3], ["1%", "2%", "3%", "4%"], fontsize=12)
    pyplot.ylabel("correction rate", fontsize=12)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks([0, 0.25, 0.5, 0.75, 1], ["0%", "25%", "50%", "75%", "100%"], fontsize=12)
    pyplot.savefig("./results/figures/2.1.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))
    x, y = [], []
    for filter_index in range(1, 13):
        x.append(sum(load("./results/data/a" + str(filter_index).zfill(2) + "[v].npy")))
        y.append(sum(task_1[filter_index, 2][:, -1]) / 2000.0)
        pyplot.scatter([log(float(x[-1])) / log(4)], [y[-1]],
                       color=gradient_colors[filter_index - 1], edgecolor="black", label=str(filter_index).zfill(2),
                       zorder=2)
    x, y = array(x), array(y)
    a, b = polyfit(log(x) / log(4), y, deg=1)
    shown_x = linspace(5, 10, 51)
    shown_y = a * shown_x + b
    pyplot.plot(shown_x, shown_y, color="black", linewidth=2, linestyle="--", zorder=1)
    pyplot.text(shown_x[0], shown_y[0] + 0.05, "Pearson: %.2f" % (corrcoef(x, y)[0, 1]),
                fontsize=12)
    pyplot.legend(loc="upper right", ncol=3, fontsize=12)
    pyplot.xlabel("valid number", fontsize=12)
    pyplot.xticks([5, 6, 7, 8, 9, 10],
                  ["$4^5$", "$4^6$", "$4^7$", "$4^8$", "$4^9$", "$4^{10}$"], fontsize=12)
    pyplot.xlim(4.9, 10.1)
    pyplot.ylabel("correction rate", fontsize=12)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=12)
    pyplot.savefig("./results/figures/2.2.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    rate_matrix = zeros(shape=(4, 4))
    statistics = [[] for _ in range(16)]
    for position_1, dna_length in enumerate(linspace(start=100, stop=400, num=4, dtype=int)):
        for position_2, error_time in enumerate((linspace(1, 4, 4) * dna_length / 100).astype(int)):
            rate_matrix[position_1, position_2] = sum(task_2[dna_length, error_time][:, -1]) / 2000.0
            statistics[error_time - 1].append(rate_matrix[position_1, position_2])
    statistics = [mean(statistics[i]) if len(statistics[i]) > 0 else -1 for i in range(16)]
    x, y = [], []
    for i, value in enumerate(statistics):
        if value >= 0:
            x.append(i + 1)
            y.append(value)

    pyplot.pcolormesh(range(4 + 1), range(4 + 1), rate_matrix.T, vmax=1, vmin=0, cmap="RdYlGn")
    for i in range(4):
        for j in range(4):
            if rate_matrix[i, j] >= 0:
                pyplot.text(x=i + 0.5, y=j + 0.5, s="%.1f" % (rate_matrix[i, j] * 100) + "%",
                            va="center", ha="center", fontsize=12)
            else:
                pyplot.text(x=i + 0.5, y=j + 0.5, s="x",
                            va="center", ha="center", fontsize=12)

    pyplot.xlabel("length of DNA string", fontsize=12)
    pyplot.xlim(0, 4)
    pyplot.xticks([0.5, 1.5, 2.5, 3.5], ["100nt", "200nt", "300nt", "400nt"], fontsize=12)
    pyplot.ylabel("introduced error rate", fontsize=12)
    pyplot.ylim(0, 4)
    pyplot.yticks([0.5, 1.5, 2.5, 3.5], ["1%", "2%", "3%", "4%"], fontsize=12)
    pyplot.savefig("./results/figures/2.3.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    statistics = [[] for _ in range(16)]
    for position_1, dna_length in enumerate(linspace(start=100, stop=400, num=4, dtype=int)):
        for position_2, error_time in enumerate((linspace(1, 4, 4) * dna_length / 100).astype(int)):
            statistics[error_time - 1].append(sum(task_2[dna_length, error_time][:, -1]) / 2000.0)
    statistics = [mean(statistics[i]) if len(statistics[i]) > 0 else -1 for i in range(16)]
    x, y = [], []
    for i, value in enumerate(statistics):
        if value >= 0:
            x.append(i + 1)
            y.append(value)

    pyplot.scatter(x, y, color="silver", edgecolor="black", zorder=2)
    pyplot.plot(x, y, color="black", linewidth=1, zorder=1)

    pyplot.xlabel("introduced error number", fontsize=12)
    pyplot.xticks(range(1, 17), range(1, 17), fontsize=12)
    pyplot.xlim(0.65, 16.35)
    pyplot.ylabel("correction rate", fontsize=12)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=12)
    pyplot.savefig("./results/figures/2.4.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    info = []
    for dna_length in linspace(start=100, stop=400, num=4, dtype=int):
        for error_time in (linspace(1, 4, 4) * dna_length / 100).astype(int):
            if (dna_length == 300 and error_time == 12) \
                    or (dna_length == 400 and error_time == 12) \
                    or (dna_length == 400 and error_time == 16):
                continue
            for data in task_2[dna_length, error_time]:
                info.append([data[1], data[-1]])
    info = array(info).T
    info[0] = info[0] / max(info[0])
    info[0] = (info[0] * 10 + 0.5).astype(int)
    shown_data = zeros(shape=(2, 11))
    for sample in info.T:
        shown_data[0, int(sample[0])] += sample[1]
        shown_data[1, int(sample[0])] += 1

    shown_data[1][shown_data[1] == 0] = 1
    shown_data = shown_data[0] / shown_data[1]
    pyplot.bar(range(11), shown_data, color="silver", edgecolor="black")
    for index, value in enumerate(shown_data):
        pyplot.text(x=index, y=value + 0.02, s=str(int(value * 100 + 0.5)) + "%", ha="center", va="bottom", fontsize=12)

    pyplot.xlabel("normalized coefficient of variation", fontsize=12)
    pyplot.xlim(-0.65, 10.65)
    pyplot.xticks(range(11),
                  ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"], fontsize=12)
    pyplot.ylabel("correction rate", fontsize=12)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=12)
    pyplot.savefig("./results/figures/2.5.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    violin_1, violin_2, points = None, None, None
    for filter_index in range(1, 13):
        formers, latters = task_1[filter_index, 1][:, 3], task_1[filter_index, 1][:, 6]
        formers, latters = log10(formers[formers > 0]), log10(latters[latters > 0])
        violin_1 = pyplot.violinplot(formers, positions=[filter_index], widths=0.8, bw_method=0.5,
                                     showmeans=False, showextrema=False, showmedians=False)
        for body in violin_1["bodies"]:
            center = mean(body.get_paths()[0].vertices[:, 0])
            body.get_paths()[0].vertices[:, 0] = clip(body.get_paths()[0].vertices[:, 0], -inf, center)
            body.set_color(colors["algo1"])
            body.set_edgecolor("black")
            body.set_alpha(1)
        pyplot.scatter([filter_index - 0.15], [mean(formers)], color="white", edgecolor="black", s=10)
        violin_2 = pyplot.violinplot(latters, positions=[filter_index], widths=0.8, bw_method=0.5,
                                     showmeans=False, showextrema=False, showmedians=False)
        for body in violin_2["bodies"]:
            center = mean(body.get_paths()[0].vertices[:, 0])
            body.get_paths()[0].vertices[:, 0] = clip(body.get_paths()[0].vertices[:, 0], center, inf)
            body.set_color(colors["yyco1"])
            body.set_edgecolor("black")
            body.set_alpha(1)
        points = pyplot.scatter([filter_index + 0.15], [mean(latters)], color="white", edgecolor="black", s=10)
    pyplot.legend([violin_1["bodies"][0], violin_2["bodies"][0], points],
                  ["search-only", "combined", "mean value"], fontsize=12)
    pyplot.xlabel("constraint set", fontsize=12)
    pyplot.xticks(range(1, 13), filter_indices, fontsize=12)
    pyplot.xlim(0.5, 12.5)
    pyplot.ylabel("number of correction candidates", fontsize=12)
    pyplot.yticks([0, 1, 2, 3], [1, 10, 100, 1000], fontsize=12)
    pyplot.ylim(-0.2, 3.2)
    pyplot.savefig("./results/figures/2.6.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    rates = []
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    for filter_index in range(1, 13):
        rates.append(sum(task_1[filter_index, 1][:, 2]) / 2000.0)
        if filter_index < 9:
            pyplot.text(filter_index, log10(rates[-1] * 10000) + 0.03, "%.1f" % (rates[-1] * 100) + "%",
                        va="bottom", ha="center", fontsize=12)
        else:
            pyplot.text(filter_index, log10(rates[-1] * 10000) + 0.03, "%.1f" % (rates[-1] * 100) + "%",
                        va="bottom", ha="center", color="red", fontsize=12)
    pyplot.bar(linspace(1, 12, 12), log10(array(rates) * 10000), width=0.6, color="silver", edgecolor="black")
    pyplot.xlabel("constraint set", fontsize=12)
    pyplot.xticks(range(1, 13), filter_indices, fontsize=12)
    pyplot.xlim(0.45, 12.55)
    pyplot.ylabel("detection rate", fontsize=12)
    pyplot.ylim(-0.25, 4.25)
    pyplot.yticks([0, 1, 2, 3, 4], ["0.01%", "0.1%", "1%", "10%", "100%"], fontsize=12)
    pyplot.savefig("./results/figures/2.7.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    lengths = linspace(100, 400, 26)
    styles = ["-", "--", ":"]
    for number in range(1, 4):
        molecules = 1 / task_3[0, number - 1]
        pyplot.plot(lengths, molecules,
                    color="black", linewidth=2, linestyle=styles[number - 1], label=str(number) + " error(s)")

    pyplot.legend(loc="upper right", fontsize=12)
    pyplot.xlabel("DNA string length", fontsize=12)
    pyplot.xlim(100, 400)
    pyplot.xticks([100, 200, 300, 400], [100, 200, 300, 400], fontsize=12)
    pyplot.ylabel("minimum correction read(s)", fontsize=12)
    pyplot.ylim(0, 4)
    pyplot.yticks([0, 1, 2, 3, 4], ["0", "1", "2", "3", "4"], fontsize=12)
    pyplot.savefig("./results/figures/2.8.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    lengths = linspace(100, 400, 26)
    styles = ["-", "--", ":"]
    for number in range(1, 4):
        pyplot.plot(lengths, task_3[1, number - 1],
                    color="black", linewidth=2, linestyle=styles[number - 1], label=str(number) + " error(s)")

    pyplot.legend(loc="upper right", fontsize=12)
    pyplot.xlabel("DNA string length", fontsize=12)
    pyplot.xlim(100, 400)
    pyplot.xticks([100, 200, 300, 400], [100, 200, 300, 400], fontsize=12)
    pyplot.ylabel("average runtime (second)", fontsize=12)
    pyplot.ylim(0, 1)
    pyplot.yticks([0, 0.25, 0.5, 0.75, 1], ["0.00", "0.25", "0.50", "0.75", "1.00"], fontsize=12)
    pyplot.savefig("./results/figures/2.9.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()


def protect():
    def estimated_equation(x, a):
        return a ** x

    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))
    reconstructions = load(file="./results/data/step_5_encrypability_reconstruction.npy")
    with open("./results/data/step_1_compatibility_capacities.pkl", "rb") as file:
        capacities, _ = pload(file)

    pyplot.figure(figsize=(10, 5.5))
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
        numbers = array(list(range(len(medians))))
        parameter = curve_fit(estimated_equation, xdata=medians[1:], ydata=numbers[1:], p0=(2,))[0][0]
        used_capacities = parameter ** percentages * 100 * capacities[index] / 8
        display_data[index] = used_capacities
    for index, locations in display_data.items():
        pyplot.plot(percentages, locations, color=gradient_colors[index], linewidth=2, zorder=2,
                    label=str(index + 1).zfill(2))
    pyplot.legend(loc="upper left", fontsize=12)
    pyplot.xlabel("graph reconstruction percentage", fontsize=12)
    pyplot.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                  ["0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"], fontsize=12)
    pyplot.xlim(0, 1)
    pyplot.ylabel("transmitted data capacity", fontsize=12)
    pyplot.yticks([0.0 * 1e8, 0.2 * 1e8, 0.4 * 1e8, 0.6 * 1e8, 0.8 * 1e8, 1.0 * 1e8,
                   1.2 * 1e8, 1.4 * 1e8, 1.6 * 1e8, 1.8 * 1e8, 2.0 * 1e8],
                  ["0.0E+8 bytes", "0.2E+8 bytes", "0.4E+8 bytes", "0.6E+8 bytes", "0.8E+8 bytes", "1.0E+8 bytes",
                   "1.2E+8 bytes", "1.4E+8 bytes", "1.6E+8 bytes", "1.8E+8 bytes", "2.0E+8 bytes"],
                  fontsize=12)
    pyplot.ylim(0, 2 * 1e8)
    pyplot.savefig("./results/figures/3.1.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
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
            pyplot.plot([index - 0.25, index + 0.25], [min(collected_data[index]), max(collected_data[index])],
                        color=gradient_colors[index], linewidth=2)
            pyplot.plot([index + 0.25, index - 0.25], [min(collected_data[index]), max(collected_data[index])],
                        color=gradient_colors[index], linewidth=2)
        else:
            pyplot.hlines(median(collected_data[index]), index - 0.25, index + 0.25,
                          color=gradient_colors[index], linewidth=2)

    pyplot.xlabel("constraint set", fontsize=12)
    pyplot.xticks(range(12), filter_indices, fontsize=12)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylabel("transmitted data capacity", fontsize=12)
    pyplot.yticks([0, 1, 2, 3], ["B", "KB", "MB", "GB"], fontsize=12)
    pyplot.ylim(0, 3)
    pyplot.savefig("./results/figures/3.2.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    bit_length = 256
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))
    follow_ups = zeros(shape=(12, 3))  # [2, 3, 4]
    for index, filter_index in enumerate(filter_indices):
        accessor = load(file="./results/data/a" + filter_index + "[g].npy")
        out_degrees, counts = array(list(Counter(where(accessor >= 0)[0]).values())), zeros(shape=(3,), dtype=float)
        for out_degree in [2, 3, 4]:
            counts[out_degree - 2] = len(where(out_degrees == out_degree)[0])
        counts /= sum(counts)
        follow_ups[index] = counts

    pyplot.hlines(bit_length, 0, 320, color="silver", linewidth=0.75, linestyle="--", zorder=1)
    pyplot.text(4, 256 + 4, "reference", va="bottom", ha="left", fontsize=12)
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
                           label="[" + filter_indices[index] + "]   %.2f" % reference_length)
        else:
            pyplot.scatter(reference_length, bit_length,
                           color=gradient_colors[index], edgecolor="black", s=20, zorder=4,
                           label="[" + filter_indices[index] + "] %.2f" % reference_length)

    pyplot.legend(loc="upper right", ncol=3, fontsize=12)
    pyplot.xlabel("DNA string length", fontsize=12)
    pyplot.xlim(0, 320)
    pyplot.xticks([0, 64, 128, 192, 256, 320], ["0nt", "64nt", "128nt", "192nt", "256nt", "320nt"], fontsize=12)
    pyplot.ylabel("equivalent key strength", fontsize=12)
    pyplot.ylim(0, 384)
    pyplot.yticks([0, 128, 256, 384], [0, 128, 256, 384], fontsize=12)

    pyplot.savefig("./results/figures/3.3.png",
                   format="png", bbox_inches="tight", dpi=600)
    pyplot.close()


def appendix():
    def create_matrix(dna_strings):
        matrix = zeros(shape=(16, 16), dtype=int)
        vertices = [dna_to_number(dna_string=vertex, is_string=False) for vertex in dna_strings]
        for former, latter in product(range(len(dna_strings)), repeat=2):
            if dna_strings[former][1] == dna_strings[latter][0]:
                matrix[vertices[former], vertices[latter]] = 1
        return matrix

    cases = [create_matrix(dna_strings=["AC",
                                        "CG",
                                        "GT",
                                        "TA"]),  # a cycle.
             create_matrix(dna_strings=["AC", "AG",
                                        "CA", "CG",
                                        "GA", "GT",
                                        "TC", "TG"]),  # GC content is 50%.
             create_matrix(dna_strings=["AC", "AG", "AT",
                                        "CA", "CG", "CT",
                                        "GA", "GC", "GT",
                                        "TA", "TC", "TG"]),  # no homopolymer.
             create_matrix(dna_strings=["AA", "AC", "AG", "AT",
                                        "CA", "CC", "CG", "CT",
                                        "GA", "GC", "GG", "GT",
                                        "TA", "TC", "TG", "TT"])]  # complete graph.
    figure = pyplot.figure(figsize=(10, 4.5), tight_layout=True)
    for index, case in enumerate(cases):
        pyplot.subplot(2, len(cases), index + 1)
        for x in range(16):
            for y in range(16):
                if case[x, y] == 1:
                    pyplot.fill_between(x=[x, x + 1], y1=[y, y], y2=[y + 1, y + 1], color=colors["trad1"])
        for value in range(16):
            pyplot.vlines(value, 0, 16, color="black", linewidth=0.5)
            pyplot.hlines(value, 0, 16, color="black", linewidth=0.5)
        pyplot.xticks(array(list(range(16))) + 0.5,
                      ["0", "", "", "3", "", "", "6", "", "", "9", "", "", "12", "", "", "15"], fontsize=12)
        pyplot.yticks(array(list(range(16))) + 0.5,
                      ["0", "", "", "3", "", "", "6", "", "", "9", "", "", "12", "", "", "15"], fontsize=12)
        pyplot.ylabel("latter index", fontsize=12)
        pyplot.xlabel("former index", fontsize=12)
        pyplot.xlim(0, 16)
        pyplot.ylim(0, 16)
        axes = pyplot.subplot(2, len(cases), index + len(cases) + 1)
        _, process = approximate_capacity(accessor=adjacency_matrix_to_accessor(case), repeats=1, need_process=True)
        if len(process) == 1:
            pyplot.scatter([0], process, color=colors["trad1"])
        pyplot.plot(range(len(process)), process, color=colors["trad1"], marker="o", linewidth=2)
        pyplot.text(len(process) - 0.8, process[-1], "%.3f" % process[-1], va="center", ha="left", fontsize=12)
        pyplot.xlim(-0.1, 2.1)
        pyplot.ylim(-0.1, 2.1)
        pyplot.yticks([0, 1, 2], ["0", "1", "2"], fontsize=12)
        pyplot.ylabel("capacity", fontsize=12)
        pyplot.xticks([0, 1, 2], [1, 2, 3], fontsize=12)
        pyplot.xlabel("iteration", fontsize=12)
        # noinspection PyUnresolvedReferences
        axes.spines["right"].set_visible(False)
        # noinspection PyUnresolvedReferences
        axes.spines["top"].set_visible(False)
    figure.align_labels()
    pyplot.savefig("./results/figures/4.1.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    with open("./results/data/step_2_reliability_detail.pkl", "rb") as file:
        data, errors = pload(file)
    errors = log10(array(errors))
    violin = pyplot.violinplot(dataset=errors, positions=[0], showextrema=False, vert=False, widths=1)
    paint_data = violin["bodies"][0].get_paths()[0].vertices
    paint_data = paint_data[len(paint_data) // 2 + 1: -1]
    pyplot.close()
    x, y = paint_data[:, 0], paint_data[:, 1]
    y = y / sum(y) * 100.0
    [v_25, v_50, v_75] = percentile(a=errors, q=[25, 50, 75])
    lower, upper = v_25 - 1.5 * (v_75 - v_25), 1.5 * (v_75 - v_25) + v_75
    pyplot.figure(figsize=(10, 5))
    pyplot.plot(x, y, color=colors["trad1"], linewidth=1.5, zorder=4)
    pyplot.fill_between(x, y, zeros_like(y), color=colors["trad2"], alpha=0.75, zorder=3)
    pyplot.vlines(x=v_50, ymin=0, ymax=6.5, color="black", linewidth=1, zorder=2)
    pyplot.vlines(x=v_50, ymin=6.9, ymax=8.5, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=8.5, xmin=v_50 - 0.3, xmax=v_50 + 0.3, color="black", linewidth=1, zorder=2)
    pyplot.text(x=v_50, y=8.6, s="median", va="bottom", ha="center", fontsize=12, zorder=2)
    pyplot.vlines(x=v_25, ymin=0, ymax=7, color="black", linewidth=1, zorder=2)
    pyplot.vlines(x=v_75, ymin=0, ymax=7, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=6.7, xmin=v_25, xmax=v_75, color="black", linewidth=1, zorder=2)
    pyplot.text(x=v_75 + 0.1, y=6.7, s="interquartile range", va="center", ha="left", fontsize=12, zorder=2)
    pyplot.vlines(x=lower, ymin=0, ymax=4, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=3.7, xmin=lower, xmax=min(errors) + 2, color="black", linewidth=1, zorder=2)
    pyplot.text(x=min(errors) + 1.9, y=3.7, s="outlier range", va="center", ha="right", fontsize=12, zorder=2)
    pyplot.vlines(x=upper, ymin=0, ymax=4, color="black", linewidth=1, zorder=2)
    pyplot.hlines(y=3.7, xmin=upper, xmax=max(errors) - 2, color="black", linewidth=1, zorder=2)
    pyplot.text(x=max(errors) - 1.9, y=3.7, s="outlier range", va="center", ha="left", fontsize=12, zorder=2)
    pyplot.xlim(-14.1, -3.9)
    pyplot.ylim(-0.3, 10.3)
    pyplot.ylabel("frequency", fontsize=12)
    pyplot.xticks([-14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4],
                  ["$10^{-14}$", "$10^{-13}$", "$10^{-12}$", "$10^{-11}$", "$10^{-10}$",
                   "$10^{-09}$", "$10^{-08}$", "$10^{-07}$", "$10^{-06}$", "$10^{-05}$", "$10^{-04}$"],
                  fontsize=12)
    pyplot.xlabel("relative error", fontsize=12)
    pyplot.yticks([0, 2, 4, 6, 8, 10], ["0%", "2%", "4%", "6%", "8%", "10%"], fontsize=12)
    pyplot.savefig("./results/figures/4.2.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    with open("./results/data/step_2_reliability_detail.pkl", "rb") as file:
        data, errors = pload(file)
    errors = log10(array(errors))
    pyplot.scatter(errors, data, color=colors["trad1"], edgecolor="black", zorder=3)
    temp_errors, temp_data = [], []
    for value_1, value_2 in sorted(zip(errors, data), key=lambda v: v[0]):
        temp_errors.append(value_1)
        temp_data.append(value_2)
    errors, data = array(temp_errors), array(temp_data)
    a, b = polyfit(errors, data, deg=1)
    estimated = a * errors + b
    bounds = abs(std(errors) * sqrt(1.0 / len(errors)
                                    + (errors - mean(errors)) ** 2 / sum((errors - mean(errors)) ** 2)))
    pyplot.fill_between(errors, estimated - bounds, estimated + bounds, color=colors["trad2"], alpha=0.75, zorder=1)
    x = linspace(min(errors), max(errors), 100)
    y = a * x + b
    pyplot.plot(x, y, color="black", linestyle="--", linewidth=1, zorder=2)

    pyplot.xlim(-14.1, -3.9)
    pyplot.ylim(-0.1, 2.1)
    pyplot.xticks([-14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4],
                  ["$10^{-14}$", "$10^{-13}$", "$10^{-12}$", "$10^{-11}$", "$10^{-10}$",
                   "$10^{-09}$", "$10^{-08}$", "$10^{-07}$", "$10^{-06}$", "$10^{-05}$", "$10^{-04}$"],
                  fontsize=12)
    pyplot.yticks([0, 0.5, 1.0, 1.5, 2.0], ["0.0", "0.5", "1.0", "1.5", "2.0"], fontsize=12)
    pyplot.xlabel("relative error", fontsize=12)
    pyplot.ylabel("capacity from SPIDER-WEB", fontsize=12)
    pyplot.savefig("./results/figures/4.3.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 5))
    with open("./results/data/step_2_reliability_extend.pkl", "rb") as file:
        whole_data = pload(file)
    pyplot.fill_between([1.5, 2.5], [-15, -15], [-3, -3], color=colors["diffs"], zorder=0)
    pyplot.fill_between([3.5, 4.5], [-15, -15], [-3, -3], color=colors["diffs"], zorder=0)
    used = []
    for length, data in whole_data.items():
        used.append(log10(data[1]))

    violin = pyplot.violinplot(dataset=used, positions=[0.9, 1.9, 2.9, 3.9, 4.9], showextrema=False, widths=0.35)
    for patch in violin["bodies"]:
        patch.set_edgecolor(colors["trad1"])
        patch.set_facecolor(colors["trad2"])
        patch.set_linewidth(1)
        patch.set_alpha(1)

    for index, data in enumerate(used):
        [v_25, v_50, v_75] = percentile(a=data, q=[25, 50, 75])
        pyplot.scatter(x=[index + 0.9], y=[v_50],
                       color="white", edgecolor=colors["trad1"], linewidth=1, s=25, zorder=4)
        lower, upper = v_25 - 1.5 * (v_75 - v_25), 1.5 * (v_75 - v_25) + v_75
        pyplot.vlines(x=index + 1.2, ymin=lower, ymax=upper,
                      color=colors["trad1"], linewidth=1, zorder=2)
        pyplot.hlines(y=lower, xmin=index + 1.2 - 0.04, xmax=index + 1.2 + 0.04,
                      color=colors["trad1"], linewidth=1, zorder=2)
        pyplot.hlines(y=upper, xmin=index + 1.2 - 0.04, xmax=index + 1.2 + 0.04,
                      color=colors["trad1"], linewidth=1, zorder=2)
        pyplot.fill_between(x=[index + 1.2 - 0.04, index + 1.2 + 0.04], y1=[v_25, v_25], y2=[v_75, v_75],
                            color=colors["trad1"], zorder=3)
        for value in data:
            if value < lower or value > upper:
                pyplot.scatter(x=[index + 1.2], y=[value],
                               color="white", edgecolor=colors["trad1"], linewidth=1, s=20, zorder=4)
    pyplot.xlim(0.5, 5.5)
    pyplot.ylim(-15, -3)
    pyplot.xticks([1, 2, 3, 4, 5],
                  [r"$4^2 \times 4^2$", r"$4^3 \times 4^3$", r"$4^4 \times 4^4$",
                   r"$4^5 \times 4^5$", r"$4^6 \times 4^6$"], fontsize=12)
    pyplot.yticks([-15, -13, -11, -9, -7, -5, -3],
                  ["$10^{-15}$", "$10^{-13}$", "$10^{-11}$", "$10^{-09}$", "$10^{-07}$",
                   "$10^{-05}$", "$10^{-03}$"], fontsize=12)
    pyplot.xlabel("size of adjacency matrix", fontsize=12)
    pyplot.ylabel("relative error", fontsize=12)
    pyplot.savefig("./results/figures/4.4.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    stable()
    repair()
    protect()
    appendix()
