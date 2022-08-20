from collections import Counter
# noinspection PyPackageRequirements
from matplotlib import pyplot, rcParams
from numpy import load, array, zeros_like, zeros, ones, arange, linspace, percentile, random, where, inf
from numpy import median, min, max, mean, sum, abs, sqrt, std, log10, clip, polyfit
# noinspection PyPackageRequirements
from openpyxl import Workbook
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
from warnings import filterwarnings

from dsw import obtain_vertices

filterwarnings("ignore")


def supp01():
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    rates = []

    primal_graphs = load("./data/primal_graphs.pkl", allow_pickle=True)
    for filter_index in filter_indices:
        accessor = primal_graphs[filter_index]
        out_degrees = sum(where(accessor >= 0, 1, 0), axis=1)
        rates.append(len(out_degrees[out_degrees == 1]) / len(out_degrees[out_degrees >= 1]))

    figure = pyplot.figure(figsize=(10, 8), tight_layout=True)
    rcParams["font.family"] = "Times New Roman"
    pyplot.subplot(2, 1, 1)
    visited_indices = []
    for index in range(12):
        if rates[index] > 0:
            pyplot.bar([len(visited_indices)], [rates[index]], color="silver", edgecolor="black", width=0.6)
            pyplot.text(len(visited_indices), rates[index] + 0.001, "%.1f" % (rates[index] * 100) + "%", fontsize=12,
                        va="bottom", ha="center")
            visited_indices.append(filter_indices[index])

    pyplot.xlabel("selected constraint set index", fontsize=12)
    pyplot.xticks(range(len(visited_indices)), visited_indices, fontsize=12)
    pyplot.xlim(-0.35, 3.35)
    pyplot.ylabel("probability of out-degree 1", fontsize=12)
    pyplot.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10], ["0%", "2%", "4%", "6%", "8%", "10%"], fontsize=12)
    pyplot.ylim(-0.005, 0.105)

    proposed_results, out_1_results = load("./data/performance_evaluation_2.pkl", allow_pickle=True)

    pyplot.subplot(2, 1, 2)
    visited_indices = []
    for index, (filter_index, out_values) in enumerate(out_1_results.items()):
        visited_indices.append(filter_index)
        value_1, value_2 = std(out_values), std(proposed_results[int(filter_index) - 1])
        if index > 0:
            pyplot.bar([index - 0.15], [value_1], color="#81B8DF", edgecolor="black", width=0.3)
            pyplot.bar([index + 0.15], [value_2], color="#FE817D", edgecolor="black", width=0.3)
        else:
            pyplot.bar([index - 0.15], [value_1],
                       color="#81B8DF", edgecolor="black", width=0.3, label="retain vertex of out-degree 1")
            pyplot.bar([index + 0.15], [value_2],
                       color="#FE817D", edgecolor="black", width=0.3, label="remove vertex of out-degree 1")
        pyplot.text(index - 0.15, value_1 + 0.001, "%.3f" % value_1, fontsize=12, va="bottom", ha="center")
        pyplot.text(index + 0.15, value_2 + 0.001, "%.3f" % value_2, fontsize=12, va="bottom", ha="center")

    pyplot.legend(loc="upper right", ncol=2, fontsize=12)
    pyplot.xlabel("selected constraint set index", fontsize=12)
    pyplot.xticks(range(len(visited_indices)), visited_indices, fontsize=12)
    pyplot.xlim(-0.35, 3.35)
    pyplot.ylabel("standard deviation of code rate", fontsize=12)
    pyplot.yticks([0, 0.02, 0.04, 0.06, 0.08, 0.10], ["0.00", "0.02", "0.04", "0.06", "0.08", "0.10"], fontsize=12)
    pyplot.ylim(-0.005, 0.105)

    figure.align_labels()
    figure.text(0.024, 1.00, "a", va="center", ha="center", fontsize=14)
    figure.text(0.024, 0.51, "b", va="center", ha="center", fontsize=14)

    pyplot.savefig("./show/supp01.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def supp02():
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    capacities, _ = load(file="./data/capacity_reference.pkl", allow_pickle=True)
    transcode_results = load(file="./data/performance_evaluation_2.pkl", allow_pickle=True)[0]

    pyplot.figure(figsize=(10, 8))
    rcParams["font.family"] = "Times New Roman"
    for index, filter_index in enumerate(filter_indices):
        capacity, code_rates = capacities[index], transcode_results[index]
        pyplot.text(x=index, y=capacity + 0.01, s="%.4f" % capacity, ha="center", va="bottom", fontsize=12)
        pyplot.hlines(capacity, index - 0.4, index + 0.4, color="#81B8DF", linewidth=1)
        if max(code_rates) - min(code_rates) < 0.01 or max(code_rates) > capacity:
            code_rate = median(code_rates)
            if code_rate > capacity:
                code_rate = capacity
            pyplot.hlines(code_rate, index - 0.25, index + 0.25, linewidths=1, edgecolors="#FE817D", zorder=3)
            pyplot.scatter([index], code_rate, color="white", edgecolor="#FE817D", linewidth=1, s=15, zorder=4)
        else:
            violin = pyplot.violinplot(dataset=code_rates, positions=[index], bw_method=0.5, showextrema=False)
            for patch in violin["bodies"]:
                patch.set_edgecolor("#FE817D")
                patch.set_facecolor("#FCBBAE")
                patch.set_linewidth(1)
                patch.set_alpha(1)
            pyplot.scatter([index], median(code_rates), color="white", edgecolor="#FE817D",
                           linewidth=1, s=15, zorder=4)
        if index % 2 != 0:
            pyplot.fill_between([index - 0.5, index + 0.5], [0.92, 0.92], [2.08, 2.08], color="#F1F1F1", zorder=0)

    pyplot.xlabel("constraint set index", fontsize=12)
    pyplot.xticks(range(12), filter_indices, fontsize=12)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylabel("code rate", fontsize=12)
    pyplot.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0],
                  ["1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0"], fontsize=12)
    pyplot.ylim(0.95, 2.05)
    pyplot.savefig("./show/supp02.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def supp03():
    random_seed, param_number = 2021, 100

    random.seed(random_seed)
    fountain_params = random.uniform(low=array([0.1, 0.01]), high=array([1, 0.1]), size=(param_number, 2)).tolist()

    random.seed(random_seed)
    total_indices = array(list(range(6144)))
    random.shuffle(total_indices)
    yinyang_params = total_indices[:param_number]
    yinyang_params = yinyang_params.tolist()

    hedges_params = []
    correct_penalties = [-0.035, -0.082, -0.127, -0.229, -0.265, -0.324]
    for pattern_index, correct_penalty in enumerate(correct_penalties):
        hedges_params.append([pattern_index + 1, correct_penalty])
    hedges_params = array(hedges_params).tolist()

    spiderweb_params = []
    coding_graphs = load(file="./data/coding_graphs.pkl", allow_pickle=True)
    for filter_index in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
        vertices = obtain_vertices(accessor=coding_graphs[filter_index])
        random.seed(random_seed)
        random.shuffle(vertices)
        chosen_indices = vertices[:param_number]
        spiderweb_params.append(chosen_indices)
    spiderweb_params = array(spiderweb_params).T.tolist()

    book = Workbook()
    sheet = book.create_sheet(title="supp03", index=0)
    sheet.append(["param index", "DNA Fountain", "", "Yin-Yang Code", "HEDGES", "",
                  "SPIDER-WEB", "", "", "", "", "", "", "", "", "", "", ""])
    sheet.append(["", "c", "delta", "rule id", "pattern", "correct penalty",
                  "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12"])
    for param_index in range(param_number):
        record = [str(param_index + 1)]
        for value in fountain_params[param_index]:
            record.append("%.2f" % value)
        record.append(str(yinyang_params[param_index]))
        if param_index < 6:
            record.append(str(int(hedges_params[param_index][0])))
            record.append("%.3f" % hedges_params[param_index][1])
        else:
            record.append("N/A")
            record.append("N/A")
        for vertex_index in spiderweb_params[param_index]:
            record.append(str(vertex_index))

        sheet.append(record)

    book.save("./show/supp03.xlsx")


def supp04():
    display_data = [[[] for _ in range(4)] for _ in range(12)]
    record = load(file="./data/performance_evaluation_1.pkl", allow_pickle=True)
    for algorithm_index, filter_index, _, _, code_rate in record:
        display_data[int(filter_index) - 1][int(algorithm_index) - 1].append(code_rate)

    book = Workbook()
    sheet = book.create_sheet(title="supp04", index=0)
    sheet.append(["constraint set index", "DNA Fountain", "Yin-Yang Code", "HEDGES", "SPIDER-WEB"])
    for index, data in enumerate(display_data):
        record = [str(index + 1).zfill(2)]
        for method_index, sub_data in enumerate(data):
            if sum(sub_data) == -len(sub_data):
                record.append("N/A")
            elif max(sub_data) - min(sub_data) < 0.001:
                record.append("%.2f" % std(sub_data))
            else:
                sub_data = array(sub_data)
                sub_data = sub_data[sub_data > 0]
                record.append("%.4f" % std(sub_data))
        sheet.append(record)

    book.save("./show/supp04.xlsx")


def supp05():
    counts = ones(shape=(12, 6)) * 100
    record = load("./data/performance_evaluation_1.pkl", allow_pickle=True)
    for algorithm_index, filter_index, param_index, _, code_rate in record:
        if algorithm_index == 3 and code_rate == -1:  # HEDGES
            counts[int(filter_index) - 1, int(param_index) - 1] -= 1
    counts /= 100.0

    book = Workbook()
    sheet = book.create_sheet(title="supp05", index=0)
    sheet.append(["constraint set index", "pattern 1", "pattern 2", "pattern 3", "pattern 4", "pattern 5", "pattern 6"])
    for index, data in enumerate(counts):
        record = [str(index + 1).zfill(2)]
        for value in data:
            record.append(str(int(value * 100)) + "%")
        sheet.append(record)

    book.save("./show/supp05.xlsx")


def supp06():
    print("main06.pdf is the \"case of path-based error correcting\", created by PowerPoint.")


def supp07():
    data = load(file="./data/correction_evaluation_2.pkl", allow_pickle=True)

    pyplot.figure(figsize=(10, 8), tight_layout=True)
    rcParams["font.family"] = "Times New Roman"
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    values, rates = [[], []], [[], []]
    for filter_index in range(1, 13):
        formers, latters = data[str(filter_index).zfill(2)]
        formers, latters = log10(formers), log10(latters)

        violin_1 = pyplot.violinplot(formers, positions=[filter_index], widths=0.8, bw_method=0.5,
                                     showmeans=False, showextrema=False, showmedians=False)
        for body in violin_1["bodies"]:
            center = mean(body.get_paths()[0].vertices[:, 0])
            body.get_paths()[0].vertices[:, 0] = clip(body.get_paths()[0].vertices[:, 0], -inf, center)
            body.set_color("#FE817D")
            body.set_edgecolor("black")
            body.set_linewidth(1.5)
            body.set_alpha(1)
        pyplot.scatter([filter_index], [median(formers)], color="#FE817D", edgecolor="black", linewidth=1.5, s=40)
        values[0].append(mean(formers))

        violin_2 = pyplot.violinplot(latters, positions=[filter_index], widths=0.8, bw_method=0.5,
                                     showmeans=False, showextrema=False, showmedians=False)
        for body in violin_2["bodies"]:
            center = mean(body.get_paths()[0].vertices[:, 0])
            body.get_paths()[0].vertices[:, 0] = clip(body.get_paths()[0].vertices[:, 0], center, inf)
            body.set_color("#00C382")
            body.set_edgecolor("black")
            body.set_linewidth(1.5)
            body.set_alpha(1)

        pyplot.scatter([filter_index], [median(latters)], color="#00C382", edgecolor="black", linewidth=1.5, s=40)
        values[1].append(mean(latters))

        if filter_index == 1:
            pyplot.legend([violin_1["bodies"][0], violin_2["bodies"][0]], ["search-only", "combined"], fontsize=12)

    pyplot.xlabel("constraint set index", fontsize=12)
    pyplot.xticks(range(1, 13), filter_indices, fontsize=12)
    pyplot.xlim(0.5, 12.5)
    pyplot.ylabel("candidate number", fontsize=12)
    pyplot.yticks([0, 1, 2, 3], ["1", "10", "100", "1000"], fontsize=12)
    pyplot.ylim(-0.2, 3.2)
    pyplot.savefig("./show/supp07.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def supp08():
    data = load(file="./data/correction_evaluation_3.pkl", allow_pickle=True)

    figure = pyplot.figure(figsize=(10, 4.5), tight_layout=True)
    rcParams["font.family"] = "Times New Roman"

    pyplot.subplot(1, 2, 1)
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))
    for index in range(12):
        pyplot.plot(range(4), data["constraint influence"][index], color=gradient_colors[index],
                    linewidth=2, marker="o", label=str(index + 1).zfill(2))
    pyplot.legend(loc="upper right", ncol=2, fontsize=12)
    pyplot.xlabel("error rate", fontsize=12)
    pyplot.xlim(-0.1, 3.1)
    pyplot.xticks([0, 1, 2, 3], ["1%", "2%", "3%", "4%"], fontsize=12)
    pyplot.ylabel("correction rate", fontsize=12)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks([0, 0.25, 0.5, 0.75, 1], ["0%", "25%", "50%", "75%", "100%"], fontsize=12)

    pyplot.subplot(1, 2, 2)
    pyplot.pcolormesh(range(4 + 1), range(4 + 1), data["length influence"].T, vmax=1, vmin=0, cmap="RdYlGn")
    for i in range(4):
        for j in range(4):
            if data["length influence"][i, j] >= 0:
                pyplot.text(x=i + 0.5, y=j + 0.5, s="%.1f" % (data["length influence"][i, j] * 100) + "%",
                            va="center", ha="center", fontsize=12)
            else:
                pyplot.text(x=i + 0.5, y=j + 0.5, s="x",
                            va="center", ha="center", fontsize=12)

    pyplot.xlabel("DNA string length", fontsize=12)
    pyplot.xlim(0, 4)
    pyplot.xticks([0.5, 1.5, 2.5, 3.5], ["100nt", "200nt", "300nt", "400nt"], fontsize=12)
    pyplot.ylabel("error rate", fontsize=12)
    pyplot.ylim(0, 4)
    pyplot.yticks([0.5, 1.5, 2.5, 3.5], ["1%", "2%", "3%", "4%"], fontsize=12)

    figure.align_labels()
    figure.text(0.02, 0.98, "a", va="center", ha="center", fontsize=14)
    figure.text(0.52, 0.98, "b", va="center", ha="center", fontsize=14)

    pyplot.savefig("./show/supp08.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def supp09():
    data = load(file="./data/correction_evaluation_1.pkl", allow_pickle=True)

    figure = pyplot.figure(figsize=(10, 11), tight_layout=True)
    rcParams["font.family"] = "Times New Roman"

    for location, establish_reads in enumerate([10, 20, 50, 100]):
        pyplot.subplot(4, 1, location + 1)
        pyplot.title("error rate = " + str(location + 1) + "% | established reads = " + str(establish_reads),
                     fontsize=12)
        obtained_set, values = data["frequency-based recovery"][location + 1][1], []
        for index, (_, count) in enumerate(obtained_set.items()):
            values.append(count)

        counts, location, firsts = Counter(values), 1, [True, True]
        maximum_size = max(list(counts.keys()))
        for value in arange(1, maximum_size + 1, 1)[::-1]:
            if value in counts:
                if location <= 10000:
                    if firsts[0]:
                        pyplot.hlines(value, log10(location), log10(location + counts[value]), color="#FE817D",
                                      linewidth=4, label="right", zorder=3)
                        firsts[0] = False
                    else:
                        pyplot.hlines(value, log10(location), log10(location + counts[value]), color="#FE817D",
                                      linewidth=4, zorder=3)
                else:
                    if firsts[1]:
                        pyplot.hlines(value, log10(location), log10(location + counts[value]), color="#81B8DF",
                                      linewidth=4, label="wrong", zorder=3)
                        firsts[1] = False
                    else:
                        pyplot.hlines(value, log10(location), log10(location + counts[value]), color="#81B8DF",
                                      linewidth=4, zorder=3)
                location += counts[value]

        pyplot.fill_between([4, 6], -maximum_size / 10, maximum_size + maximum_size / 10,
                            color="#EEEEEE", label="ignored range", zorder=0)
        pyplot.legend(loc="upper right", ncol=2, fontsize=12)
        pyplot.xlabel("priority of repaired DNA strings", fontsize=12)
        pyplot.xlim(-0.2, 6.2)
        pyplot.xticks([0, 1, 2, 3, 4, 5, 6], ["1E+0", "1E+1", "1E+2", "1E+3", "1E+4", "1E+5", "1E+6"], fontsize=12)
        pyplot.ylabel("frequency", fontsize=12)
        pyplot.ylim(-maximum_size / 10, maximum_size + maximum_size / 10)
        pyplot.yticks(linspace(0, maximum_size, 6), linspace(0, maximum_size, 6, dtype=int), fontsize=12)

    figure.align_labels()
    figure.text(0.025, 0.98, "a", va="center", ha="center", fontsize=14)
    figure.text(0.025, 0.73, "b", va="center", ha="center", fontsize=14)
    figure.text(0.025, 0.48, "c", va="center", ha="center", fontsize=14)
    figure.text(0.025, 0.24, "d", va="center", ha="center", fontsize=14)

    pyplot.savefig("./show/supp09.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def supp10():
    records = load(file="./data/correction_evaluation_4.pkl", allow_pickle=True)
    matrix_1, matrix_2, values = zeros(shape=(12, 6)), zeros(shape=(12, 6)), []
    for filter_index in range(1, 13):
        for pattern_index in range(1, 7):
            values.append(records[(filter_index, pattern_index)][:, 1].tolist())
            matrix_1[filter_index - 1, pattern_index - 1] = min(records[(filter_index, pattern_index)][:, 1])
            matrix_2[filter_index - 1, pattern_index - 1] = max(records[(filter_index, pattern_index)][:, 1])
    values = array(values)

    print("%.2f" % mean(values) + " seconds per sequence.")

    book = Workbook()
    sheet = book.create_sheet(title="supp10", index=0)
    sheet.append(["constraint set index", "pattern 1", "pattern 2", "pattern 3", "pattern 4", "pattern 5", "pattern 6"])
    for filter_index in range(12):
        record = [str(filter_index + 1).zfill(2)]
        for pattern_index in range(6):
            record.append("%.2f" % matrix_1[filter_index - 1, pattern_index - 1] + "~" +
                          "%.2f" % matrix_2[filter_index - 1, pattern_index - 1])
        sheet.append(record)

    book.save("./show/supp10.xlsx")


def supp11():
    def estimated_equation(x, a):
        return a ** x

    reconstructions = load(file="./data/reconstruction_evaluation.npy")
    capacities, _ = load(file="./data/capacity_reference.pkl", allow_pickle=True)
    accessors = load(file="./data/coding_graphs.pkl", allow_pickle=True)

    display_data, percentages = {}, linspace(0, 1, 101)
    for index, data in enumerate(reconstructions):
        medians = [0]
        for iteration_data in data.T[1:]:
            used_data = iteration_data[iteration_data > 0]
            if len(used_data) == 0:
                break
            medians.append(median(used_data))

        number = len(obtain_vertices(accessors[str(index + 1).zfill(2)]))
        medians = array(medians) / number

        numbers = array(list(range(len(medians))))
        parameter = curve_fit(estimated_equation, xdata=medians[1:], ydata=numbers[1:], p0=(2,))[0][0]
        used_capacities = parameter ** percentages * 100 * capacities[index] / 8
        display_data[index] = used_capacities

    pyplot.figure(figsize=(10, 8))
    rcParams["font.family"] = "Times New Roman"
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))

    for index, locations in display_data.items():
        pyplot.plot(percentages, locations, color=gradient_colors[index], linewidth=2, zorder=2,
                    label="constraint set " + str(index + 1).zfill(2))

    pyplot.legend(loc="upper left", fontsize=12)
    pyplot.xlabel("graph reconstruction percentage", fontsize=12)
    pyplot.xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                  ["0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"], fontsize=12)
    pyplot.xlim(0, 1)
    pyplot.ylabel("transmitted file size", fontsize=12)
    pyplot.yticks([0.0 * 1e8, 0.2 * 1e8, 0.4 * 1e8, 0.6 * 1e8, 0.8 * 1e8, 1.0 * 1e8,
                   1.2 * 1e8, 1.4 * 1e8, 1.6 * 1e8, 1.8 * 1e8, 2.0 * 1e8],
                  ["0.0E+8 bytes", "0.2E+8 bytes", "0.4E+8 bytes", "0.6E+8 bytes", "0.8E+8 bytes", "1.0E+8 bytes",
                   "1.2E+8 bytes", "1.4E+8 bytes", "1.6E+8 bytes", "1.8E+8 bytes", "2.0E+8 bytes"],
                  fontsize=12)
    pyplot.ylim(0, 2 * 1e8)

    pyplot.savefig("./show/supp11.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def supp12():
    record = load(file="./data/capacity_evaluation.pkl", allow_pickle=True)

    data, errors = record["detail"]
    errors = log10(array(errors))
    x = linspace(min(errors), max(errors), 101)
    y = gaussian_kde(errors)(x)
    y = y / sum(y) * 100.0

    v_25, v_50, v_75 = percentile(a=errors, q=[25, 50, 75])
    lower, upper = v_25 - 1.5 * (v_75 - v_25), 1.5 * (v_75 - v_25) + v_75

    figure = pyplot.figure(figsize=(10, 11), tight_layout=True)
    rcParams["font.family"] = "Times New Roman"

    pyplot.subplot(3, 1, 1)

    pyplot.plot(x, y, color="#81B8DF", linewidth=1.5, zorder=4)
    pyplot.fill_between(x, y, zeros_like(y), color="#B1CCDF", alpha=0.75, zorder=3)

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
    pyplot.ylim(-0.5, 10.5)
    pyplot.ylabel("frequency", fontsize=12)
    pyplot.xticks([-14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4],
                  ["1E\N{MINUS SIGN}14", "1E\N{MINUS SIGN}13", "1E\N{MINUS SIGN}12", "1E\N{MINUS SIGN}11",
                   "1E\N{MINUS SIGN}10", "1E\N{MINUS SIGN}09", "1E\N{MINUS SIGN}08", "1E\N{MINUS SIGN}07",
                   "1E\N{MINUS SIGN}06", "1E\N{MINUS SIGN}05", "1E\N{MINUS SIGN}04"],
                  fontsize=12)
    pyplot.xlabel("relative error", fontsize=12)
    pyplot.yticks([0, 2, 4, 6, 8, 10], ["0%", "2%", "4%", "6%", "8%", "10%"], fontsize=12)

    pyplot.subplot(3, 1, 2)
    pyplot.scatter(errors, data, color="#81B8DF", edgecolor="black", zorder=3)

    temp_errors, temp_data = [], []
    for value_1, value_2 in sorted(zip(errors, data), key=lambda v: v[0]):
        temp_errors.append(value_1)
        temp_data.append(value_2)
    errors, data = array(temp_errors), array(temp_data)

    a, b = polyfit(errors, data, deg=1)
    estimated = a * errors + b
    bounds = abs(std(errors) * sqrt(1.0 / len(errors)
                                    + (errors - mean(errors)) ** 2 / sum((errors - mean(errors)) ** 2)))
    pyplot.fill_between(errors, estimated - bounds, estimated + bounds, color="#B1CCDF", alpha=0.75, zorder=1)
    x = linspace(min(errors), max(errors), 100)
    y = a * x + b
    pyplot.plot(x, y, color="black", linestyle="--", linewidth=1, zorder=2)

    pyplot.xlabel("relative error", fontsize=12)
    pyplot.ylabel("capacity from SPIDER-WEB", fontsize=12)
    pyplot.xlim(-14.1, -3.9)
    pyplot.ylim(-0.1, 2.1)
    pyplot.xticks([-14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4],
                  ["1E\N{MINUS SIGN}14", "1E\N{MINUS SIGN}13", "1E\N{MINUS SIGN}12", "1E\N{MINUS SIGN}11",
                   "1E\N{MINUS SIGN}10", "1E\N{MINUS SIGN}09", "1E\N{MINUS SIGN}08", "1E\N{MINUS SIGN}07",
                   "1E\N{MINUS SIGN}06", "1E\N{MINUS SIGN}05", "1E\N{MINUS SIGN}04"],
                  fontsize=12)
    pyplot.yticks([0, 0.5, 1.0, 1.5, 2.0], ["0.0", "0.5", "1.0", "1.5", "2.0"], fontsize=12)
    whole_data = record["extend"]

    pyplot.subplot(3, 1, 3)

    pyplot.fill_between([1.5, 2.5], [-15, -15], [-3, -3], color="#F1F1F1", zorder=0)
    pyplot.fill_between([3.5, 4.5], [-15, -15], [-3, -3], color="#F1F1F1", zorder=0)

    used = []
    for length, data in whole_data.items():
        used.append(log10(data[1]))

    violin = pyplot.violinplot(dataset=used, positions=[0.9, 1.9, 2.9, 3.9, 4.9], showextrema=False, widths=0.35)
    for patch in violin["bodies"]:
        patch.set_edgecolor("#81B8DF")
        patch.set_facecolor("#B1CCDF")
        patch.set_linewidth(1)
        patch.set_alpha(1)

    for index, data in enumerate(used):
        [v_25, v_50, v_75] = percentile(a=data, q=[25, 50, 75])
        pyplot.scatter(x=[index + 0.9], y=[v_50],
                       color="white", edgecolor="#81B8DF", linewidth=1, s=25, zorder=4)
        lower, upper = v_25 - 1.5 * (v_75 - v_25), 1.5 * (v_75 - v_25) + v_75
        pyplot.vlines(x=index + 1.2, ymin=lower, ymax=upper,
                      color="#81B8DF", linewidth=1, zorder=2)
        pyplot.hlines(y=lower, xmin=index + 1.2 - 0.04, xmax=index + 1.2 + 0.04,
                      color="#81B8DF", linewidth=1, zorder=2)
        pyplot.hlines(y=upper, xmin=index + 1.2 - 0.04, xmax=index + 1.2 + 0.04,
                      color="#81B8DF", linewidth=1, zorder=2)
        pyplot.fill_between(x=[index + 1.2 - 0.04, index + 1.2 + 0.04], y1=[v_25, v_25], y2=[v_75, v_75],
                            color="#81B8DF", zorder=3)
        for value in data:
            if value < lower or value > upper:
                pyplot.scatter(x=[index + 1.2], y=[value],
                               color="white", edgecolor="#81B8DF", linewidth=1, s=20, zorder=4)
    pyplot.xlabel("size of adjacency matrix", fontsize=12)
    pyplot.ylabel("relative error", fontsize=12)
    pyplot.xlim(0.5, 5.5)
    pyplot.ylim(-15, -3)
    pyplot.xticks([1, 2, 3, 4, 5],
                  ["4^2", "4^3", "4^4", "4^5", "4^6"],
                  fontsize=12)
    pyplot.yticks([-15, -13, -11, -9, -7, -5, -3],
                  ["1E\N{MINUS SIGN}15", "1E\N{MINUS SIGN}13", "1E\N{MINUS SIGN}11", "1E\N{MINUS SIGN}09",
                   "1E\N{MINUS SIGN}07", "1E\N{MINUS SIGN}05", "1E\N{MINUS SIGN}03"],
                  fontsize=12)

    figure.text(0.024, 0.99, "a", va="center", ha="center", fontsize=14)
    figure.text(0.024, 0.66, "b", va="center", ha="center", fontsize=14)
    figure.text(0.024, 0.34, "c", va="center", ha="center", fontsize=14)

    figure.align_labels()
    pyplot.savefig("./show/supp12.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    # supp01()
    # supp02()
    # supp03()
    # supp04()
    # supp05()
    # supp06()
    # supp07()
    # supp08()
    # supp09()
    supp10()
    # supp11()
    # supp12()
