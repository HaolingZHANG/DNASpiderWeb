from collections import Counter
# noinspection PyPackageRequirements
from matplotlib import pyplot, rcParams, patches
from numpy import load, zeros, array, linspace, arange, argmax, where
from numpy import maximum, min, mean, median, max, sum, log, log2, log10
# noinspection PyPackageRequirements
from openpyxl import Workbook
from scipy.optimize import curve_fit
from warnings import filterwarnings

from dsw import obtain_vertices
from experiments import local_bio_filters

filterwarnings("ignore")


def main01():
    print("main01.pdf is the \"illustration of SPIDER-WEB\", created by PowerPoint.")


def main02():
    capacities, _ = load(file="./data/capacity_reference.pkl", allow_pickle=True)

    book = Workbook()
    sheet = book.create_sheet(title="main03", index=0)
    sheet.append(["constraint set index", "long homopolymer", "regionalized GC content", "undesired motifs",
                  "capacity"])

    for constraint_index, bio_filter in local_bio_filters.items():
        save_data = [constraint_index]
        row_info = constraint_index.ljust(20)
        if bio_filter.max_homopolymer_runs is not None:
            row_info += " | " + str(bio_filter.max_homopolymer_runs).ljust(20)
            save_data.append(str(bio_filter.max_homopolymer_runs))
        else:
            row_info += " | " + "N/A".ljust(20)
            save_data.append("N/A")

        if bio_filter.gc_range is not None:
            if bio_filter.gc_range[0] != bio_filter.gc_range[1]:
                gc_range = [str(int(bio_filter.gc_range[0] * 100)) + "%", str(int(bio_filter.gc_range[1] * 100)) + "%"]
                save_data.append(gc_range[0] + " ~ " + gc_range[1])
                row_info += " | " + (gc_range[0] + " ~ " + gc_range[1]).ljust(25)
            else:
                row_info += " | " + (str(int(bio_filter.gc_range[0] * 100)) + "%").ljust(25)
                save_data.append(str(int(bio_filter.gc_range[0] * 100)) + "%")
        else:
            row_info += " | " + "N/A".ljust(25)
            save_data.append("N/A")

        if bio_filter.undesired_motifs is not None:
            row_info += " | " + str(bio_filter.undesired_motifs)[1:-1].replace("\"", "").replace("\'", "").ljust(100)
            save_data.append(str(bio_filter.undesired_motifs)[1:-1].replace("\"", "").replace("\'", ""))
        else:
            row_info += " | " + "N/A".ljust(100)
            save_data.append("N/A")

        row_info += " | " + "%.4f" % capacities[int(constraint_index) - 1]
        save_data.append("%.4f" % capacities[int(constraint_index) - 1])

        sheet.append(save_data)

    book.save("./show/main02.xlsx")


def main03():
    display_data = [[[] for _ in range(4)] for _ in range(12)]
    record = load(file="./data/performance_evaluation_1.pkl", allow_pickle=True)
    for algorithm_index, filter_index, _, _, code_rate in record:
        display_data[int(filter_index - 1)][int(algorithm_index - 1)].append(code_rate)
    book = Workbook()
    sheet = book.create_sheet(title="main04", index=0)
    sheet.append(["constraint set index", "DNA Fountain", "Yin-Yang Code", "HEDGES", "SPIDER-WEB"])
    for index, data in enumerate(display_data):
        record = [str(index + 1).zfill(2)]
        for method_index, sub_data in enumerate(data):
            if sum(sub_data) == -len(sub_data):
                record.append("N/A")
            elif max(sub_data) - min(sub_data) < 0.001:
                record.append("%.2f" % mean(sub_data))
            else:
                sub_data = array(sub_data)
                sub_data = sub_data[sub_data > 0]
                record.append("%.2f ~ %.2f" % (min(sub_data), max(sub_data)))
        sheet.append(record)

    book.save("./show/main03.xlsx")


def main04():
    figure = pyplot.figure(figsize=(10, 6), tight_layout=True)
    rcParams["font.family"] = "Times New Roman"

    data = load(file="./data/correction_evaluation_1.pkl", allow_pickle=True)

    pyplot.subplot(2, 2, 1)
    for error_times in [1, 2, 3, 4]:
        info = maximum(data["correction rate"][error_times][:, -1], 1e-3)
        print("error = " + str(error_times) + "%")
        print("%.2f" % mean(data["correction rate"][error_times][:, -1]) + " on average")
        print("%.2f" % max(data["correction rate"][error_times][:, -1]) + " on maximum")
        print()

        info[info == 0] = 1e-3
        pyplot.boxplot(log10(info), positions=[error_times - 1], widths=0.4, showfliers=False, whis=(0, 100),
                       boxprops=dict(color="#FE817D", facecolor="#FCBBAE", linewidth=1),
                       medianprops=dict(color="black", linewidth=3), patch_artist=True)
        pyplot.text(error_times - 1, max(log10(info)), "%.2f" % max(info), va="bottom", ha="center", fontsize=12)

    pyplot.xlabel("error rate", fontsize=12)
    pyplot.ylabel("log10(seconds)", fontsize=12)
    pyplot.xlim(-0.5, 3.5)
    pyplot.ylim(-4, 1)
    pyplot.xticks([0, 1, 2, 3], ["1%", "2%", "3%", "4%"], fontsize=12)
    pyplot.yticks([-4, -3, -2, -1, 0, 1],
                  ["\N{MINUS SIGN}4", "\N{MINUS SIGN}3", "\N{MINUS SIGN}2", "\N{MINUS SIGN}1", "+0", "+1"], fontsize=12)

    pyplot.subplot(2, 2, 2)
    values = [sum(data["correction rate"][error_times][:, -2]) / 10000.0 for error_times in [1, 2, 3, 4]]
    pyplot.bar(arange(4), values, color="#FCBBAE", edgecolor="#FE817D", linewidth=1)
    for error_times in range(1, 4):
        pyplot.text(error_times, values[error_times], "%.2f" % (values[error_times] * 100) + "%",
                    va="bottom", ha="center", fontsize=12)
    pyplot.xlabel("error rate", fontsize=12)
    pyplot.ylabel("correction rate", fontsize=12)
    pyplot.xlim(-0.6, 3.6)
    pyplot.ylim(0, 1)
    pyplot.xticks([0, 1, 2, 3], ["1%", "2%", "3%", "4%"], fontsize=12)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=12)

    pyplot.subplot(2, 2, 3)
    locations = [1.0, 0.5, 0.2, 0.06]
    terminals = []
    for index, error_times in enumerate([1, 2, 3, 4]):
        info = data["minimum reads"][error_times]
        values = sorted(Counter(info).items(), key=lambda kv: (kv[0], kv[1]))
        values = array(values).T.astype(float)
        x, y = values[0], values[1] / sum(values[1]) / locations[index] * 0.7
        pyplot.plot(x, y + index + 0.1, color="#FE817D", zorder=2)
        pyplot.fill_between(x, index + 0.1, y + index + 0.1, color="#FE817D", linewidth=0, alpha=0.5, zorder=2)
        terminals.append(int(max(x)))

        pyplot.fill_between([1, 90], index + 0.1, index + 0.8, color="#EEEEEE", zorder=1)
        pyplot.vlines(max(x), 0.1, index + 0.1, color="black", linewidth=0.75, linestyle="--", zorder=1)
        pyplot.text(88, index + 0.65, "error rate = " + str(error_times) + "%", va="top", ha="right", fontsize=12,
                    zorder=4, bbox=dict(boxstyle="round", ec="black", fc="white"))

    pyplot.xlabel("pretreatment-free reads for recovery", fontsize=12)
    pyplot.ylabel("proportion", fontsize=12)
    pyplot.xlim(1, 90)
    pyplot.ylim(0.1, 3.8)
    pyplot.xticks(terminals, terminals, fontsize=12)
    pyplot.yticks([0.1, 0.8, 1.1, 1.8, 2.1, 2.8, 3.1, 3.8], ["0%", "100%", "0%", "50%", "0", "20%", "0%", "6%"],
                  fontsize=12)

    pyplot.subplot(2, 2, 4)
    obtained_set, values = data["frequency-based recovery"][2][1], []
    for index, (_, count) in enumerate(obtained_set.items()):
        values.append(count)
    counts, location, firsts = Counter(values), 1, [True, True]
    for value in arange(1, 21, 1)[::-1]:
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
                    pyplot.hlines(value, log10(location + 1), log10(location + counts[value]), color="#81B8DF",
                                  linewidth=4, zorder=3)
            location += counts[value]
    pyplot.legend(loc="upper right", fontsize=12)
    pyplot.xlabel("priority of repaired DNA strings", fontsize=12)
    pyplot.xlim(-0.2, 5.2)
    pyplot.xticks([0, 1, 2, 3, 4, 5], [1, 10, 100, 1000, 10000, 100000], fontsize=12)
    pyplot.ylabel("frequency", fontsize=12)
    pyplot.ylim(-1, 21)
    pyplot.yticks([0, 4, 8, 12, 16, 20], [0, 4, 8, 12, 16, 20], fontsize=12)

    figure.align_labels()
    figure.text(0.023, 0.99, "a", va="center", ha="center", fontsize=14)
    figure.text(0.507, 0.99, "b", va="center", ha="center", fontsize=14)
    figure.text(0.023, 0.49, "c", va="center", ha="center", fontsize=14)
    figure.text(0.507, 0.49, "d", va="center", ha="center", fontsize=14)

    figure.patches.extend([patches.Rectangle((0.08, 0.1775), 0.02, 0.02, fill=True, facecolor="white", zorder=10,
                                             transform=figure.transFigure, figure=figure)])
    figure.patches.extend([patches.Rectangle((0.48, 0.1775), 0.02, 0.02, fill=True, facecolor="white", zorder=10,
                                             transform=figure.transFigure, figure=figure)])
    figure.patches.extend([patches.Rectangle((0.08, 0.2790), 0.02, 0.02, fill=True, facecolor="white", zorder=10,
                                             transform=figure.transFigure, figure=figure)])
    figure.patches.extend([patches.Rectangle((0.48, 0.2790), 0.02, 0.02, fill=True, facecolor="white", zorder=10,
                                             transform=figure.transFigure, figure=figure)])
    figure.patches.extend([patches.Rectangle((0.08, 0.3805), 0.02, 0.02, fill=True, facecolor="white", zorder=10,
                                             transform=figure.transFigure, figure=figure)])
    figure.patches.extend([patches.Rectangle((0.48, 0.3805), 0.02, 0.02, fill=True, facecolor="white", zorder=10,
                                             transform=figure.transFigure, figure=figure)])

    pyplot.savefig("./show/main04.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def main05():
    # bit length or key strength is 256.
    # British Standards Institution (2020). Cryptographic Mechanisms: Recommendations and Key Lengths.
    bit_length = 256

    def estimated_equation(x, a):
        return a ** x

    reconstructions = load(file="./data/reconstruction_evaluation.npy")
    capacities, _ = load("./data/capacity_reference.pkl", allow_pickle=True)
    accessors = load("./data/coding_graphs.pkl", allow_pickle=True)

    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))

    figure = pyplot.figure(figsize=(10, 4), tight_layout=True)
    rcParams["font.family"] = "Times New Roman"

    pyplot.subplot(1, 2, 1)

    collected_data = zeros(shape=(12, 100))
    for index, data in enumerate(reconstructions):
        number = len(obtain_vertices(accessors[str(index + 1).zfill(2)]))
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

    pyplot.xlabel("constraint set index", fontsize=12)
    pyplot.xticks(range(12), filter_indices, fontsize=12)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylabel("transmitted file size", fontsize=12)
    pyplot.yticks([0, 1, 2, 3, 4], ["B", "KB", "MB", "GB", "TB"], fontsize=12)
    pyplot.ylim(0, 4)

    pyplot.subplot(1, 2, 2)
    follow_ups = zeros(shape=(12, 3))  # [2, 3, 4]
    for index, filter_index in enumerate(filter_indices):
        accessor = accessors[filter_index]
        out_degrees, counts = array(list(Counter(where(accessor >= 0)[0]).values())), zeros(shape=(3,), dtype=float)
        for out_degree in [2, 3, 4]:
            counts[out_degree - 2] = len(where(out_degrees == out_degree)[0])
        counts /= sum(counts)
        follow_ups[index] = counts

    pyplot.hlines(bit_length, 0, 320, color="silver", linewidth=0.75, linestyle="--", zorder=1)
    pyplot.text(4, 256 + 4, "AES-256", va="bottom", ha="left", fontsize=12)
    for index, filter_index in enumerate(filter_indices):
        lengths = linspace(0, 320, 321)
        rate = follow_ups[index][0] * log2(2.0) + follow_ups[index][1] * log2(6.0) + follow_ups[index][2] * log2(24.0)
        pyplot.plot(lengths, lengths * rate, color=gradient_colors[index], alpha=0.5, linewidth=2, zorder=2)

        reference_length = bit_length / rate
        pyplot.vlines(reference_length, 0, bit_length,
                      color="silver", linewidth=0.75, linestyle="--", zorder=1)
        if index > 0:
            pyplot.scatter(reference_length, bit_length,
                           color=gradient_colors[index], edgecolor="black", s=20, zorder=4,
                           label="[" + str(filter_index) + "]   %.1f" % reference_length)
        else:
            pyplot.scatter(reference_length, bit_length,
                           color=gradient_colors[index], edgecolor="black", s=20, zorder=4,
                           label="[" + str(filter_index) + "] %.1f" % reference_length)

    pyplot.legend(loc="upper right", ncol=2, fontsize=12)
    pyplot.xlabel("DNA string length", fontsize=12)
    pyplot.xlim(0, 320)
    pyplot.xticks([0, 64, 128, 192, 256, 320], ["0nt", "64nt", "128nt", "192nt", "256nt", "320nt"], fontsize=12)
    pyplot.ylabel("equivalent key strength", fontsize=12)
    pyplot.ylim(0, 512)
    pyplot.yticks([0, 128, 256, 384, 512], [0, 128, 256, 384, 512], fontsize=12)

    figure.align_labels()
    figure.text(0.024, 0.98, "a", va="center", ha="center", fontsize=14)
    figure.text(0.504, 0.98, "b", va="center", ha="center", fontsize=14)

    pyplot.savefig("./show/main05.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    main01()
    main02()
    main03()
    main04()
    main05()
