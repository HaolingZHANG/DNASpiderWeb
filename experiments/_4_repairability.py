# noinspection PyPackageRequirements
from matplotlib import pyplot, rcParams
from numpy import load, save, zeros, array, linspace, polyfit, corrcoef, where
from numpy import sum, max, mean, median, clip, inf, log, log10
from os import path
from pickle import load as pload
from pickle import dump as psave

from experiments import colors, create_folders
from experiments.code_repair import show_single_examples, show_multiple_examples
from experiments.code_repair import evaluate_repair_multiple_errors


def multiple_evaluation(task_seed, repeats):
    if not path.exists("./results/data/step_4_repairability_multiple_errors.pkl"):
        filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
        task_1, task_2, task_3, algorithms = {}, {}, {}, zeros(shape=(2, 26, 4))
        for filter_index in filter_indices:
            algorithms[filter_index] = (load("./results/data/a" + filter_index + "[g].npy"),
                                        where(load("./results/data/a" + filter_index + "[v].npy") == 1)[0])

        # task I: constraint set versus correction rate.
        for filter_index in filter_indices:
            for error_time in [1, 2, 3, 4]:
                save_path = "./results/temp/multiple." + filter_index + ".100." + str(error_time).zfill(2) + ".npy"
                if not path.exists(save_path):
                    records, _ = evaluate_repair_multiple_errors(random_seed=task_seed, error_times=error_time,
                                                                 accessor=algorithms[filter_index][0],
                                                                 vertices=algorithms[filter_index][1],
                                                                 observed_length=10, repeats=repeats,
                                                                 dna_length=100, check_iterations=error_time + 1)
                    save(file=save_path, arr=records)
                else:
                    records = load(file=save_path)

                task_1[(int(filter_index), error_time)] = records

        # task II: DNA length and error rate versus correction rate.
        for dna_length in linspace(start=100, stop=400, num=4, dtype=int):
            for error_time in (linspace(1, 4, 4) * dna_length / 100).astype(int):
                save_path = "./results/temp/multiple.01." \
                            + str(dna_length).zfill(3) + "." + str(error_time).zfill(2) + ".npy"
                if not path.exists(save_path):
                    records, _ = evaluate_repair_multiple_errors(random_seed=task_seed, error_times=error_time,
                                                                 accessor=algorithms["01"][0],
                                                                 vertices=algorithms["01"][1],
                                                                 observed_length=10, repeats=repeats,
                                                                 dna_length=dna_length, check_iterations=error_time + 1)
                    save(file=save_path, arr=records)
                else:
                    records = load(file=save_path)

                task_2[(dna_length, error_time)] = records

        # task III: detailed task II for theoretical comparison about average correction rate and average runtime.
        for dna_length in linspace(100, 400, 26, dtype=int):
            for error_time in [1, 2, 3, 4]:
                save_path = "./results/temp/detailed.01." \
                            + str(dna_length).zfill(3) + "." + str(error_time).zfill(2) + ".npy"
                if not path.exists(save_path):
                    records, times = evaluate_repair_multiple_errors(random_seed=task_seed, error_times=error_time,
                                                                     accessor=algorithms["01"][0],
                                                                     vertices=algorithms["01"][1],
                                                                     observed_length=10, repeats=repeats,
                                                                     dna_length=dna_length,
                                                                     check_iterations=error_time + 1)
                    records = records[:, -1]
                    rate = sum(records) / len(records)
                    save(file=save_path, arr=array([rate, times]))
                else:
                    records = load(file=save_path)

                task_3[0, dna_length, error_time] = records[0]
                task_3[1, dna_length, error_time] = records[1]

        with open("./results/data/step_4_repairability_multiple_errors.pkl", "wb") as file:
            psave((task_1, task_2, task_3), file)


def draw_main():
    with open("./results/data/step_4_repairability_multiple_errors.pkl", "rb") as file:
        task_1, task_2, task_3 = pload(file)

    figure = pyplot.figure(figsize=(10, 7), tight_layout=True)
    rcParams["font.family"] = "Linux Libertine"

    pyplot.subplot(2, 2, 1)
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))
    valid_numbers, observed_rates = [], []
    for filter_index in range(1, 13):
        valid_numbers.append(sum(load("./results/data/a" + str(filter_index).zfill(2) + "[v].npy")))
        rates = [sum(task_1[filter_index, error_time][:, -1]) / 2000.0 for error_time in range(1, 5)]
        observed_rates.append(rates[1])
        pyplot.plot(range(4), rates, color=gradient_colors[filter_index - 1],
                    linewidth=2, marker="o", label=str(filter_index).zfill(2))

    pyplot.legend(loc="upper right", ncol=2, fontsize=12)
    pyplot.xlabel("introduced error rate", fontsize=14)
    pyplot.xlim(-0.075, 3.075)
    pyplot.xticks([0, 1, 2, 3], ["1%", "2%", "3%", "4%"], fontsize=14)
    pyplot.ylabel("correction rate", fontsize=14)
    pyplot.ylim(-0.05, 1.05)
    pyplot.yticks([0, 0.25, 0.5, 0.75, 1], ["0%", "25%", "50%", "75%", "100%"], fontsize=14)

    print("A", corrcoef(valid_numbers, observed_rates)[0, 1])

    pyplot.subplot(2, 2, 2)
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

    pyplot.xlabel("DNA string length", fontsize=14)
    pyplot.xlim(0, 4)
    pyplot.xticks([0.5, 1.5, 2.5, 3.5], ["100nt", "200nt", "300nt", "400nt"], fontsize=14)
    pyplot.ylabel("introduced error rate", fontsize=14)
    pyplot.ylim(0, 4)
    pyplot.yticks([0.5, 1.5, 2.5, 3.5], ["1%", "2%", "3%", "4%"], fontsize=14)

    print("B", corrcoef(x, y)[0, 1])

    pyplot.subplot(2, 2, 3)
    lengths = linspace(100, 400, 26)
    styles = ["-", "--", ":"]
    maximum_molecules = []
    for number in range(1, 4):
        molecules = 1 / task_3[0, number - 1]
        pyplot.plot(lengths, molecules,
                    color="black", linewidth=1.5, linestyle=styles[number - 1], label=str(number) + " error(s)")
        maximum_molecules.append(max(molecules))

    pyplot.legend(loc="upper right", fontsize=12)
    pyplot.xlabel("DNA string length", fontsize=14)
    pyplot.xlim(100, 400)
    pyplot.xticks([100, 200, 300, 400], ["100nt", "200nt", "300nt", "400nt"], fontsize=14)
    pyplot.ylabel("minimum read(s) for correction", fontsize=14)
    pyplot.ylim(0, 4)
    pyplot.yticks([0, 1, 2, 3, 4], ["0", "1", "2", "3", "4"], fontsize=14)

    print("C", ["%.2f" % maximum_molecule for maximum_molecule in maximum_molecules])

    pyplot.subplot(2, 2, 4)
    styles = ["-", "--", ":"]
    maximum_seconds = []
    for number in range(1, 4):
        maximum_seconds.append(max(task_3[1, number - 1]))
        pyplot.plot(lengths, task_3[1, number - 1],
                    color="black", linewidth=1.5, linestyle=styles[number - 1], label=str(number) + " error(s)")

    pyplot.legend(loc="upper right", fontsize=12)
    pyplot.xlabel("DNA string length", fontsize=14)
    pyplot.xlim(100, 400)
    pyplot.xticks([100, 200, 300, 400], ["100nt", "200nt", "300nt", "400nt"], fontsize=14)
    pyplot.ylabel("average runtime (second)", fontsize=14)
    pyplot.ylim(0, 1)
    pyplot.yticks([0, 0.25, 0.5, 0.75, 1], ["0.00", "0.25", "0.50", "0.75", "1.00"], fontsize=14)

    print("D", maximum_seconds)

    figure.align_labels()

    figure.text(0.024, 0.99, "A", va="center", ha="center", fontsize=14)
    figure.text(0.520, 0.99, "B", va="center", ha="center", fontsize=14)
    figure.text(0.024, 0.51, "C", va="center", ha="center", fontsize=14)
    figure.text(0.520, 0.51, "D", va="center", ha="center", fontsize=14)

    pyplot.savefig("./results/figures/Figure M5.pdf",
                   format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def draw_corr():
    with open("./results/data/step_4_repairability_multiple_errors.pkl", "rb") as file:
        task_1, task_2, task_3 = pload(file)

    figure = pyplot.figure(figsize=(10, 10), tight_layout=True)
    rcParams["font.family"] = "Linux Libertine"

    pyplot.subplot(2, 1, 1)
    numbers, rates = [], []
    gradient_colors = pyplot.get_cmap(name="rainbow")(linspace(0, 1, 12))
    x, y = [], []
    for filter_index in range(1, 13):
        x.append(sum(load("./results/data/a" + str(filter_index).zfill(2) + "[v].npy")))
        y.append(sum(task_1[filter_index, 2][:, -1]) / 2000.0)
        pyplot.scatter([log(float(x[-1])) / log(4)], [y[-1]],
                       color=gradient_colors[filter_index - 1], edgecolor="black", label=str(filter_index).zfill(2),
                       s=40, zorder=2)
    numbers = array(numbers)
    rates = array(rates)
    pyplot.scatter(log(numbers) / log(4), rates, color="silver", edgecolor="black")
    x, y = array(x), array(y)
    a, b = polyfit(log(x) / log(4), y, deg=1)
    shown_x = linspace(5, 10, 51)
    shown_y = a * shown_x + b
    pyplot.plot(shown_x, shown_y, color="black", linewidth=2, linestyle="--", zorder=1)
    pyplot.legend(loc="upper right", ncol=3, fontsize=14)
    pyplot.xlabel("valid number", fontsize=14)
    pyplot.xticks([5, 6, 7, 8, 9, 10],
                  ["4^5", "4^6", "4^7", "4^8", "4^9", "4^10"], fontsize=14)
    pyplot.xlim(4.95, 10.05)
    pyplot.ylabel("correction rate", fontsize=14)
    pyplot.ylim(-0.06, 1.06)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=14)

    pyplot.subplot(2, 1, 2)
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

    pyplot.xlabel("introduced error number", fontsize=14)
    pyplot.xticks(range(1, 17), range(1, 17), fontsize=14)
    pyplot.xlim(0.85, 16.15)
    pyplot.ylabel("correction rate", fontsize=14)
    pyplot.ylim(-0.06, 1.06)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["0%", "20%", "40%", "60%", "80%", "100%"], fontsize=14)

    figure.align_labels()
    figure.text(0.024, 0.99, "A", va="center", ha="center", fontsize=14)
    figure.text(0.024, 0.51, "B", va="center", ha="center", fontsize=14)

    pyplot.savefig("./results/figures/Figure S3.pdf",
                   format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()

    pyplot.figure(figsize=(10, 8), tight_layout=True)
    rcParams["font.family"] = "Linux Libertine"
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    values, rates = [[], []], [[], []]
    for filter_index in range(1, 13):
        indices = where(task_1[filter_index, 1][:, 2] == 1)[0]
        formers, latters = task_1[filter_index, 1][:, 3], task_1[filter_index, 1][:, 6]
        formers, latters = log10(formers[indices]), log10(latters[indices])

        violin_1 = pyplot.violinplot(formers, positions=[filter_index], widths=0.8, bw_method=0.5,
                                     showmeans=False, showextrema=False, showmedians=False)
        for body in violin_1["bodies"]:
            center = mean(body.get_paths()[0].vertices[:, 0])
            body.get_paths()[0].vertices[:, 0] = clip(body.get_paths()[0].vertices[:, 0], -inf, center)
            body.set_color(colors["algo1"])
            body.set_edgecolor("black")
            body.set_linewidth(1.5)
            body.set_alpha(1)
        pyplot.scatter([filter_index], [median(formers)], color=colors["algo1"], edgecolor="black", linewidth=1.5, s=40)
        values[0].append(mean(formers))

        violin_2 = pyplot.violinplot(latters, positions=[filter_index], widths=0.8, bw_method=0.5,
                                     showmeans=False, showextrema=False, showmedians=False)
        for body in violin_2["bodies"]:
            center = mean(body.get_paths()[0].vertices[:, 0])
            body.get_paths()[0].vertices[:, 0] = clip(body.get_paths()[0].vertices[:, 0], center, inf)
            body.set_color(colors["yyco1"])
            body.set_edgecolor("black")
            body.set_linewidth(1.5)
            body.set_alpha(1)

        pyplot.scatter([filter_index], [median(latters)], color=colors["yyco1"], edgecolor="black", linewidth=1.5, s=40)
        values[1].append(mean(latters))

        if filter_index == 1:
            pyplot.legend([violin_1["bodies"][0], violin_2["bodies"][0]], ["search-only", "combined"], fontsize=12)

    pyplot.xlabel("constraint set", fontsize=14)
    pyplot.xticks(range(1, 13), filter_indices, fontsize=14)
    pyplot.xlim(0.5, 12.5)
    pyplot.ylabel("candidate number", fontsize=14)
    pyplot.yticks([0, 1, 2, 3], ["1", "10", "100", "1000"], fontsize=14)
    pyplot.ylim(-0.2, 3.2)
    pyplot.savefig("./results/figures/Figure S4.pdf",
                   format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()

    show_single_examples()
    show_multiple_examples()

    multiple_evaluation(task_seed=2021, repeats=2000)

    draw_main()
    draw_corr()
