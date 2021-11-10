__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from matplotlib import pyplot
from numpy import log10, zeros_like, percentile, sum, min, max


colors = {
    "trad1": "#81B8DF",
    "trad2": "#B1CCDF",
}


def evaluate_practicability(data, errors):
    errors = log10(errors)

    violin = pyplot.violinplot(dataset=errors, positions=[0], showextrema=False, vert=False, widths=1)
    paint_data = violin["bodies"][0].get_paths()[0].vertices
    paint_data = paint_data[len(paint_data) // 2 + 1: -1]
    pyplot.close()

    x, y = paint_data[:, 0], paint_data[:, 1]
    y = y / sum(y) * 100.0

    pyplot.figure(figsize=(10, 4))
    pyplot.rc("font", family="Times New Roman")
    pyplot.subplot(2, 1, 2)

    pyplot.plot([-16, -2], [0, 0], color="black", linewidth=0.75)

    pyplot.plot(x, y, color=colors["trad1"], linewidth=2, zorder=1)
    pyplot.fill_between(x, y, zeros_like(y), color=colors["trad2"], zorder=0)

    [v_25, v_50, v_75] = percentile(a=errors, q=[25, 50, 75])

    pyplot.vlines(x=v_50, ymin=0, ymax=6.5, color="black", linewidth=0.75, zorder=2)
    pyplot.vlines(x=v_50, ymin=6.9, ymax=8.5, color="black", linewidth=0.75, zorder=2)
    pyplot.hlines(y=8.5, xmin=v_50 - 0.25, xmax=v_50 + 0.25, color="black", linewidth=0.75, zorder=2)
    pyplot.text(x=v_50, y=8.6, s="median", va="bottom", ha="center", zorder=2)

    pyplot.vlines(x=v_25, ymin=0, ymax=7, color="black", linewidth=0.75, zorder=2)
    pyplot.vlines(x=v_75, ymin=0, ymax=7, color="black", linewidth=0.75, zorder=2)
    pyplot.hlines(y=6.7, xmin=v_25, xmax=v_75, color="black", linewidth=0.75, zorder=2)
    pyplot.text(x=v_75 + 0.1, y=6.7, s="interquartile range", va="center", ha="left", zorder=2)

    lower, upper = v_25 - 1.5 * (v_75 - v_25), 1.5 * (v_75 - v_25) + v_75
    pyplot.vlines(x=lower, ymin=0, ymax=4, color="black", linewidth=0.75, zorder=2)
    pyplot.vlines(x=min(errors), ymin=0, ymax=4, color="black", linewidth=0.75, zorder=2)
    pyplot.hlines(y=3.7, xmin=lower, xmax=min(errors), color="black", linewidth=0.75, zorder=2)
    pyplot.text(x=min(errors) - 0.1, y=3.7, s="outlier range", va="center", ha="right", zorder=2)
    pyplot.vlines(x=upper, ymin=0, ymax=4, color="black", linewidth=0.75, zorder=2)
    pyplot.vlines(x=max(errors), ymin=0, ymax=4, color="black", linewidth=0.75, zorder=2)
    pyplot.hlines(y=3.7, xmin=upper, xmax=max(errors), color="black", linewidth=0.75, zorder=2)
    pyplot.text(x=max(errors) + 0.1, y=3.7, s="outlier range", va="center", ha="left", zorder=2)

    pyplot.xlim(-16, -2)
    pyplot.ylim(-0.5, 10.5)
    pyplot.xlabel("systematic error")
    pyplot.ylabel("frequency")
    pyplot.xticks([-15, -13, -11, -9, -7, -5, -3],
                  ["$10^{-15}$", "$10^{-13}$", "$10^{-11}$", "$10^{-9}$", "$10^{-7}$", "$10^{-5}$", "$10^{-3}$"])
    pyplot.yticks([0, 2, 4, 6, 8, 10], ["0%", "2%", "4%", "6%", "8%", "10%"])

    pyplot.subplot(2, 1, 1)
    pyplot.scatter(errors, data, color=colors["trad2"], edgecolor=colors["trad1"])

    pyplot.xlim(-16, -2)
    pyplot.ylim(-0.1, 2.1)
    pyplot.xticks([-15, -13, -11, -9, -7, -5, -3], ["", "", "", "", "", "", ""])
    pyplot.yticks([0, 0.5, 1.0, 1.5, 2.0], ["0.0  ", "0.5  ", "1.0  ", "1.5  ", "2.0  "])
    pyplot.ylabel("information density")

    pyplot.savefig("../results/errors.png", format="png",
                   bbox_inches="tight", transparent=True, dpi=600)
    pyplot.close()


if __name__ == "__main__":
    # obtained from TestRandom class in test/test_upper_bound.py
    results = [1.06E+00, 1.16E+00, 9.63E-01, 8.06E-01, 1.15E+00, 1.03E+00, 4.65E-01, 7.65E-01, 8.11E-01, 1.38E+00,
               4.06E-01, 1.02E+00, 1.16E+00, 1.00E+00, 1.11E+00, 1.17E+00, 6.61E-01, 1.09E+00, 1.09E+00, 9.42E-01,
               7.67E-01, 9.64E-01, 3.26E-05, 1.03E+00, 1.07E+00, 8.41E-01, 1.03E+00, 1.15E+00, 9.32E-01, 1.12E+00,
               1.09E+00, 8.39E-01, 9.89E-01, 8.22E-01, 1.15E+00, 1.19E+00, 9.66E-01, 6.61E-01, 1.00E+00, 8.99E-01,
               9.34E-01, 6.41E-01, 1.32E+00, 1.19E+00, 7.98E-01, 4.70E-01, 1.08E+00, 1.01E+00, 8.19E-01, 1.23E+00,
               1.07E+00, 1.11E+00, 1.20E+00, 6.94E-01, 1.02E+00, 1.23E+00, 1.08E+00, 1.15E+00, 8.57E-01, 3.17E-01,
               1.09E+00, 1.08E+00, 1.36E+00, 9.61E-01, 2.56E-01, 1.04E+00, 9.35E-01, 1.11E+00, 8.79E-01, 9.04E-01,
               1.10E+00, 1.00E+00, 9.96E-01, 1.02E+00, 1.44E-05, 8.87E-01, 9.87E-01, 9.83E-01, 1.07E+00, 9.99E-01,
               1.20E+00, 5.93E-01, 4.80E-01, 1.21E+00, 7.45E-01, 8.75E-01, 1.04E+00, 1.02E+00, 9.06E-01, 1.07E+00,
               1.04E+00, 1.09E+00, 6.94E-01, 7.37E-01, 1.19E+00, 9.04E-01, 1.02E+00, 1.18E+00, 1.01E+00, 9.18E-01]
    differs = [2.46E-11, 1.24E-11, 9.80E-11, 2.91E-11, 9.96E-12, 1.71E-11, 7.54E-10, 1.45E-10, 1.73E-10, 6.13E-12,
               8.49E-11, 3.33E-11, 4.26E-11, 6.66E-11, 4.57E-11, 5.80E-11, 1.07E-11, 3.84E-11, 3.84E-11, 4.05E-11,
               4.90E-11, 2.70E-11, 3.26E-05, 4.57E-11, 9.72E-11, 1.56E-11, 8.03E-11, 1.87E-11, 1.07E-11, 2.49E-11,
               3.34E-12, 2.56E-11, 2.12E-11, 7.67E-12, 9.21E-12, 1.52E-11, 2.20E-11, 3.25E-10, 3.47E-11, 3.19E-11,
               1.05E-10, 2.05E-10, 3.92E-11, 7.50E-12, 2.90E-10, 1.32E-11, 1.89E-11, 7.10E-11, 1.08E-11, 1.68E-11,
               1.66E-11, 5.75E-12, 1.55E-12, 3.04E-11, 4.56E-11, 7.99E-12, 1.16E-11, 1.44E-11, 8.54E-11, 3.49E-10,
               4.60E-13, 1.33E-12, 2.54E-13, 4.66E-12, 1.74E-10, 6.07E-11, 9.17E-11, 2.13E-11, 8.05E-11, 9.18E-12,
               2.19E-11, 4.71E-10, 1.46E-10, 3.21E-11, 1.44E-05, 9.18E-11, 2.53E-11, 2.08E-10, 1.06E-10, 3.18E-11,
               5.17E-11, 6.22E-11, 1.12E-11, 1.89E-11, 1.79E-10, 4.93E-11, 1.22E-10, 2.09E-11, 4.78E-11, 2.70E-11,
               2.02E-11, 2.42E-12, 1.16E-10, 2.85E-09, 1.48E-11, 8.83E-12, 1.48E-11, 3.10E-12, 2.89E-12, 1.11E-11]

    evaluate_practicability(data=results, errors=differs)
