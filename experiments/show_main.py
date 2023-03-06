from logging import getLogger, CRITICAL
from matplotlib import pyplot, rcParams, patches, lines
from numpy import zeros, array, random, linspace, arange, max, cumsum, deg2rad, sum, sin, cos, log10, where, nan
from warnings import filterwarnings

from experiments import load_data

filterwarnings("ignore")

getLogger("matplotlib").setLevel(CRITICAL)

rcParams["font.family"] = "Arial"
rcParams["mathtext.fontset"] = "custom"
rcParams["mathtext.rm"] = "Linux Libertine"
rcParams["mathtext.cal"] = "Lucida Calligraphy"
rcParams["mathtext.it"] = "Linux Libertine:italic"
rcParams["mathtext.bf"] = "Linux Libertine:bold"


# noinspection PyUnresolvedReferences
def main01():
    pyplot.figure(figsize=(10, 7), tight_layout=True)

    ax = pyplot.subplot(3, 1, 1)
    pyplot.text(0.03, 0.5, "generating process", color="#0070C0",
                va="center", ha="center", rotation="vertical", fontsize=12)
    pyplot.vlines(0.1, 0, 1, color="#0070C0", lw=1.5)
    ax.add_patch(patches.Circle(xy=(0.5, 0.5), radius=0.200, facecolor="#EEEEEE", edgecolor="#000000", lw=0.75))
    pyplot.text(0.5, 0.5, r"$\mathcal{D}_4^\ell$", va="center", ha="center", fontsize=14)
    pyplot.text(0.5, 0.1, "quaternary de Bruijn graph\nof order " + r"$\ell$",
                va="center", ha="center", fontsize=9)
    pyplot.text(1.0, 0.6, "screen vertices\nby constraint set " + r"$\mathrm{\mathbb{C}}$",
                va="center", ha="center", color="#0070C0", fontsize=9)
    pyplot.annotate("", xy=(1.2, 0.5), xytext=(0.8, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2))
    ax.add_patch(patches.Circle(xy=(1.5, 0.5), radius=0.175, facecolor="#FFE699", edgecolor="#FFC30E", lw=0.75))
    pyplot.text(1.5, 0.5, r"$\mathcal{D}_\mathrm{\mathbb{C}}^\ell$", va="center", ha="center", fontsize=14)
    pyplot.text(1.5, 0.1, "a screened subgraph\nof " + r"$\mathcal{D}_4^\ell$", va="center", ha="center", fontsize=9)
    pyplot.text(2.0, 0.6, "trim vertices\nthat outgoing arcs < 2",
                va="center", ha="center", color="#0070C0", fontsize=9)
    pyplot.annotate("", xy=(2.2, 0.5), xytext=(1.8, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2))
    ax.add_patch(patches.Circle(xy=(2.5, 0.5), radius=0.150, facecolor="#F6BE98", edgecolor="#EE833A", lw=0.75))
    pyplot.text(2.5, 0.1, "a trimmed subgraph\nof " + r"$\mathcal{D}_\mathrm{\mathbb{C}}^\ell$",
                va="center", ha="center", fontsize=9)
    pyplot.text(3.5, 0.9, "bind digits for arcs\nby a predetermined partial order",
                va="center", ha="center", color="#0070C0", fontsize=9)
    pyplot.annotate("", xy=(4.25, 0.6), xytext=(2.75, 0.6),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2, connectionstyle="arc3,rad=-0.3"))
    pyplot.hlines(0.69, 3.18, 3.82, lw=0.75)
    pyplot.hlines(0.62, 3.18, 3.82, lw=0.75)
    pyplot.hlines(0.07, 3.18, 3.82, lw=0.75)
    pyplot.vlines(3.50, 0.07, 0.69, lw=0.75)
    pyplot.text(3.35, 0.65, "nucleotide", va="center", ha="center", fontsize=8)
    pyplot.text(3.65, 0.65, "digit", va="center", ha="center", fontsize=8)
    color_groups = [["#FBE5D6", "#F4B183", "#ED7D31", "#A24A0E"], ["#E2F0D9", "#A9D18E", "#70AD47", "#4D7731"]]
    for index, bases in enumerate(["ACXX", "AXGX", "AXXT", "XCGX", "XCXT", "XXGT",
                                   "ACGX", "ACXT", "AXGT", "XCGT", "ACGT"]):
        digit = 0
        for location, base in enumerate(bases):
            if base != "X":
                pyplot.fill_between([3.5 - 0.04 - (3 - location) * 0.08 - 0.04,
                                     3.5 - 0.04 - (3 - location) * 0.08 + 0.04],
                                    0.59 - index * 0.05 - 0.02, 0.59 - index * 0.05 + 0.03,
                                    color=color_groups[0][{"A": 0, "C": 1, "G": 2, "T": 3}[base]], zorder=0, lw=0)
                if base != "T":
                    pyplot.text(3.5 - 0.04 - (3 - location) * 0.08, 0.59 - index * 0.05,
                                base, va="center", ha="center", fontsize=8)
                else:
                    pyplot.text(3.5 - 0.04 - (3 - location) * 0.08, 0.59 - index * 0.05,
                                base, va="center", ha="center", fontsize=8, color="white")
                pyplot.fill_between([3.5 + 0.04 + location * 0.08 - 0.04, 3.5 + 0.04 + location * 0.08 + 0.04],
                                    0.59 - index * 0.05 - 0.02, 0.59 - index * 0.05 + 0.03,
                                    color=color_groups[1][digit], zorder=0, lw=0)
                if digit < 3:
                    pyplot.text(3.5 + 0.04 + location * 0.08, 0.59 - index * 0.05,
                                str(digit), va="center", ha="center", fontsize=8, zorder=1)
                else:
                    pyplot.text(3.5 + 0.04 + location * 0.08, 0.59 - index * 0.05,
                                str(digit), va="center", ha="center", fontsize=8, color="white", zorder=1)
                digit += 1
    ax.add_patch(patches.Circle(xy=(4.5, 0.5), radius=0.15, facecolor="#FCBBAE", edgecolor="#ff8D8A", lw=0.75))
    pyplot.text(4.5, 0.485, r"$\mathcal{D}$", va="center", ha="center", fontsize=14)
    pyplot.text(4.5, 0.1, "coding digraph\nof " + r"$\ell$" + " and " + r"$\mathrm{\mathbb{C}}$",
                va="center", ha="center", fontsize=9)
    pyplot.axis("off")
    pyplot.xlim(0, 5)
    pyplot.ylim(0, 1)

    pyplot.subplot(3, 1, 2)
    pyplot.text(0.03, 0.5, "coding process", color="#0070C0",
                va="center", ha="center", rotation="vertical", fontsize=12)
    pyplot.vlines(0.1, 0, 1, color="#0070C0", lw=1.5)
    # pyplot.text(3.6, 0.9, "graph-based encoding", va="center", ha="center", color="#0070C0", fontsize=12)
    pyplot.text(0.5, 0.92, "part of coding digraph", va="center", ha="center", fontsize=10)
    points = array([[0.5, 0.2, 0.4, 0.6, 0.8, 0.5, 0.7, 0.4, 0.6, 0.4, 0.6, 0.8],
                    [0.9, 0.7, 0.7, 0.7, 0.7, 0.5, 0.5, 0.3, 0.3, 0.1, 0.1, 0.1]])
    points[1] -= 0.05
    pyplot.scatter(points[0], points[1], color="white", edgecolor="black", lw=0.75, s=30)
    for former, latter, info, mark, bias in [(0, 1, "{A|0}", 0, -0.07), (0, 2, "{C|1}", 0, -0.03),
                                             (0, 3, "{G|2}", 1, +0.03), (0, 4, "{T|3}", 0, +0.07),
                                             (3, 5, "{A|0}", 1, -0.03), (3, 6, "{G|1}", 0, +0.03),
                                             (5, 7, "{A|0}", 0, -0.03), (5, 8, "{C|1}", 1, +0.03),
                                             (8, 9, "{C|0}", 0, -0.05), (8, 10, "{G|1}", 0, 0),
                                             (8, 11, "{T|2}", 1, + 0.05)]:
        point_1, point_2 = (points[0][former], points[1][former]), (points[0][latter], points[1][latter])
        connection = None
        if points[0][former] < points[0][latter]:
            connection = "arc3,rad=0.2"
        elif points[0][former] > points[0][latter]:
            connection = "arc3,rad=-0.2"
        if mark:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="#000000", connectionstyle=connection,
                                            shrinkA=3, shrinkB=5, lw=0.75))
            x, y = (point_1[0] + point_2[0]) / 2 + bias, (point_1[1] + point_2[1]) / 2 + 0.01
            pyplot.plot([x - 0.05, x - 0.05, x + 0.05, x + 0.05, x - 0.05],
                        [y - 0.025, y + 0.025, y + 0.025, y - 0.025, y - 0.025], color="k", lw=0.75)
            pyplot.fill_between([x - 0.05, x], y - 0.025, y + 0.025,
                                fc=color_groups[0][{"A": 0, "C": 1, "G": 2, "T": 3}[info[1]]], lw=0, zorder=2)
            if info[1] != "T":
                pyplot.text(x - 0.025, y - 0.005, info[1], color="k", va="center", ha="center", fontsize=7.5, zorder=3)
            else:
                pyplot.text(x - 0.025, y - 0.005, info[1], color="w", va="center", ha="center", fontsize=7.5, zorder=3)
            pyplot.fill_between([x, x + 0.05], y - 0.025, y + 0.025, fc=color_groups[1][int(info[3])], lw=0, zorder=2)
            if info[3] != "3":
                pyplot.text(x + 0.025, y - 0.005, info[3], color="k", va="center", ha="center", fontsize=7.5, zorder=3)
            else:
                pyplot.text(x + 0.025, y - 0.005, info[3], color="w", va="center", ha="center", fontsize=7.5, zorder=3)
        else:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="#666666", connectionstyle=connection,
                                            shrinkA=3, shrinkB=5, lw=0.75))
    for former, latter in [(1, 2), (2, 3), (3, 4), (5, 6), (7, 8), (9, 10), (10, 11)]:
        pyplot.hlines(points[1][former], points[0][former] + 0.05, points[0][latter] - 0.05, lw=0.75, ls="--")
    pyplot.text(1.2, 0.75, "out-degree", va="center", ha="center", fontsize=10)
    pyplot.hlines(points[1][4], points[0][4] + 0.05, 1.15, lw=0.75, ls="--")
    pyplot.text(1.2, points[1][4], "4", va="center", ha="center", fontsize=9)
    pyplot.hlines(points[1][6], points[0][6] + 0.05, 1.15, lw=0.75, ls="--")
    pyplot.text(1.2, points[1][6], "2", va="center", ha="center", fontsize=9)
    pyplot.hlines(points[1][8], points[0][8] + 0.05, 1.15, lw=0.75, ls="--")
    pyplot.text(1.2, points[1][8], "2", va="center", ha="center", fontsize=9)
    pyplot.hlines(points[1][11], points[0][11] + 0.05, 1.15, lw=0.75, ls="--")
    pyplot.text(1.2, points[1][11], "3", va="center", ha="center", fontsize=9)
    pyplot.fill_between([1.45, 1.55], 0.16, 0.55, facecolor="w", edgecolor="k", linewidth=0.75)
    pyplot.fill_between([2.03, 2.17], 0.46, 0.55, facecolor="w", edgecolor="k", linewidth=0.75)
    pyplot.fill_between([2.35, 2.45], 0.16, 0.55, facecolor="w", edgecolor="k", linewidth=0.75)
    for index, info in enumerate([("$4$", r"$\times$", "$10$", r"$=$", "$42$", r"$-$", "$2$"),
                                  ("$2$", r"$\times$", "$5$", r"$=$", "$10$", r"$-$", "$0$"),
                                  ("$2$", r"$\times$", "$2$", r"$=$", "$5$", r"$-$", "$1$"),
                                  ("$3$", r"$\times$", "$0$", r"$=$", "$2$", r"$-$", "$2$")]):
        for location, string in enumerate(info):
            pyplot.text(1.5 + location * 0.15, 0.5 - index * 0.1, string, va="center", ha="center", fontsize=10)
    pyplot.annotate("", xy=(1.25, points[1][4]), xytext=(1.45, 0.5),
                    arrowprops=dict(arrowstyle="-", color="k", lw=0.75, ls="--",
                                    connectionstyle="arc3,rad=0.1"))
    pyplot.annotate("", xy=(1.25, points[1][6]), xytext=(1.45, 0.4),
                    arrowprops=dict(arrowstyle="-", color="k", lw=0.75, ls="--",
                                    connectionstyle="arc3,rad=0.05"))
    pyplot.annotate("", xy=(1.25, points[1][8]), xytext=(1.45, 0.3),
                    arrowprops=dict(arrowstyle="-", color="k", lw=0.75, ls="--",
                                    connectionstyle="arc3,rad=-0.05"))
    pyplot.annotate("", xy=(1.25, points[1][11]), xytext=(1.45, 0.2),
                    arrowprops=dict(arrowstyle="-", color="k", lw=0.75, ls="--",
                                    connectionstyle="arc3,rad=-0.1"))
    pyplot.text(3.6, 0.8, "[101010]", va="center", ha="center", fontsize=10)
    pyplot.text(3.8, 0.8, "binary message", va="center", ha="left", fontsize=10)
    pyplot.text(2.55, 0.65, r"$2^1+2^3+2^5$", va="center", ha="right", fontsize=10)
    pyplot.hlines(0.65, 2.60, 2.88, lw=0.75, ls="--")
    pyplot.text(3.55, 0.65, "decimal conversion", va="center", ha="right", color="#0070C0", fontsize=10)
    pyplot.annotate("", xy=(3.6, 0.55), xytext=(3.6, 0.75),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2))
    pyplot.text(3.6, 0.5, "42", va="center", ha="center", fontsize=10)
    pyplot.text(3.8, 0.5, "decimal number", va="center", ha="left", fontsize=10)
    pyplot.hlines(0.35, 2.45, 2.73, lw=0.75, ls="--")
    pyplot.text(3.55, 0.35, "graph-based conversion", va="center", ha="right", color="#0070C0", fontsize=10)
    pyplot.annotate("", xy=(3.6, 0.25), xytext=(3.6, 0.45),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2))
    pyplot.text(3.6, 0.20, "[2012]", va="center", ha="center", fontsize=10)
    pyplot.text(3.8, 0.20, "graph-based vector", va="center", ha="left", fontsize=10)
    pyplot.text(3.6, 0.12, r"$=$", va="center", ha="center", rotation="vertical", fontsize=10)
    pyplot.text(3.6, 0.04, "[GACT]", va="center", ha="center", fontsize=10)
    pyplot.text(3.8, 0.04, "DNA sequence", va="center", ha="left", fontsize=10)
    pyplot.annotate("", xy=(4.45, 0.8), xytext=(4.45, 0.04),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2, connectionstyle="arc3,rad=0.2"))
    pyplot.text(4.8, 0.42, "decode by\nreverse steps", va="center", ha="center", color="#0070C0", fontsize=10)
    pyplot.axis("off")
    pyplot.xlim(0, 5)
    pyplot.ylim(0, 1)

    pyplot.subplot(3, 1, 3)
    pyplot.text(0.03, 0.5, "error correcting process", color="#0070C0",
                va="center", ha="center", rotation="vertical", fontsize=12)
    pyplot.vlines(0.1, 0, 1, color="#0070C0", lw=1.5)
    pyplot.fill_between([0.15, 3.50], 0.05, 0.95, color="#EFEFEF", lw=0)
    # pyplot.text(3.0, 0.85, "path-based error correcting", va="center", ha="center", color="#0070C0", fontsize=12)
    # pyplot.text(4.6, 0.50, "pretreatment-free\nretrieval mechanism",
    #             va="center", ha="center", color="#0070C0", fontsize=12)
    points = array([[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 0.4, 0.6, 0.8, 1.0, 1.2, 0.4, 0.6, 0.8, 1.0, 0.4, 0.6, 0.8],
                    [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.55, 0.55, 0.55, 0.55, 0.55, 0.4, 0.4, 0.4, 0.4, 0.25, 0.25, 0.25]])
    pyplot.scatter(points[0], points[1], color="white", edgecolor="black", lw=0.75, s=30, zorder=2)
    pyplot.annotate("", xy=(1.1, 0.7), xytext=(1.1, 0.83),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2))
    pyplot.text(1.1, 0.85, "find error", va="center", ha="center", color="#0070C0", fontsize=9)
    for former, latter, mark in [(0, 1, 0), (1, 2, 0), (2, 3, 0), (3, 4, 0), (4, 5, 1),
                                 (0, 6, 0), (6, 7, 0), (7, 8, 0), (8, 9, 1), (9, 10, 0),
                                 (0, 11, 0), (11, 12, 0), (12, 13, 1), (13, 14, 0),
                                 (0, 15, 0), (15, 16, 1), (16, 17, 0)]:
        point_1, point_2 = (points[0][former], points[1][former]), (points[0][latter], points[1][latter])
        if latter == former + 1:
            if mark:
                pyplot.fill_between([point_1[0] - 0.05, point_2[0] + 0.05], point_1[1] - 0.05, point_1[1] + 0.05,
                                    facecolor="#FCBBAE", edgecolor="#ff8D8A", linewidth=0.75, zorder=0)
                pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                                arrowprops=dict(arrowstyle="<|-", color="red", shrinkA=3, shrinkB=4, lw=0.75))
            else:
                pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                                arrowprops=dict(arrowstyle="<|-", color="k", shrinkA=3, shrinkB=4, lw=0.75))
        else:
            if mark:
                pyplot.fill_between([point_1[0] - 0.05, point_2[0] + 0.05], point_1[1] - 0.05, point_1[1] + 0.05,
                                    facecolor="#FCBBAE", edgecolor="#ff8D8A", linewidth=0.75, zorder=0)
                pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                                arrowprops=dict(arrowstyle="<|-", color="red", connectionstyle="arc3,rad=-0.2",
                                                shrinkA=3, shrinkB=4, lw=0.75))
            else:
                pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                                arrowprops=dict(arrowstyle="<|-", color="k", connectionstyle="arc3,rad=-0.2",
                                                shrinkA=3, shrinkB=4, lw=0.75))
    pyplot.annotate("", xy=(0.85, 0.25), xytext=(1.1, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2, connectionstyle="arc3,rad=-0.3"))
    pyplot.text(0.9, 0.12, "local exhaustive\nreverse search", va="center", ha="center", color="#0070C0", fontsize=9)
    pyplot.hlines(0.70, 1.30, 1.35, color="k", lw=0.75, ls="--")
    pyplot.plot([1.40, 1.35, 1.35, 1.40], [0.80, 0.80, 0.08, 0.08], color="k", lw=0.75, ls="--")
    points = array([[1.50, 1.65, 2.40, 1.65, 1.80, 1.95, 2.10, 1.65, 1.80, 1.95],
                    [0.80, 0.80, 0.80, 0.72, 0.72, 0.72, 0.72, 0.64, 0.64, 0.64]])
    pyplot.scatter(points[0], points[1], color="white", edgecolor="black", lw=0.75, s=30, zorder=2)
    pyplot.fill_between([points[0][0] - 0.05, points[0][1] + 0.05], points[1][0] - 0.05, points[1][0] + 0.05,
                        facecolor="#FCBBAE", edgecolor="#FF8D8A", linewidth=0.75, zorder=1)
    pyplot.text((points[0][1] + points[0][2]) / 2, (points[1][1] + points[1][2]) / 2 - 0.01, r"$\cdots$", fontsize=9,
                va="center", ha="center", bbox=dict(boxstyle="round", fc="#EFEFEF", lw=0, alpha=1), zorder=2)
    pyplot.text(points[0][1] + 0.07, points[1][0], "substitution", va="bottom", ha="left", fontsize=9, zorder=3)
    for former, latter in [(0, 1), (1, 2), (0, 3), (3, 4), (4, 5), (5, 6), (6, 2), (0, 7), (7, 8), (8, 9)]:
        point_1, point_2 = (points[0][former], points[1][former]), (points[0][latter], points[1][latter])
        if (former, latter) in [(0, 1), (8, 9)]:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="red", shrinkA=3, shrinkB=4, lw=0.75))
        elif latter != former + 1:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="k", connectionstyle="arc3,rad=-0.2",
                                            shrinkA=3, shrinkB=4, lw=0.75))
        else:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="k", shrinkA=3, shrinkB=4, lw=0.75))
    points = array([[1.50, 1.80, 2.40, 1.65, 1.80, 1.95, 2.10, 1.65, 1.80, 1.95],
                    [0.52, 0.52, 0.52, 0.44, 0.44, 0.44, 0.44, 0.36, 0.36, 0.36]])
    pyplot.scatter(points[0], points[1], color="white", edgecolor="black", lw=0.75, s=30, zorder=2)
    pyplot.fill_between([points[0][0] - 0.05, points[0][1] + 0.05], points[1][0] - 0.05, points[1][0] + 0.05,
                        facecolor="#FCBBAE", edgecolor="#FF8D8A", linewidth=0.75, zorder=1)
    pyplot.text((points[0][1] + points[0][2]) / 2, (points[1][1] + points[1][2]) / 2 - 0.01, r"$\cdots$", fontsize=9,
                va="center", ha="center", bbox=dict(boxstyle="round", fc="#EFEFEF", lw=0, alpha=1), zorder=2)
    pyplot.text(points[0][1] + 0.07, points[1][0], "insertion", va="bottom", ha="left", fontsize=9, zorder=3)
    for former, latter in [(0, 1), (1, 2), (0, 3), (3, 4), (4, 5), (5, 6), (6, 2), (0, 7), (7, 8), (8, 9)]:
        point_1, point_2 = (points[0][former], points[1][former]), (points[0][latter], points[1][latter])
        if (former, latter) in [(0, 1), (8, 9)]:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="red", shrinkA=3, shrinkB=4, lw=0.75))
        elif latter != former + 1:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="k", connectionstyle="arc3,rad=-0.2",
                                            shrinkA=3, shrinkB=4, lw=0.75))
        else:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="k", shrinkA=3, shrinkB=4, lw=0.75))
    points = array([[1.50, 1.65, 1.80, 2.40, 1.80, 1.95, 2.10, 2.25, 1.80, 1.95, 2.10],
                    [0.24, 0.24, 0.24, 0.24, 0.16, 0.16, 0.16, 0.16, 0.08, 0.08, 0.08]])
    pyplot.scatter(points[0], points[1], color="white", edgecolor="black", lw=0.75, s=30, zorder=2)
    pyplot.fill_between([points[0][0] - 0.05, points[0][2] + 0.05], points[1][0] - 0.05, points[1][0] + 0.05,
                        facecolor="#FCBBAE", edgecolor="#FF8D8A", linewidth=0.75, zorder=1)
    pyplot.text((points[0][2] + points[0][3]) / 2, (points[1][2] + points[1][3]) / 2 - 0.01, r"$\cdots$", fontsize=9,
                va="center", ha="center", bbox=dict(boxstyle="round", fc="#EFEFEF", lw=0, alpha=1), zorder=2)
    pyplot.text(points[0][2] + 0.07, points[1][0], "deletion", va="bottom", ha="left", fontsize=9, zorder=3)
    for former, latter in [(0, 1), (1, 2), (2, 3), (0, 4), (4, 5), (5, 6), (6, 7), (7, 3), (0, 8), (8, 9), (9, 10)]:
        point_1, point_2 = (points[0][former], points[1][former]), (points[0][latter], points[1][latter])
        if (former, latter) in [(0, 1), (9, 10)]:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="red", shrinkA=3, shrinkB=4, lw=0.75))
        elif latter != former + 1:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="k", connectionstyle="arc3,rad=-0.2",
                                            shrinkA=3, shrinkB=4, lw=0.75))
        else:
            pyplot.annotate(s="", xy=point_1, xytext=point_2, zorder=1,
                            arrowprops=dict(arrowstyle="<|-", color="k", shrinkA=3, shrinkB=4, lw=0.75))
    pyplot.annotate("", xy=(2.75, 0.5), xytext=(2.45, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2))
    pyplot.text(2.6, 0.6, "generate\nVT-check", va="center", ha="center", color="#0070C0", fontsize=9)
    random.seed(2021)
    orders = arange(8)
    random.shuffle(orders)
    for index, order in enumerate(orders):
        pyplot.fill_between([2.80, 3.00], 0.66 - index * 0.08, 0.70 - index * 0.08, fc="w", ec="k", lw=0.75)
        pyplot.hlines(0.68 - index * 0.08, 3.00, 3.07, color="k", lw=0.75)
        pyplot.fill_between([3.07, 3.15], 0.66 - index * 0.08, 0.70 - index * 0.08, fc="w", ec="k", lw=0.75)
        if order != 0:
            pyplot.text(3.18, 0.675 - index * 0.08, r"$\times$", color="r", va="center", ha="center", fontsize=12)
            x = random.uniform(low=2.82, high=2.98, size=(random.randint(1, 3),))
            y = zeros(shape=(len(x),)) + 0.68 - index * 0.08
            pyplot.scatter(x, y, color="r", s=10, zorder=2)
        else:
            pyplot.text(3.35, 0.8 - index * 0.08, "output\nsolution\ncandidate",
                        va="center", ha="center", color="#0070C0", fontsize=9)
            pyplot.annotate("", xy=(3.5, 0.68 - index * 0.08), xytext=(3.2, 0.68 - index * 0.08),
                            arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2))
    curve = (1.0 / linspace(1, 2, 20) - 0.5) * 0.9
    pyplot.fill_between(linspace(3.5, 3.7, 20), 0.5 - curve, 0.5 + curve, color="#EFEFEF", lw=0)
    pyplot.annotate("", xy=(3.8, 0.3), xytext=(3.8, 0.7), zorder=1,
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2, connectionstyle="arc3,rad=0.45"))
    pyplot.scatter([3.7], [0.5], s=50, color="#EFEFEF", edgecolor="k", lw=0.75, zorder=2)
    pyplot.text(4, 0.4, "solution\ncandidate", va="center", ha="center", fontsize=9)
    pyplot.fill_between([3.85, 4.15], 0.28, 0.32, facecolor="w", edgecolor="k", lw=0.75)
    pyplot.text(4, 0.8, "raw DNA\nsequence", va="center", ha="center", fontsize=9)
    pyplot.fill_between([3.85, 4.15], 0.68, 0.72, facecolor="w", edgecolor="k", lw=0.75)
    pyplot.scatter([3.90, 4.03], [0.7, 0.7], color="r", s=10, zorder=2)

    pyplot.annotate("", xy=(4.4, 0.22), xytext=(4.2, 0.3), zorder=1,
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2, connectionstyle="arc3,rad=-0.2"))
    pyplot.text(4.42, 0.31, "count", va="center", ha="center", color="#0070C0", fontsize=10)
    pyplot.text(4.27, 0.16, "frequency", va="center", ha="right", fontsize=9)
    pyplot.plot([4.3, 4.3, 4.9], [0.24, 0.08, 0.08], lw=0.75, color="k")
    for index, high in enumerate([1.0, 1.0, 1.0, 0.9, 0.8, 0.2, 0.1, 0.1, 0.1, 0.1]):
        pyplot.fill_between([4.32 + index * 0.06, 4.36 + index * 0.06], 0.1, 0.1 + high * 0.1,
                            fc="lightblue", ec="k", lw=0.75)
    pyplot.hlines(0.06, 4.32, 4.60, color="#0070C0", lw=1.5)
    pyplot.text(4.46, 0.02, "retrieve", va="center", color="#0070C0", ha="center", fontsize=10)
    pyplot.axis("off")
    pyplot.xlim(0, 5)
    pyplot.ylim(0, 1)

    pyplot.savefig("./show/main01.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


# noinspection PyUnresolvedReferences
def main02():
    figure = pyplot.figure(figsize=(10, 6), tight_layout=True)
    grid = pyplot.GridSpec(2, 2)

    pyplot.subplot(grid[0, :])

    def draw_graph(ps, s, w, se=None, h=None):
        for i, (xi, yi) in enumerate(ps):
            pyplot.scatter(xi, yi, color="white", edgecolor="black", s=s, lw=w)
        for f, l in [(1, 2), (2, 4), (2, 5), (3, 1), (3, 8), (4, 3), (5, 6), (6, 1), (6, 8), (7, 4), (7, 5), (8, 7)]:
            if se is None:
                pyplot.annotate(s="", xy=tuple(ps[l - 1]), xytext=tuple(ps[f - 1]),
                                arrowprops=dict(arrowstyle="-|>", color="black",
                                                shrinkA=s / 5, shrinkB=s / 5, lw=w), zorder=1)
            elif (f, l) in se:
                pyplot.annotate(s="", xy=tuple(ps[l - 1]), xytext=tuple(ps[f - 1]),
                                arrowprops=dict(arrowstyle="-|>", color=("black" if (f, l) not in h else "red"),
                                                shrinkA=s / 5, shrinkB=s / 5, lw=w), zorder=1)
        for f, l in [(1, 3), (5, 7)]:
            if se is None:
                pyplot.annotate(s="", xy=tuple(ps[f - 1]), xytext=tuple(ps[l - 1]), zorder=1,
                                arrowprops=dict(arrowstyle="<|-", color="black",
                                                connectionstyle="arc3,rad=0.9", shrinkA=s / 5, shrinkB=s / 5, lw=w))
            elif (f, l) in se:
                pyplot.annotate(s="", xy=tuple(ps[f - 1]), xytext=tuple(ps[l - 1]), zorder=1,
                                arrowprops=dict(arrowstyle="<|-", color="black" if (f, l) not in h else "red",
                                                connectionstyle="arc3,rad=0.9", shrinkA=s / 5, shrinkB=s / 5, lw=w))
        for f, l in [(4, 2), (8, 6)]:
            if se is None:
                pyplot.annotate(s="", xy=tuple(ps[l - 1]), xytext=tuple(ps[f - 1]), zorder=1,
                                arrowprops=dict(arrowstyle="-|>", color="black",
                                                connectionstyle="arc3,rad=0.9", shrinkA=s / 5, shrinkB=s / 5, lw=w))
            elif (f, l) in se:
                pyplot.annotate(s="", xy=tuple(ps[l - 1]), xytext=tuple(ps[f - 1]), zorder=1,
                                arrowprops=dict(arrowstyle="-|>", color="black" if (f, l) not in h else "red",
                                                connectionstyle="arc3,rad=0.9", shrinkA=s / 5, shrinkB=s / 5, lw=w))

    length, bias, points = 0.3, 22.5, []
    for angle in arange(0, 360, 45):
        radio = deg2rad(angle + bias)
        points.append([sin(radio) * length + 0.45, cos(radio) * length + 0.5])
    points = array(points)
    pyplot.text(0.45, 0.85, "GC-balanced coding digraph", va="center", ha="center", fontsize=10)
    draw_graph(ps=points, s=20, w=1)
    color_groups = [["#FBE5D6", "#F4B183", "#ED7D31", "#A24A0E"], ["#E2F0D9", "#A9D18E", "#70AD47", "#4D7731"]]
    x, y = 0.45, 0.5
    for location, info in zip([(x + 0.38, y + 0.15), (x - 0.38, y + 0.15), (x + 0.38, y - 0.15), (x - 0.38, y - 0.15),
                               (x + 0.10, y + 0.07), (x + 0.10, y - 0.07), (x - 0.10, y + 0.07), (x - 0.10, y - 0.07),
                               (x + 0.22, y + 0.21), (x - 0.22, y + 0.21), (x + 0.22, y - 0.21), (x - 0.22, y - 0.21),
                               (x + 0.15, y + 0.14), (x + 0.15, y - 0.14), (x - 0.15, y + 0.14), (x - 0.15, y - 0.14)],
                              ["{G|1}", "{G|1}", "{C|0}", "{C|0}", "{T|1}", "{T|1}", "{A|0}", "{A|0}",
                               "{C|0}", "{C|0}", "{G|1}", "{G|1}", "{A|0}", "{A|0}", "{T|1}", "{T|1}"]):
        xs = [location[0] - 0.050, location[0] - 0.050, location[0] + 0.050, location[0] + 0.050, location[0] - 0.050]
        ys = [location[1] - 0.025, location[1] + 0.025, location[1] + 0.025, location[1] - 0.025, location[1] - 0.025]
        pyplot.plot(xs, ys, color="k", lw=0.75)
        pyplot.vlines(location[0], location[1] - 0.025, location[1] + 0.025, lw=0.75, zorder=3)
        pyplot.fill_between([location[0] - 0.05, location[0]], location[1] - 0.025, location[1] + 0.025,
                            fc=color_groups[0][{"A": 0, "C": 1, "G": 2, "T": 3}[info[1]]], lw=0, zorder=2)
        if info[1] != "T":
            pyplot.text(location[0] - 0.025, location[1] - 0.005, info[1],
                        color="k", va="center", ha="center", fontsize=7.5, zorder=3)
        else:
            pyplot.text(location[0] - 0.025, location[1] - 0.005, info[1],
                        color="w", va="center", ha="center", fontsize=7.5, zorder=3)
        pyplot.fill_between([location[0], location[0] + 0.05], location[1] - 0.025, location[1] + 0.025,
                            fc=color_groups[1][int(info[3])], lw=0, zorder=2)
        if info[3] != "3":
            pyplot.text(location[0] + 0.025, location[1] - 0.005, info[3],
                        color="k", va="center", ha="center", fontsize=7.5, zorder=3)
        else:
            pyplot.text(location[0] + 0.025, location[1] - 0.005, info[3],
                        color="w", va="center", ha="center", fontsize=7.5, zorder=3)

    pyplot.scatter([0.1], [0.1], color="white", edgecolor="black", s=20, lw=1)
    pyplot.text(0.13, 0.1, "vertex", va="center", ha="left", fontsize=8)
    pyplot.annotate(s="", xy=(0.3, 0.14), xytext=(0.4, 0.14), zorder=1,
                    arrowprops=dict(arrowstyle="<|-", color="black", lw=1))
    pyplot.annotate(s="", xy=(0.3, 0.06), xytext=(0.4, 0.06), zorder=1,
                    arrowprops=dict(arrowstyle="<|-", color="red", lw=1))
    pyplot.text(0.4, 0.14, "normal arc", va="center", ha="left", fontsize=8)
    pyplot.text(0.4, 0.06, "access arc (error found)", va="center", ha="left", fontsize=8)
    pyplot.fill_between([0.7 - 0.017, 0.7 + 0.017], 0.14 - 0.015, 0.14 + 0.018, lw=0, color="lightblue")
    pyplot.text(0.7, 0.14, "N", va="center", ha="center", fontsize=8)
    pyplot.text(0.73, 0.14, "adjust (substitute/insert/delete) by N", va="center", ha="left", fontsize=8)
    pyplot.text(0.95, 0.06, "N", va="center", ha="center", color="r", fontsize=8)
    pyplot.text(0.98, 0.06, "nucleotide error", va="center", ha="left", fontsize=8)

    pyplot.text(1.0, 0.58, "produce a\nsequence", va="center", ha="center", color="#0070C0", fontsize=9)
    pyplot.annotate("", xy=(1.15, 0.5), xytext=(0.85, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2))
    length = 0.1
    for center, select, highlight, sequence in zip([(1.25, 0.35), (1.6, 0.35)],
                                                   [[(3, 8), (8, 7), (7, 5), (5, 6), (6, 1), (1, 2)],
                                                    [(3, 8), (8, 7), (7, 5), (5, 6), (5, 7)]],
                                                   [[], [(5, 6), (5, 7)]],
                                                   ["TCTGAC", "TCTAAC"]):
        points = []
        for angle in arange(0, 360, 45):
            radio = deg2rad(angle + bias)
            points.append([sin(radio) * length + center[0], cos(radio) * length + center[1]])
        points = array(points)
        draw_graph(ps=points, s=10, w=0.5, se=select, h=highlight)
        start = center[0] - 0.03 * 2.5
        for index, nucleotide in enumerate(sequence):
            if index == 3 and nucleotide == "A":
                pyplot.text(start + index * 0.03, 0.5, nucleotide, color="r", va="center", ha="center", fontsize=8)
            else:
                pyplot.text(start + index * 0.03, 0.5, nucleotide, va="center", ha="center", fontsize=8)
    pyplot.text(1.53, 0.69, "find error\nand correct", va="center", ha="center", color="#0070C0", fontsize=9)
    pyplot.annotate("", xy=(1.75, 0.675), xytext=(1.6, 0.55),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2, connectionstyle="arc3,rad=-0.2"))
    pyplot.text(1.57, 1.1, "G", va="center", ha="center", fontsize=8)
    pyplot.text(1.60, 1.1, "C", va="center", ha="center", fontsize=8)
    pyplot.text(1.63, 1.1, "T", va="center", ha="center", fontsize=8)
    pyplot.annotate("", xy=(1.54, 1.1), xytext=(1.25, 0.55),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2, connectionstyle="arc3,rad=-0.3"))
    pyplot.text(1.2, 1.0, "generate\nVT-check",
                va="center", ha="center", color="#0070C0", fontsize=9)
    pyplot.annotate("", xy=(3.73, 1.03), xytext=(1.65, 1.1),
                    arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2, connectionstyle="arc3,rad=-0.03"))
    pyplot.text(3.55, 1.1, "sieve solutions", va="center", ha="center", color="#0070C0", fontsize=9)
    pyplot.fill_between([1.75, 4.00], 0, 1.0, color="#EEEEEE", zorder=0)
    pyplot.text(2.10, 0.95, "exhaustive", va="center", ha="center", fontsize=9)
    pyplot.text(1.88, 0.82, "local\nand\nreverse", va="center", ha="center", fontsize=9)
    pyplot.annotate("", xy=(2.0, 0.9), xytext=(2.2, 0.9),
                    arrowprops=dict(arrowstyle="<|-", color="black", shrinkA=0, shrinkB=0, lw=0.75))
    pyplot.annotate("", xy=(2.0, 0.9), xytext=(2.0, 0.73),
                    arrowprops=dict(arrowstyle="<|-", color="black", shrinkA=0, shrinkB=0, lw=0.75))
    pyplot.hlines(0.86, 2.15, 2.65, color="black", lw=0.75)
    pyplot.hlines(0.86, 2.75, 3.25, color="black", lw=0.75)
    pyplot.hlines(0.86, 3.35, 3.55, color="black", lw=0.75)
    pyplot.text(2.40, 0.88, "substitute", va="center", ha="center", fontsize=9)
    pyplot.text(3.00, 0.88, "insert", va="center", ha="center", fontsize=9)
    pyplot.text(3.45, 0.88, "delete", va="center", ha="center", fontsize=9)
    pyplot.vlines(2.08, 0.575, 0.775, lw=0.75)
    pyplot.vlines(2.08, 0.310, 0.510, lw=0.75)
    pyplot.vlines(2.08, 0.025, 0.225, lw=0.75)
    pyplot.text(2.06, 0.675, "(+0) position", va="center", ha="right", fontsize=9)
    pyplot.text(2.06, 0.410, "(\N{MINUS SIGN}1) position", va="center", ha="right", fontsize=9)
    pyplot.text(2.06, 0.125, "(\N{MINUS SIGN}2) position", va="center", ha="right", fontsize=9)
    for center, select, highlight, sequence, change in zip([(2.25, 0.675), (2.25, 0.410), (2.25, 0.125), (2.55, 0.675),
                                                            (2.85, 0.675), (2.85, 0.410), (2.85, 0.125),
                                                            (3.15, 0.675), (3.15, 0.410), (3.15, 0.125),
                                                            (3.45, 0.675), (3.45, 0.410), (3.45, 0.125)],
                                                           [[(3, 8), (8, 7), (7, 5), (5, 7), (7, 4), (4, 2)],
                                                            [(3, 8), (8, 7), (7, 4), (4, 3), (4, 2)],
                                                            [(3, 8), (8, 6), (6, 8), (8, 7)],
                                                            [(3, 8), (8, 7), (7, 5), (5, 6), (6, 1), (1, 2)],
                                                            [(3, 8), (8, 7), (7, 5), (5, 7), (7, 4), (4, 3), (4, 2)],
                                                            [(3, 8), (8, 7), (7, 4), (4, 3), (4, 2)],
                                                            [(3, 8), (8, 7), (7, 5), (7, 4)],
                                                            [(3, 8), (8, 7), (7, 5), (5, 6), (6, 1), (1, 2), (1, 3)],
                                                            [(3, 8), (8, 7), (7, 5), (5, 6), (5, 7)],
                                                            [(3, 8), (8, 6), (6, 1), (6, 8)],
                                                            [(3, 8), (8, 7), (7, 5), (5, 6), (5, 7)],
                                                            [(3, 8), (8, 7), (7, 4), (4, 3), (4, 2)],
                                                            [(3, 8), (8, 7), (8, 6)]],
                                                           [[], [(4, 3), (4, 2)], [(6, 8), (8, 7)], [],
                                                            [(4, 3), (4, 2)], [(4, 3), (4, 2)], [(7, 5), (7, 4)],
                                                            [(1, 2), (1, 3)], [(5, 6), (5, 7)], [(6, 1), (6, 8)],
                                                            [(5, 6), (5, 7)], [(4, 2), (4, 3)], [(8, 7), (8, 6)]],
                                                           ["TCTCAC", "TCAAAC", "TGTAAC", "TCTGAC",
                                                            "TCTCAAC", "TCATAAC", "TCCTAAC",
                                                            "TCTGAAC", "TCTTAAC", "TGCTAAC",
                                                            "TCTAC", "TCAAC", "TTAAC"],
                                                           [4, 3, 2, 4, 4, 3, 2, 4, 3, 2, 4, 3, 2]):
        points = []
        for angle in arange(0, 360, 45):
            radio = deg2rad(angle + bias)
            points.append([sin(radio) * length + center[0], cos(radio) * length + center[1]])
        points = array(points)
        draw_graph(ps=points, s=10, w=0.5, se=select, h=highlight)
        if sequence not in ["TCTCAC", "TCTGAC"]:
            for location, nucleotide in enumerate(sequence):
                if location == change - 1:
                    pyplot.fill_between([center[0] - 0.08 + location * 0.03 - 0.016,
                                         center[0] - 0.08 + location * 0.03 + 0.016],
                                        center[1] + 0.155, center[1] + 0.155 - 0.05, lw=0., color="lightblue")
                if location == len(select) - 3:
                    pyplot.text(center[0] - 0.08 + location * 0.03, center[1] + 0.125,
                                nucleotide, color="red", va="center", ha="center", fontsize=8)
                elif location <= change - 1:
                    pyplot.text(center[0] - 0.08 + location * 0.03, center[1] + 0.125,
                                nucleotide, va="center", ha="center", fontsize=8)
                else:
                    pyplot.scatter([center[0] - 0.08 + location * 0.03], [center[1] + 0.125], marker=".",
                                   s=8, color="black")
        else:
            for location, nucleotide in enumerate(sequence):
                pyplot.text(center[0] - 0.08 + location * 0.03, center[1] + 0.125,
                            nucleotide, va="center", ha="center", fontsize=8)
                if location == change - 1:
                    pyplot.fill_between([center[0] - 0.08 + location * 0.03 - 0.016,
                                         center[0] - 0.08 + location * 0.03 + 0.016],
                                        center[1] + 0.155, center[1] + 0.155 - 0.05, lw=0., color="lightblue")
    pyplot.vlines(3.6, 0, 1, lw=0.75, ls="--")
    for index, (sequence_1, sequence_2) in enumerate(zip(["TCTCAC", "TCTGAC"], ["CCT", "GCT"])):
        for location, nucleotide in enumerate(sequence_1):
            pyplot.text(3.65 + location * 0.03, 0.6 - (index * 0.4),
                        nucleotide, va="center", ha="center", fontsize=8)
        pyplot.annotate("", xy=(3.73, 0.9 - (index * 0.4)), xytext=(3.73, 0.6 - (index * 0.4)),
                        arrowprops=dict(arrowstyle="-|>", color="#0070C0", lw=2, shrinkA=9, shrinkB=9))
        pyplot.text(3.88, 0.75 - (index * 0.4), "generate\nVT-check",
                    va="center", ha="center", color="#0070C0", fontsize=9)
        for location, nucleotide in enumerate(sequence_2):
            if index == 0 and location == 0:
                pyplot.text(3.7 + location * 0.03, 0.9 - (index * 0.4),
                            nucleotide, va="center", ha="center", color="red", fontsize=8)
            else:
                pyplot.text(3.7 + location * 0.03, 0.9 - (index * 0.4),
                            nucleotide, va="center", ha="center", fontsize=8)
    pyplot.plot([3.63, 3.82, 3.82, 3.63, 3.63], [0.24, 0.24, 0.17, 0.17, 0.24], color="#0070C0", lw=0.75)
    pyplot.text(3.73, 0.1, "solution", va="center", ha="center", color="#0070C0", fontsize=9)
    pyplot.axis("off")
    pyplot.xlim(0.0, 4.0)
    pyplot.ylim(-0.05, 1.15)

    data = load_data(load_path="./raw/correction_evaluation_1.pkl")

    ax = pyplot.subplot(grid[1, 0])
    correction_rates = []
    for index, values in enumerate(data.values()):
        counts = values[0][:, 3]
        counts = counts[where(values[0][:, 4] > 0)]
        correction_rates.append(len(counts[counts == 1]) / 10000.0)
    pyplot.bar(arange(len(correction_rates))[:9], correction_rates[:9], fc="#A9D18E", ec="k", lw=0.75)
    for index in range(9):
        pyplot.text(index, correction_rates[index] + 0.01, "%.1f" % (correction_rates[index] * 100) + "%",
                    va="center", ha="center", fontsize=9)
    pyplot.xlabel("error rate", fontsize=10)
    pyplot.ylabel("correction rate", fontsize=10)
    pyplot.xticks(arange(0, 10), ["%.1f" % (value / 200.0 * 100.0) + "%" for value in arange(0, 10)], fontsize=10)
    pyplot.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0], ["50%", "60%", "70%", "80%", "90%", "100%"], fontsize=10)
    pyplot.xlim(-0.6, 8.6)
    pyplot.ylim(0.50, 1.01)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    data = load_data(load_path="./raw/correction_evaluation_2.pkl")

    used_colors = pyplot.get_cmap("RdYlBu_r")(linspace(0.2, 0.8, 8))
    ax = pyplot.subplot(grid[1, 1])
    locations = [1]
    for error_rate, values in data.items():
        if error_rate <= 0.040:
            counts, proportions, index = linspace(1, max(values), max(values)), [], int(error_rate * 200) - 1
            for count in counts:
                proportions.append(len(where(values == count)[0]))
            proportions = array(proportions, dtype=float)
            proportions /= sum(proportions)
            proportions = cumsum(proportions)
            locations.append(counts[-1])
            pyplot.plot(log10(counts), proportions, color=used_colors[index], lw=3, zorder=-index,
                        label="%.1f" % (error_rate * 100))
            pyplot.scatter([log10(max(values))], [proportions[-1]], s=20, fc=used_colors[index],
                           ec="k", lw=0.75, zorder=len(data) + 1)
            pyplot.vlines(log10(max(values)), 0, 1, ls="--", lw=0.75, zorder=0)

    locations.append(20)
    locations = array(locations)
    pyplot.legend(loc="lower right", ncol=2, fontsize=9, framealpha=1, title="% of edit errors", title_fontsize=9)
    pyplot.xlabel("reads number", fontsize=10)
    pyplot.ylabel("correction rate", fontsize=10)
    pyplot.xticks(log10(locations), locations.astype(int), fontsize=10)
    pyplot.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0], ["50%", "60%", "70%", "80%", "90%", "100%"], fontsize=10)
    pyplot.xlim(0, log10(20))
    pyplot.ylim(0.50, 1.01)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    figure.align_labels()
    figure.text(0.020, 0.99, "a", va="center", ha="center", fontsize=14)
    figure.text(0.020, 0.54, "b", va="center", ha="center", fontsize=14)
    figure.text(0.510, 0.54, "c", va="center", ha="center", fontsize=14)

    pyplot.savefig("./show/main02.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def main03():
    figure = pyplot.figure(figsize=(10, 6), tight_layout=True)

    ax = pyplot.subplot(2, 2, 1)

    matrix = load_data("./raw/correction_evaluation_3.pkl")["loss"]

    used_colors = pyplot.get_cmap("RdYlBu_r")(linspace(0.2, 0.8, 8))
    for index, values in enumerate(matrix):
        values = log10(values[values > 0])
        pyplot.scatter(arange(len(values)) + 0.5, values + 0.5, s=40, fc=used_colors[index], ec="k",
                       lw=0.75, label=("%.1f" % ((index + 1) / 2)), zorder=2)
        pyplot.plot(arange(len(values)) + 0.5, values + 0.5, lw=0.75, ls="--", color="k", zorder=1)
    pyplot.legend(loc="upper right", ncol=2, framealpha=1, title="% of edit errors", fontsize=9, title_fontsize=9)
    pyplot.xlabel("reads number", fontsize=10)
    pyplot.ylabel("maximum loss after retrieval", fontsize=10)
    pyplot.xticks(arange(10) + 0.5, (arange(10) + 1) * 5, fontsize=10)
    pyplot.yticks(arange(6) + 0.5, ["1E+" + str(v) for v in arange(6)], fontsize=10)
    pyplot.xlim(0, 10)
    pyplot.ylim(0.2, 5.8)
    # noinspection PyUnresolvedReferences
    ax.spines["top"].set_visible(False)
    # noinspection PyUnresolvedReferences
    ax.spines["right"].set_visible(False)

    ax = pyplot.subplot(2, 2, 2)

    matrix = load_data("./raw/correction_evaluation_3.pkl")["gap"].T

    show_data = matrix.copy().astype(float)
    show_data[show_data <= 0] = nan

    # noinspection PyUnresolvedReferences
    ax.add_patch(patches.FancyBboxPatch((0.25, 5.85), width=3.3, height=1.9,
                                        fc="w", ec="silver", boxstyle="Round,pad=0.1", lw=0.75))
    pyplot.text(1.9, 7.3, "smallest gap between\ncorrect-incorrect reads", va="center", ha="center", fontsize=9,
                zorder=2)

    for index, color in enumerate(pyplot.get_cmap("gist_earth_r")(linspace(0, 0.5, 40))):
        pyplot.fill_between([0.4 + index * (3.0 / 40.0), 0.4 + (index + 1) * (3.0 / 40.0)], 6.5, 6.8,
                            fc=color, lw=0, zorder=2)
    for index, number in enumerate([0, 10, 20, 30, 40]):
        pyplot.vlines(0.4 + 0.75 * index, 6.5, 6.3, lw=0.75, zorder=2)
        pyplot.text(0.4 + 0.75 * index, 6.0, str(number), va="center", ha="center", fontsize=9, zorder=2)
    pyplot.plot([0.4, 3.4, 3.4, 0.4, 0.4], [6.8, 6.8, 6.5, 6.5, 6.8], color="k", lw=0.75, zorder=3)

    pyplot.pcolormesh(arange(11), arange(9), show_data.T, vmin=0, vmax=70, cmap="gist_earth_r")
    for x in range(10):
        for y in range(8):
            if matrix[x, y] > 0:
                pyplot.text(x + 0.5, y + 0.5, matrix[x, y], va="center", ha="center", fontsize=9)

    pyplot.xlabel("reads number", fontsize=10)
    pyplot.ylabel("error rate", fontsize=10)
    pyplot.xticks(arange(10) + 0.5, (arange(10) + 1) * 5, fontsize=10)
    pyplot.yticks(arange(8) + 0.5, ["%.1f" % ((v + 1) / 2) + "%" for v in arange(8)], fontsize=10)
    pyplot.xlim(0, 10)
    pyplot.ylim(0, 8)
    # noinspection PyUnresolvedReferences
    ax.spines["top"].set_visible(False)
    # noinspection PyUnresolvedReferences
    ax.spines["right"].set_visible(False)

    ax = pyplot.subplot(2, 2, 3)

    data = load_data("./raw/correction_evaluation_5.pkl")
    used_colors = pyplot.get_cmap("RdYlBu_r")(linspace(0.2, 0.8, 8))
    used_colors = [used_colors[-1], used_colors[0]]
    used_styles = ["-", "--", ":"]

    for error_rate, color in zip([0.04, 0.005], used_colors):
        for retrieval_rate, style in zip([1.0, 0.999, 0.99], used_styles):
            x_values, y_values = data[(error_rate, retrieval_rate)]
            info_1 = "%.1f" % (error_rate * 100) + "%"
            info_2 = "%.1f" % (retrieval_rate * 100) + "%" if retrieval_rate < 1 else "lossless"
            pyplot.plot(log10(x_values), y_values, color=color, lw=2, ls=style, zorder=1,
                        label=info_1 + " | " + info_2)

    pyplot.legend(loc="upper left", ncol=2, fontsize=9, title_fontsize=9,
                  title="  minimum reads number\nfor \"error rate | retrieval rate\"")
    pyplot.xlabel("sequence diversity", fontsize=10)
    pyplot.ylabel("reads number", fontsize=10)
    pyplot.xticks(arange(0, 13),
                  ["1E+" + str(value).zfill(2) if value % 3 == 0 else "" for value in arange(0, 13)], fontsize=10)
    pyplot.yticks(arange(0, 281, 40), arange(0, 281, 40), fontsize=10)
    pyplot.xlim(0, 12)
    pyplot.ylim(0, 280)
    # noinspection PyUnresolvedReferences
    ax.spines["top"].set_visible(False)
    # noinspection PyUnresolvedReferences
    ax.spines["right"].set_visible(False)

    ax = pyplot.subplot(2, 2, 4)
    used_styles = ["-", "--", ":"]

    # create by symbolic regression technique (4.0%).
    def function_4_0(diversity, depth):
        return (log10(depth + 1) * log10(diversity)) / 0.346

    # create by symbolic regression technique (0.5%).
    def function_0_5(diversity, depth):
        return (log10(depth + 1) * log10(diversity)) / 1.601

    for selected_diversity, diversity_name, style in zip([1e12, 1e9, 1e6], ["1E+12", "1E+09", "1E+06"], used_styles):
        x_values = linspace(1, 280, 280)
        y_values = array([function_4_0(selected_diversity, selected_depth) for selected_depth in x_values])
        locations = where(x_values > y_values)
        x_values, y_values = x_values[locations], y_values[locations]
        pyplot.plot(x_values, y_values, color=used_colors[0], lw=2, ls=style, label="4.0% | " + diversity_name)

    for selected_diversity, diversity_name, style in zip([1e12, 1e9, 1e6], ["1E+12", "1E+09", "1E+06"], used_styles):
        x_values = linspace(1, 280, 280)
        y_values = array([function_0_5(selected_diversity, selected_depth) for selected_depth in x_values])
        locations = where(x_values > y_values)
        x_values, y_values = x_values[locations], y_values[locations]
        pyplot.plot(x_values, y_values, color=used_colors[1], lw=2, ls=style, label="0.5% | " + diversity_name)

    pyplot.legend(loc="upper left", ncol=2, fontsize=9, title_fontsize=9,
                  title="    non-blocking reads threshold\nfor \"error rate | sequence diversity\"")
    pyplot.xlabel("reads number", fontsize=10)
    pyplot.ylabel("reads threshold", fontsize=10)
    pyplot.xticks(arange(0, 281, 40), arange(0, 281, 40), fontsize=10)
    pyplot.yticks(arange(0, 141, 20), arange(0, 141, 20), fontsize=10)
    pyplot.xlim(0, 280)
    pyplot.ylim(0, 140)
    # noinspection PyUnresolvedReferences
    ax.spines["top"].set_visible(False)
    # noinspection PyUnresolvedReferences
    ax.spines["right"].set_visible(False)

    figure.align_labels()
    figure.text(0.020, 0.99, "a", va="center", ha="center", fontsize=14)
    figure.text(0.518, 0.99, "b", va="center", ha="center", fontsize=14)
    figure.text(0.020, 0.49, "c", va="center", ha="center", fontsize=14)
    figure.text(0.518, 0.49, "d", va="center", ha="center", fontsize=14)

    pyplot.savefig("./show/main03.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


def main04():
    figure = pyplot.figure(figsize=(10, 6), tight_layout=True)

    data = load_data(load_path="./raw/correction_evaluation_1.pkl")
    correction_speeds = []
    for key, values in data.items():
        if key <= 0.04:
            correction_speeds.append(200.0 / (values[1] / 10000.0))
    correction_speeds = array(correction_speeds)

    ax = pyplot.subplot(2, 1, 1)

    # outsize for legend.
    pyplot.bar([-1], [-1], fc="#E2F0D9", ec="k", lw=0.75, width=0.3, label="high-throughput sequencing")

    for index, speed in enumerate(correction_speeds):
        height = (speed - 100000) / 250 + 1600
        pyplot.vlines(index - 0.3, 0, 1300, lw=0.75, zorder=0)
        pyplot.vlines(index, 0, 1300, lw=0.75, zorder=0)
        pyplot.plot([index - 0.3, index - 0.3, index, index], [1500, height, height, 1500],
                    lw=0.75, color="k", zorder=0)
        pyplot.text(index - 0.15, height + 20, int(speed), va="bottom", ha="center", fontsize=9)
        pyplot.fill_between([index - 0.3, index], 0, 1300, fc="#E2F0D9", ec="w", lw=0, zorder=-1)
        pyplot.fill_between([index - 0.3, index], 1500, height, fc="#E2F0D9", ec="w", lw=0, zorder=-1)
    pyplot.text(-0.5, 1400, " ", bbox=dict(fc="w", ec="w"), va="center", ha="center", fontsize=1, zorder=5)
    figure.add_artist(lines.Line2D([0.100, 0.108], [0.703, 0.713], lw=1, color="k"))
    figure.add_artist(lines.Line2D([0.100, 0.108], [0.728, 0.738], lw=1, color="k"))
    pyplot.bar(arange(len(correction_speeds)) + 0.15, correction_speeds / 450, fc="#A9D18E", ec="k", lw=0.75, width=0.3,
               label="single-molecule sequencing")
    for h_location, value in enumerate(correction_speeds / 450):
        pyplot.text(h_location + 0.15, value + 20, int(value), va="bottom", ha="center", fontsize=9)
    pyplot.legend(loc="upper right", fontsize=9)
    pyplot.xlabel("error rate", fontsize=10)
    pyplot.ylabel("performance\n(correction speed / sequencing speed)", fontsize=10)
    pyplot.xticks(arange(len(correction_speeds)),
                  ["%.1f" % (value / 2) + "%" for value in arange(len(correction_speeds))], fontsize=10)
    pyplot.yticks([0, 400, 800, 1200, 1600, 2000, 2400, 2800, 3200, 3600],
                  [0, 400, 800, 1200, 100000, 200000, 300000, 400000, 500000, 600000], fontsize=10)
    pyplot.xlim(-0.5, len(correction_speeds) - 0.5)
    pyplot.ylim(0, 3600)
    # noinspection PyUnresolvedReferences
    ax.spines["top"].set_visible(False)
    # noinspection PyUnresolvedReferences
    ax.spines["right"].set_visible(False)

    def previous(sequence_diversity):
        return 200 * (sequence_diversity ** 2) + 48000200 * sequence_diversity

    def propose40(sequence_diversity):
        a = 7267.904 * sequence_diversity * (log10(sequence_diversity) ** 2)
        b = 76672.229 * sequence_diversity
        return a + b

    def propose05(sequence_diversity):
        a = 183.178 * sequence_diversity * (log10(sequence_diversity) ** 2)
        b = 50508.157 * sequence_diversity
        return a + b

    values_1, values_2, values_3 = [], [], []
    for n_i, number in enumerate(10 ** linspace(0, 20, 21)):
        value_1 = previous(number)
        value_2 = propose40(number)
        value_3 = propose05(number)
        values_1.append(value_1)
        values_2.append(value_2)
        values_3.append(value_3)

    ax = pyplot.subplot(2, 1, 2)

    info = 200 * 0.581

    locations = [0]
    for scale in range(6):
        location = log10(8 * 1024 ** (scale + 1) / info)
        locations.append(location)
        pyplot.vlines(location, 0, 80, lw=0.75, ls="--", zorder=1)

    locations.append(20)

    for index, unit in enumerate(["~ KB", "KB ~ MB", "MB ~ GB", "GB ~ TB",
                                  "TB ~ PB", "PB ~ EB", "EB ~"]):
        if unit != "EB ~":
            pyplot.annotate("", xy=(locations[index], 48), xytext=(locations[index + 1], 48),
                            arrowprops=dict(arrowstyle="<|-|>", color="k", lw=0.75), zorder=1)
        else:
            pyplot.annotate("", xy=(locations[index], 48), xytext=(locations[index + 1], 48),
                            arrowprops=dict(arrowstyle="-|>", color="k", lw=0.75), zorder=1)

        pyplot.text((locations[index] + locations[index + 1]) / 2, 47.5, unit, zorder=2,
                    bbox=dict(fc="w", ec="w"), va="center", ha="center", fontsize=9)

    pyplot.plot(linspace(0, 20, 21), log10(values_1), color="#999999", lw=2,
                label="conventional pipeline", zorder=2)
    pyplot.plot(linspace(0, 20, 21), log10(values_2), color="#4D7731", lw=2,
                label="our pipeline (4.0% of edit errors)", zorder=2)
    pyplot.plot(linspace(0, 20, 21), log10(values_3), color="#4D7731", lw=2, ls="--",
                label="our pipeline (0.5% of edit errors)", zorder=2)
    pyplot.legend(loc="lower right", framealpha=1, fontsize=9)
    pyplot.xlabel("sequence diversity", fontsize=10)
    pyplot.ylabel("approximated operations", fontsize=10, labelpad=8)
    pyplot.xticks(linspace(0, 20, 11), ["1E+" + str(int(i)).zfill(2) for i in linspace(0, 20, 11)], fontsize=10)
    pyplot.yticks(linspace(0, 50, 6), ["1E+" + str(int(i)).zfill(2) for i in linspace(0, 50, 6)], fontsize=10)
    pyplot.xlim(0, 20)
    pyplot.ylim(0, 50)
    # noinspection PyUnresolvedReferences
    ax.spines["top"].set_visible(False)
    # noinspection PyUnresolvedReferences
    ax.spines["right"].set_visible(False)

    figure.align_labels()
    figure.text(0.03, 0.99, "a", va="center", ha="center", fontsize=14)
    figure.text(0.03, 0.50, "b", va="center", ha="center", fontsize=14)

    pyplot.savefig("./show/main04.pdf", format="pdf", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    main01()
    # main02()
    # main03()
    # main04()
