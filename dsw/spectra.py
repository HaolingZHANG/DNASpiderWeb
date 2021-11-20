__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from numpy import random, zeros_like, ones, where, max, all, median, abs, log2


# noinspection PyArgumentList
def calculate_capacity(graph, replay=1, need_process=False):
    """
    Calculate capacity of the specific graph through Perron–Frobenius theorem: lambda_max(graph).
    |cite| Oskar Perron (1907) Mathematische Annalen
    |cite| Georg Frobenius (1912) Sitzungsberichte der Königlich Preussischen Akademie der Wissenschaften
    |cite| Brian H. Marcus et al. (2001) Lecture notes

    :param graph: graph of DNA Spider-Web.
    :type graph: numpy.ndarray

    :param replay: replay time for approximating the upper bound.
    :type replay: int

    :param need_process: need eigenvalue in the process.
    :type need_process: bool

    :return: capacity of this graph (and eigenvalues in the process).
    :rtype: float
    """
    if all(graph == -1):
        if need_process:
            return (0.0, [0.0]) if replay == 1 else (0.0, [[0.0] for _ in range(replay)])
        else:
            return 0.0

    results, process = [], []
    for _ in range(replay):
        process.append([])
        if replay > 1:
            last_eigenvector, last_eigenvalue = abs(random.random(size=(len(graph),))), 0.0
        else:
            last_eigenvector, last_eigenvalue = ones(shape=(len(graph),)), 0.0

        while True:
            eigenvector = zeros_like(last_eigenvector)
            for positions in graph.T:
                eigenvector[where(positions >= 0)] += last_eigenvector[positions[positions >= 0]]
            eigenvalue = max(eigenvector)
            eigenvector = eigenvector / eigenvalue
            process[-1].append(log2(eigenvalue))
            if abs(eigenvalue - last_eigenvalue) < 1e-10:
                results.append(log2(eigenvalue))
                break

            last_eigenvalue, last_eigenvector = eigenvalue, eigenvector

    if need_process:
        return (median(results), process[0]) if replay == 1 else (median(results), process)
    else:
        return median(results)
