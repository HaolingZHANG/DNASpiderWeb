__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import math
import numpy
from dsw import n_system


def approximate_upper_bound(graph, replay):
    """
    Approximate capacity of the specific graph through Perron–Frobenius theorem: lambda_max(graph).
    |cite| Oskar Perron (1907) Mathematische Annalen
    |cite| Georg Frobenius (1912) Sitzungsberichte der Königlich Preussischen Akademie der Wissenschaften
    |cite| Brian H. Marcus et al. (2001) Lecture notes

    :param graph: graph of DNA Spider-Web.
    :type graph: numpy.ndarray

    :param replay: replay time for approximating the upper bound.
    :type replay: int

    :return: capacity of this graph.
    :rtype: float
    """
    if numpy.all(graph == -1):
        return 0.0

    results = []
    for _ in range(replay):
        last_eigenvector, last_eigenvalue = numpy.random.random(size=(len(graph),)), 0.0
        while True:
            eigenvector = numpy.zeros_like(last_eigenvector)
            for positions in graph.T:
                eigenvector[numpy.where(positions >= 0)] += last_eigenvector[positions[positions >= 0]]
            eigenvalue = numpy.max(eigenvector)
            eigenvector = eigenvector / eigenvalue
            if abs(eigenvalue - last_eigenvalue) < 1e-10:
                results.append(math.log(eigenvalue, len(n_system)))
                break

            last_eigenvalue, last_eigenvector = eigenvalue, eigenvector

    return numpy.mean(results)
