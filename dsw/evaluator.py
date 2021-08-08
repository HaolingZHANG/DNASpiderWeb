__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import math
import numpy
from dsw import n_system


def approximate_upper_bound(graph):
    """
    Approximate capacity of the specific graph through Perron–Frobenius theorem: lambda_max(graph).
    |cite| Oskar Perron (1907) Mathematische Annalen
    |cite| Georg Frobenius (1912) Sitzungsberichte der Königlich Preussischen Akademie der Wissenschaften
    |cite| Brian H. Marcus et al. (2001) Lecture notes

    :param graph: graph of DNA Spider-Web.
    :type graph: numpy.ndarray

    :return: capacity of this graph.
    :rtype: float
    """
    # calculate bound for an empty graph.
    if numpy.all(graph == -1):
        return 0.0

    # approximate dominant (or maximum) eigenvalue by Von Mises iteration (also known as the power method).
    # |cite| Richard von Mises et al. (1929) Zeitschrift für Angewandte Mathematik und Mechanik
    # |cite| Qingkai Kong et al. (2020) Academic Press
    last_eigenvector, last_eigenvalue = numpy.ones(shape=(len(graph),), dtype=numpy.float), -math.inf
    while True:
        eigenvector = numpy.zeros_like(last_eigenvector)
        for vertex_index, vertex in enumerate(graph):
            eigenvector[vertex_index] = numpy.sum(last_eigenvector[vertex[vertex >= 0]])
        eigenvalue = numpy.max(eigenvector)
        eigenvector = eigenvector / eigenvalue
        if last_eigenvalue is not None:
            if abs(eigenvalue - last_eigenvalue) < 1e-12:
                # logarithm of the maximum output degree (if degenerate nucleotides are not considered, it is 4).
                return math.log(eigenvalue, len(n_system))

        last_eigenvalue, last_eigenvector = eigenvalue, eigenvector
