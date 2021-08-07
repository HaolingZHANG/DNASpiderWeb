__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import numpy
from unittest import TestCase
from dsw.generator import to_graph
from dsw.evaluator import approximate_upper_bound


class TestTerminals(TestCase):

    def setUp(self):
        self.unconstrained_graph = numpy.array([[1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],    # AA
                                                [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],    # AC
                                                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],    # AG
                                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],    # AT
                                                [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],    # CA
                                                [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],    # CC
                                                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],    # CG
                                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],    # CT
                                                [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],    # GA
                                                [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],    # GC
                                                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],    # GG
                                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],    # GT
                                                [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],    # TA
                                                [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],    # TC
                                                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],    # TG
                                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]])   # TT

        self.empty_graph = numpy.zeros(shape=(16, 16), dtype=numpy.int)

    def test(self):
        self.assertEqual(approximate_upper_bound(to_graph(self.unconstrained_graph)), 1.0)
        self.assertEqual(approximate_upper_bound(to_graph(self.empty_graph)), 0.0)


class TestArtificialExamples(TestCase):

    def setUp(self):
        # TODO
        self.artificial_graphs = [(numpy.array([1]), 1),
                                  (numpy.array([1]), 1),
                                  (numpy.array([1]), 1),
                                  (numpy.array([1]), 1),
                                  (numpy.array([1]), 1)]

    def test(self):
        pass
