__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import numpy
from unittest import TestCase
from dsw.algorithm import encode, decode


class TestTranscode(TestCase):

    def setUp(self):
        numpy.random.seed(2021)
        vertices = numpy.load(file="../outputs/screen1[vertices].npy")
        self.length = 20
        self.number = 20
        self.start_index = numpy.random.choice(numpy.where(vertices == 1)[0])
        self.graph = numpy.load(file="../outputs/screen1[graph].npy")
        self.matrix = numpy.random.randint(low=0, high=2, size=(self.number, self.length))

    def test(self):
        for bits in self.matrix:
            oligo = encode(bits=bits, graph=self.graph, start_index=self.start_index)
            recovered_bits = decode(oligo=oligo, bit_length=self.length, graph=self.graph, start_index=self.start_index)
            self.assertEqual(bits.tolist(), recovered_bits)
