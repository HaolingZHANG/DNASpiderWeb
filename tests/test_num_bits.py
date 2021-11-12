__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from random import seed, randint
from unittest import TestCase

from dsw.coder import to_number, to_bits


class TestNumber(TestCase):

    def setUp(self):
        seed(2021)
        self.test_bit_groups = [[randint(0, 1) for _ in range(12)] for _ in range(10)]
        self.verify_numbers = [3294, 167, 3187, 2567, 2599, 2775, 394, 212, 887, 1858]

    def test(self):
        results = [to_number(bits) for bits in self.test_bit_groups]
        self.assertEqual(results, self.verify_numbers)


class TestBits(TestCase):

    def setUp(self):
        seed(2021)
        self.bit_length = 12
        self.test_numbers = [randint(0, 2 ** 12) for _ in range(10)]
        self.verify_bit_groups = [[1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1],
                                  [1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0],
                                  [0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0],
                                  [0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1],
                                  [1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0],
                                  [1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0],
                                  [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1],
                                  [1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1],
                                  [1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1]]

    def test(self):
        results = [to_bits(number, self.bit_length) for number in self.test_numbers]
        self.assertEqual(results, self.verify_bit_groups)


class TestTransform(TestCase):

    def setUp(self):
        self.bit_length = 10

    def test(self):
        for number in range(2 ** (self.bit_length * 2)):
            bits = to_bits(number, self.bit_length)
            transformed_number = to_number(bits)
            self.assertEqual(number, transformed_number)
