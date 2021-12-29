__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from random import seed, randint
from unittest import TestCase

from dsw import bit_to_number, number_to_bit


class TestNumber(TestCase):

    def setUp(self):
        seed(2021)
        self.test_binary_messages = [[randint(0, 1) for _ in range(12)] for _ in range(10)]
        self.verify_numbers = [3294, 167, 3187, 2567, 2599, 2775, 394, 212, 887, 1858]

    def test(self):
        results = [int(bit_to_number(bits)) for bits in self.test_binary_messages]
        self.assertEqual(results, self.verify_numbers)
        results = [bit_to_number(bits, is_string=False) for bits in self.test_binary_messages]
        self.assertEqual(results, self.verify_numbers)


class TestBits(TestCase):

    def setUp(self):
        seed(2021)
        self.bit_length = 12
        self.test_numbers = [randint(0, 2 ** 12) for _ in range(10)]
        self.verify_binary_messages = [[1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1],
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
        results = [number_to_bit(decimal_number=str(number), bit_length=self.bit_length)
                   for number in self.test_numbers]
        self.assertEqual(results, self.verify_binary_messages)
        results = [number_to_bit(decimal_number=number, bit_length=self.bit_length)
                   for number in self.test_numbers]
        self.assertEqual(results, self.verify_binary_messages)


class TestTransform(TestCase):

    def setUp(self):
        self.bit_length = 8

    def test(self):
        for requested in range(2 ** (self.bit_length * 2)):
            bits = number_to_bit(decimal_number=str(requested), bit_length=self.bit_length)
            predicted_1 = int(bit_to_number(bit_array=bits))
            predicted_2 = bit_to_number(bit_array=bits, is_string=False)
            self.assertEqual(requested, predicted_1)
            self.assertEqual(requested, predicted_2)
