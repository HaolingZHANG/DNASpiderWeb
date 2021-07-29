__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import random
from unittest import TestCase
from dsw.generate import obtain_number, obtain_oligo


class TestNumber(TestCase):

    def setUp(self):
        random.seed(2021)
        self.oligo_length = 10
        self.test_oligo_groups = [["ACGT"[random.randint(0, 3)] for _ in range(self.oligo_length)] for _ in range(10)]
        self.verify_numbers = [937770, 807052, 502551, 991128, 88040, 531838, 908835, 656257, 889866, 144398]

    def test(self):
        results = [obtain_number(oligo=oligo) for oligo in self.test_oligo_groups]
        self.assertEqual(results, self.verify_numbers)


class TestOligo(TestCase):

    def setUp(self):
        random.seed(2021)
        self.oligo_length = 10
        self.test_numbers = [random.randint(0, 4 ** self.oligo_length) for _ in range(10)]
        self.verify_oligo_groups = ["TATGTTCAAA", "GATCGGTCAC", "CTTGGGAGCC", "ACACGTTTAG", "TGATAATTGG",
                                    "TTAGGCCAGT", "AGAAGGAGGG", "GGAAATAAAA", "GAGGATCGAG", "GCCTGCGACG"]

    def test(self):
        results = [obtain_oligo(decimal_number=number, length=self.oligo_length) for number in self.test_numbers]
        self.assertEqual(results, self.verify_oligo_groups)


class TestTransform(TestCase):

    def setUp(self):
        self.oligo_length = 10

    def test(self):
        for number in range(4 ** self.oligo_length):
            oligo = obtain_oligo(decimal_number=number, length=self.oligo_length)
            transformed_number = obtain_number(oligo=oligo)
            self.assertEqual(number, transformed_number)