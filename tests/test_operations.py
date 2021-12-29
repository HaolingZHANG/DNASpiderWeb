__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from unittest import TestCase

from dsw import calculus_addition, calculus_subtraction, calculus_multiplication, calculus_division


class TestAddition(TestCase):

    def setUp(self):
        self.numbers = list(range(1000))
        self.bases = [0, 1, 2, 3, 4]

    def test(self):
        for number in self.numbers:
            for base in self.bases:
                requested = str(number + base)
                predicted = calculus_addition(number=str(number), base=str(base))
                self.assertEqual(requested, predicted)


class TestSubtraction(TestCase):

    def setUp(self):
        self.numbers = list(range(4, 1000))
        self.bases = [0, 1, 2, 3, 4]

    def test(self):
        for number in self.numbers:
            for base in self.bases:
                requested = str(number - base)
                predicted = calculus_subtraction(number=str(number), base=str(base))
                self.assertEqual(requested, predicted)


class TestMultiplication(TestCase):

    def setUp(self):
        self.numbers = list(range(1000))
        self.bases = [1, 2, 3, 4]

    def test(self):
        for number in self.numbers:
            for base in self.bases:
                requested = str(number * base)
                predicted = calculus_multiplication(number=str(number), base=str(base))
                self.assertEqual(requested, predicted)


class TestDivision(TestCase):

    def setUp(self):
        self.numbers = list(range(4, 1000))
        self.bases = [2, 3, 4]

    def test(self):
        for number in self.numbers:
            for base in self.bases:
                requested = (str(number // base), str(number % base))
                predicted = calculus_division(number=str(number), base=str(base))
                self.assertEqual(requested, predicted)
