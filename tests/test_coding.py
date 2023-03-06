from numpy import array, all
from unittest import TestCase

from dsw import encode, decode


class TestNormalEncode(TestCase):

    def setUp(self):
        self.accessor = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])
        self.start_index = 1
        self.vt_length = 4
        self.binary_message = array([0, 1, 0, 1, 0, 1, 0, 1])

    def test(self):
        dna_sequence = encode(accessor=self.accessor, binary_message=self.binary_message, start_index=self.start_index)
        self.assertEqual("TCTCTCT", dna_sequence)
        dna_sequence, vt_list = encode(accessor=self.accessor, binary_message=self.binary_message,
                                       start_index=self.start_index, vt_length=self.vt_length)
        self.assertEqual("TCTCTCT", dna_sequence)
        self.assertEqual("TAGC", vt_list)


class TestNormalDecode(TestCase):

    def setUp(self):
        self.accessor = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])
        self.start_index = 1
        self.vt_length = 4
        self.dna_sequence = "TCTCTCT"

    def test(self):
        binary_message = decode(accessor=self.accessor, dna_sequence=self.dna_sequence, start_index=1, bit_length=8)
        self.assertEqual(all(binary_message == array([0, 1, 0, 1, 0, 1, 0, 1])), True)


class TestNormalMap(TestCase):

    def setUp(self):
        self.accessor = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])
        self.seed = 2021


class TestFasterEncode(TestCase):

    def setUp(self):
        self.accessor = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])
        self.start_index = 1
        self.vt_length = 4
        self.binary_message = array([0, 1, 0, 1, 0, 1, 0, 1])

    def test(self):
        dna_sequence = encode(accessor=self.accessor, binary_message=self.binary_message, start_index=self.start_index,
                              is_faster=True)
        self.assertEqual("AGAGAGAG", dna_sequence)
        dna_sequence, vt_list = encode(accessor=self.accessor, binary_message=self.binary_message,
                                       start_index=self.start_index, vt_length=self.vt_length, is_faster=True)
        self.assertEqual("AGAGAGAG", dna_sequence)
        self.assertEqual("AATA", vt_list)


class TestFasterDecode(TestCase):

    def setUp(self):
        self.accessor = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])
        self.start_index = 1
        self.vt_length = 4
        self.dna_sequence = "AGAGAGAG"

    def test(self):
        binary_message = decode(accessor=self.accessor, dna_sequence=self.dna_sequence, start_index=1, bit_length=8,
                                is_faster=True)
        self.assertEqual(all(binary_message == array([0, 1, 0, 1, 0, 1, 0, 1])), True)
