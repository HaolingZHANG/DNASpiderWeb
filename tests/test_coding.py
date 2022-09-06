from numpy import array, all
from unittest import TestCase

from dsw import encode, decode, set_vt, repair_dna


class TestEncode(TestCase):

    def setUp(self):
        self.accessor = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])
        self.start_index = 1
        self.vt_length = 4
        self.binary_message = array([0, 1, 0, 1, 0, 1, 0, 1])

    def test(self):
        dna_string = encode(accessor=self.accessor, binary_message=self.binary_message, start_index=self.start_index)
        self.assertEqual("TCTCTCT", dna_string)
        dna_string, vt_list = encode(accessor=self.accessor, binary_message=self.binary_message,
                                     start_index=self.start_index, vt_length=self.vt_length)
        self.assertEqual("TCTCTCT", dna_string)
        self.assertEqual("TAGC", vt_list)


class TestDecode(TestCase):

    def setUp(self):
        self.accessor = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])
        self.start_index = 1
        self.vt_length = 4
        self.dna_string = "TCTCTCT"

    def test(self):
        binary_message = decode(accessor=self.accessor, dna_string=self.dna_string, start_index=1, bit_length=8)
        self.assertEqual(all(binary_message == array([0, 1, 0, 1, 0, 1, 0, 1])), True)


class TestRepair(TestCase):

    def setUp(self):
        self.accessor = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])
        self.dna_string = "TCTCTATCTCTC"  # "TCTCTCTCTCTC" is original DNA string
        self.vt_check = set_vt(dna_string="TCTCTCTCTCTC", vt_length=4)  # check list of Varshamov-Tenengolts code.

    def test(self):
        repaired_dna_strings, additions = repair_dna(dna_string=self.dna_string,
                                                     accessor=self.accessor, start_index=1,
                                                     observed_length=2, has_indel=True)
        self.assertEqual(repaired_dna_strings, ["TCTCTCTCTCTC", "TCTCTGTCTCTC"])
        self.assertEqual(additions, (1, False, 2, 15))

        repaired_dna_strings, additions = repair_dna(dna_string=self.dna_string, vt_check=self.vt_check,
                                                     accessor=self.accessor, start_index=1,
                                                     observed_length=2, has_indel=True)
        print(repaired_dna_strings, additions)
        self.assertEqual(repaired_dna_strings, ["TCTCTCTCTCTC"])
        self.assertEqual(additions, (1, True, 2, 15))
