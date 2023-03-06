from numpy import array
from unittest import TestCase

from dsw import set_vt, repair_dna


class TestRepair(TestCase):

    def setUp(self):
        self.accessor = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                               [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])
        self.dna_sequence = "TCTCTATCTCTC"  # "TCTCTCTCTCTC" is original DNA sequence
        self.vt_check = set_vt(dna_sequence="TCTCTCTCTCTC", vt_length=4)  # check list of Varshamov-Tenengolts code.

    def test(self):
        repaired_dna_sequences, additions = repair_dna(dna_sequence=self.dna_sequence,
                                                       accessor=self.accessor, start_index=1,
                                                       observed_length=2, has_indel=True)
        self.assertEqual(repaired_dna_sequences, ["TCTCTCTCTCTC", "TCTCTGTCTCTC"])
        self.assertEqual(additions, (1, False, 2, 14))

        repaired_dna_sequences, additions = repair_dna(dna_sequence=self.dna_sequence, vt_check=self.vt_check,
                                                       accessor=self.accessor, start_index=1,
                                                       observed_length=2, has_indel=True)
        self.assertEqual(repaired_dna_sequences, ["TCTCTCTCTCTC"])
        self.assertEqual(additions, (1, True, 2, 14))
