from numpy import random, array, all
from unittest import TestCase

from dsw import create_random_shuffles, get_complete_accessor, encode, decode


class TestShuffles(TestCase):

    def setUp(self):
        random.seed(2021)
        self.repeats = 20
        self.bit_length = 100
        self.shuffle_count = 10
        self.shuffle_group = []
        for _ in range(self.shuffle_count):
            self.shuffle_group.append(create_random_shuffles(observed_length=5))
        self.accessor = get_complete_accessor(observed_length=5)
        self.start_indices = array(list(range(4 ** 5)))
        random.shuffle(self.start_indices)
        self.start_indices = self.start_indices[: self.repeats]
        self.bit_matrix = random.randint(low=0, high=2, size=(self.repeats, self.bit_length), dtype=int)

    def test(self):
        for source in self.bit_matrix:
            for shuffles in self.shuffle_group:
                for start_index in self.start_indices:
                    oligo = encode(binary_message=source, accessor=self.accessor,
                                   shuffles=shuffles, start_index=start_index)
                    target = decode(dna_sequence=oligo, accessor=self.accessor, bit_length=self.bit_length,
                                    shuffles=shuffles, start_index=start_index)
                    self.assertEqual(all(source == target), True)
