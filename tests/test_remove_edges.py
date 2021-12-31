from numpy import random, array, all
from unittest import TestCase

from dsw import remove_random_edges, get_complete_accessor, encode, decode


class TestRemoveEdges(TestCase):

    def setUp(self):
        self.random_seed = 2021
        self.repeats = 20
        self.bit_length = 100
        self.accessor = get_complete_accessor(observed_length=5)
        self.remove_range = [5, 20]
        self.edge_group = []
        for remove_number in range(self.remove_range[0], self.remove_range[1] + 1):
            self.edge_group.append(remove_random_edges(accessor=self.accessor, remove_number=remove_number,
                                                       random_seed=self.random_seed))
        self.start_indices = array(list(range(4 ** 5)))
        random.shuffle(self.start_indices)
        self.start_indices = self.start_indices[: self.repeats]
        self.bit_matrix = random.randint(low=0, high=2, size=(self.repeats, self.bit_length), dtype=int)

    def test(self):
        for source in self.bit_matrix:
            for removed_edges in self.edge_group:
                for start_index in self.start_indices:
                    oligo = encode(binary_message=source, accessor=self.accessor,
                                   removed_edges=removed_edges, start_index=start_index)
                    target = decode(dna_string=oligo, accessor=self.accessor, bit_length=self.bit_length,
                                    removed_edges=removed_edges, start_index=start_index)
                    self.assertEqual(all(source == target), True)
