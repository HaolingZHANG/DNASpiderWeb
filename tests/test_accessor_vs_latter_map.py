from numpy import all
from unittest import TestCase

from dsw import latter_map_to_accessor, accessor_to_latter_map, get_complete_accessor


class TestConvert(TestCase):

    def setUp(self):
        self.sizes = list(range(2, 10))
        self.checks = []
        for size in self.sizes:
            self.checks.append((size, get_complete_accessor(observed_length=size)))

    def test(self):
        for observed_length, accessor in self.checks:
            latter_map = accessor_to_latter_map(accessor=accessor)
            recovered_accessor = latter_map_to_accessor(latter_map=latter_map, observed_length=observed_length)
            self.assertEqual(all(accessor == recovered_accessor), True)
