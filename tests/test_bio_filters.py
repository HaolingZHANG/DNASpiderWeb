__author__ = "Zhang, Haoling [hlzchn@gmail.com]"

from numpy import random
from unittest import TestCase

from dsw import LocalBioFilter


class TestLocalBioFilterTotally(TestCase):

    def setUp(self):
        self.random_seed = 2021
        self.single_repeats = 1000
        self.dna_length = 80
        self.observed_length = 10
        self.reference = "ACGT" * (self.dna_length // 4)
        self.bio_filters = [LocalBioFilter(observed_length=self.observed_length, max_homopolymer_runs=2),
                            LocalBioFilter(observed_length=self.observed_length, gc_range=[0.5, 0.5]),
                            LocalBioFilter(observed_length=self.observed_length, undesired_motifs=["GCC"])]
        self.dna_string_group = [[], [], []]

        random.seed(self.random_seed)

        # for normal DNA strings.
        for _ in range(self.single_repeats):
            nucleotides = random.choice(["AAA", "CCC", "GGG", "TTT"])
            location = random.randint(0, self.dna_length - 2)
            dna_string = list(self.reference)
            dna_string[int(location)], dna_string[int(location) + 1], dna_string[int(location) + 2] = nucleotides
            self.dna_string_group[0].append((False, "".join(dna_string)))

            location = random.randint(0, self.dna_length)
            mapping = {"A": "C", "C": "A", "G": "T", "T": "G"}
            dna_string = list(self.reference)
            dna_string[location] = mapping[dna_string[location]]
            self.dna_string_group[1].append((False, "".join(dna_string)))

            location = random.randint(0, self.dna_length - 2)
            dna_string = list(self.reference)
            dna_string[int(location)], dna_string[int(location) + 1], dna_string[int(location) + 2] = ["G", "C", "C"]
            self.dna_string_group[2].append((False, "".join(dna_string)))

    def test(self):
        for dna_strings, bio_filter in zip(self.dna_string_group, self.bio_filters):
            for flag, dna_string in dna_strings:
                self.assertEqual(flag, bio_filter.valid(dna_string=dna_string, only_last=False))


class TestLocalBioFilterOnlyLast(TestCase):

    def setUp(self):
        self.random_seed = 2021
        self.single_repeats = 1000
        self.dna_length = 80
        self.observed_length = 10
        self.reference = "ACGT" * (self.dna_length // 4)
        self.bio_filters = [LocalBioFilter(observed_length=self.observed_length, max_homopolymer_runs=2),
                            LocalBioFilter(observed_length=self.observed_length, gc_range=[0.5, 0.5]),
                            LocalBioFilter(observed_length=self.observed_length, undesired_motifs=["GCC"])]
        self.dna_string_group = [[], [], []]

        random.seed(self.random_seed)

        # for normal DNA strings.
        for _ in range(self.single_repeats):
            nucleotides = random.choice(["AAA", "CCC", "GGG", "TTT"])
            location = random.randint(0, self.dna_length - 2)
            dna_string = list(self.reference)
            dna_string[int(location)], dna_string[int(location) + 1], dna_string[int(location) + 2] = nucleotides
            dna_string = "".join(dna_string)

            flag = True
            for pattern in ["AAA", "CCC", "GGG", "TTT"]:
                if pattern in dna_string[-self.observed_length:]:
                    flag = False
            self.dna_string_group[0].append((flag, dna_string))

            location = random.randint(0, self.dna_length)
            mapping = {"A": "C", "C": "A", "G": "T", "T": "G"}
            dna_string = list(self.reference)
            dna_string[location] = mapping[dna_string[location]]
            dna_string = "".join(dna_string)

            gc_count = dna_string[-self.observed_length:].count("C") + dna_string[-self.observed_length:].count("G")
            flag = gc_count == (self.observed_length // 2)
            self.dna_string_group[1].append((flag, dna_string))

            location = random.randint(0, self.dna_length - 2)
            dna_string = list(self.reference)
            dna_string[int(location)], dna_string[int(location) + 1], dna_string[int(location) + 2] = ["G", "C", "C"]
            dna_string = "".join(dna_string)

            flag = "GCC" not in dna_string[-self.observed_length:]
            self.dna_string_group[2].append((flag, dna_string))

    def test(self):
        for dna_strings, bio_filter in zip(self.dna_string_group, self.bio_filters):
            for flag, dna_string in dna_strings:
                self.assertEqual(flag, bio_filter.valid(dna_string=dna_string, only_last=True))
