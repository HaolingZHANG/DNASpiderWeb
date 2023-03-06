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
        self.dna_sequence_group = [[], [], []]

        random.seed(self.random_seed)

        # for normal DNA sequences.
        for _ in range(self.single_repeats):
            nucleotides = random.choice(["AAA", "CCC", "GGG", "TTT"])
            location = random.randint(0, self.dna_length - 2)
            dna_sequence = list(self.reference)
            dna_sequence[int(location)], dna_sequence[int(location) + 1], dna_sequence[int(location) + 2] = nucleotides
            self.dna_sequence_group[0].append((False, "".join(dna_sequence)))

            location = random.randint(0, self.dna_length)
            mapping = {"A": "C", "C": "A", "G": "T", "T": "G"}
            dna_sequence = list(self.reference)
            dna_sequence[location] = mapping[dna_sequence[location]]
            self.dna_sequence_group[1].append((False, "".join(dna_sequence)))

            location = random.randint(0, self.dna_length - 2)
            dna_sequence = list(self.reference)
            dna_sequence[int(location)], dna_sequence[int(location) + 1], dna_sequence[int(location) + 2] = ["G", "C", "C"]
            self.dna_sequence_group[2].append((False, "".join(dna_sequence)))

    def test(self):
        for dna_sequences, bio_filter in zip(self.dna_sequence_group, self.bio_filters):
            for flag, dna_sequence in dna_sequences:
                self.assertEqual(flag, bio_filter.valid(dna_sequence=dna_sequence, only_last=False))


class TestLocalBioFilterOnlyLast(TestCase):

    # noinspection PyUnresolvedReferences
    def setUp(self):
        self.random_seed = 2021
        self.single_repeats = 1000
        self.dna_length = 80
        self.observed_length = 10
        self.reference = "ACGT" * (self.dna_length // 4)
        self.bio_filters = [LocalBioFilter(observed_length=self.observed_length, max_homopolymer_runs=2),
                            LocalBioFilter(observed_length=self.observed_length, gc_range=[0.5, 0.5]),
                            LocalBioFilter(observed_length=self.observed_length, undesired_motifs=["GCC"])]
        self.dna_sequence_group = [[], [], []]

        random.seed(self.random_seed)

        # for normal DNA sequences.
        for _ in range(self.single_repeats):
            nucleotides = random.choice(["AAA", "CCC", "GGG", "TTT"])
            location = random.randint(0, self.dna_length - 2)
            dna_sequence = list(self.reference)
            dna_sequence[int(location)], dna_sequence[int(location) + 1], dna_sequence[int(location) + 2] = nucleotides
            dna_sequence = "".join(dna_sequence)

            flag = True
            for pattern in ["AAA", "CCC", "GGG", "TTT"]:
                if pattern in dna_sequence[-self.observed_length:]:
                    flag = False
            self.dna_sequence_group[0].append((flag, dna_sequence))

            location = random.randint(0, self.dna_length)
            mapping = {"A": "C", "C": "A", "G": "T", "T": "G"}
            dna_sequence = list(self.reference)
            dna_sequence[location] = mapping[dna_sequence[location]]
            dna_sequence = "".join(dna_sequence)

            gc_count = dna_sequence[-self.observed_length:].count("C") + dna_sequence[-self.observed_length:].count("G")
            flag = gc_count == (self.observed_length // 2)
            self.dna_sequence_group[1].append((flag, dna_sequence))

            location = random.randint(0, self.dna_length - 2)
            dna_sequence = list(self.reference)
            dna_sequence[int(location) + 0] = "G"
            dna_sequence[int(location) + 1] = "C"
            dna_sequence[int(location) + 2] = "C"
            dna_sequence = "".join(dna_sequence)

            flag = "GCC" not in dna_sequence[-self.observed_length:]
            self.dna_sequence_group[2].append((flag, dna_sequence))

    def test(self):
        for dna_sequences, bio_filter in zip(self.dna_sequence_group, self.bio_filters):
            for flag, dna_sequence in dna_sequences:
                self.assertEqual(flag, bio_filter.valid(dna_sequence=dna_sequence, only_last=True))
