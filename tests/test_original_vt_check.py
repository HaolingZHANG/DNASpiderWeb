from numpy import random, array, sum, all
from unittest import TestCase


def repair_by_vt(wrong_dna_values, wrong_vt_value, right_vt_value):
    results = []
    for assume_location, nucleotide_value in enumerate(wrong_dna_values):
        for repaired_value in list(filter(lambda n: n != nucleotide_value, [0, 1, 2, 3])):
            if (wrong_vt_value - right_vt_value) == (nucleotide_value - repaired_value) * (assume_location + 1):
                repaired_dna_values = wrong_dna_values.copy()
                repaired_dna_values[assume_location] = repaired_value
                results.append(repaired_dna_values)

    for assume_location, nucleotide_value in enumerate(wrong_dna_values):
        assume_difference = nucleotide_value * (assume_location + 1) + sum(wrong_dna_values[assume_location + 1:])
        if wrong_vt_value - right_vt_value == assume_difference:
            repaired_dna_values = wrong_dna_values.tolist()
            del repaired_dna_values[assume_location]
            results.append(array(repaired_dna_values))

    for assume_location in range(len(wrong_dna_values)):
        for nucleotide_value in [0, 1, 2, 3]:
            assume_difference = - nucleotide_value * (assume_location + 1) - sum(wrong_dna_values[assume_location:])
            if wrong_vt_value - right_vt_value == assume_difference:
                repaired_dna_values = wrong_dna_values.tolist()
                repaired_dna_values.insert(assume_location, nucleotide_value)
                results.append(array(repaired_dna_values))

    return results


class TestVT(TestCase):

    def setUp(self):
        random.seed(2021)
        self.repeats = 100
        self.dna_length = 20
        self.dna_value_groups = random.randint(0, 4, size=(self.repeats, self.dna_length))

    def test(self):
        for right_dna_values in self.dna_value_groups:
            right_vt_value = 0
            for index, value in enumerate(right_dna_values):
                right_vt_value += value * (index + 1)
            right_vt_value %= 4 ** 10
            error_location = len(right_dna_values) // 2

            # 3 substitution error.
            for error_value in list(filter(lambda n: n != int(right_dna_values[error_location]), [0, 1, 2, 3])):
                wrong_dna_values = right_dna_values.copy()
                wrong_dna_values[len(right_dna_values) // 2] = error_value

                wrong_vt_value = 0
                for index, value in enumerate(wrong_dna_values):
                    wrong_vt_value += value * (index + 1)
                wrong_vt_value %= 4 ** 10

                repair_dna_value_groups, found = repair_by_vt(wrong_dna_values, wrong_vt_value, right_vt_value), False
                for repair_dna_values in repair_dna_value_groups:
                    if len(repair_dna_values) == len(right_dna_values) and all(repair_dna_values == repair_dna_values):
                        found = True
                        break

                self.assertEqual(found, True)

            # 4 insertion error.
            for error_value in [0, 1, 2, 3]:
                wrong_dna_values = right_dna_values.tolist()
                wrong_dna_values.insert(error_location, error_value)
                wrong_dna_values = array(wrong_dna_values)

                wrong_vt_value = 0
                for index, value in enumerate(wrong_dna_values):
                    wrong_vt_value += value * (index + 1)
                wrong_vt_value %= 4 ** 10

                repair_dna_value_groups, found = repair_by_vt(wrong_dna_values, wrong_vt_value, right_vt_value), False
                for repair_dna_values in repair_dna_value_groups:
                    if len(repair_dna_values) == len(right_dna_values) and all(repair_dna_values == repair_dna_values):
                        found = True
                        break

                self.assertEqual(found, True)

            # 1 deletion error.
            wrong_dna_values = right_dna_values.tolist()
            del wrong_dna_values[error_location]
            wrong_dna_values = array(wrong_dna_values)

            wrong_vt_value = 0
            for index, value in enumerate(wrong_dna_values):
                wrong_vt_value += value * (index + 1)
            wrong_vt_value %= 4 ** 10

            repair_dna_value_groups, found = repair_by_vt(wrong_dna_values, wrong_vt_value, right_vt_value), False
            for repair_dna_values in repair_dna_value_groups:
                if len(repair_dna_values) == len(right_dna_values) and all(repair_dna_values == repair_dna_values):
                    found = True
                    break

            self.assertEqual(found, True)
