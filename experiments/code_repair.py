from itertools import product
from time import time
from numpy import random, array, arange, zeros, abs, min, max, log, where, all, longlong
from warnings import filterwarnings

from dsw import set_vt, repair_dna, obtain_vertices, bit_to_number
from dsw import Monitor, accessor_to_latter_map, approximate_capacity, remove_nasty_arc

filterwarnings("ignore", category=RuntimeWarning)


def introduce_errors(right_dna_sequence, errors):
    if errors == 0:
        return right_dna_sequence

    nucleotides = ["A", "C", "G", "T"]

    total_length = len(right_dna_sequence)
    while True:
        wrong_dna_sequence = list(right_dna_sequence)  # introduce error under most dense.
        locations = []

        for _ in range(errors):
            wrong_type, location, success = random.randint(0, 3), None, False
            candidates = arange(11, len(wrong_dna_sequence) - 11, dtype=int)
            random.shuffle(candidates)
            for location in candidates:
                # uniform distribution detection
                if len(locations) == 0 or min(abs(array(locations) - location)) > total_length / (errors + 2):
                    success = True
                    locations.append(location)
                    break

            if not success:
                wrong_dna_sequence, locations, wrong_type = list(right_dna_sequence), [], -1
                break

            if wrong_type == 0:  # substitution
                nucleotide = wrong_dna_sequence[location]
                wrong_dna_sequence[location] = random.choice(list(filter(lambda n: n != nucleotide, nucleotides)))
            elif wrong_type == 1:  # insertion
                wrong_dna_sequence.insert(location, random.choice(nucleotides))
            else:  # deletion
                del wrong_dna_sequence[location]

        wrong_dna_sequence = "".join(wrong_dna_sequence)

        if wrong_dna_sequence != right_dna_sequence:
            return wrong_dna_sequence


def evaluate_repair_errors(random_seed, accessor, vertices, observed_length, repeats, dna_length, errors):
    nucleotides = ["A", "C", "G", "T"]

    random.seed(random_seed)

    records, sequence_data, count, monitor = [], [], 0, Monitor()
    for current_repeat in range(repeats):
        start_index = random.choice(vertices)

        # create ideal DNA string randomly.
        vertex_index, right_dna_sequence = start_index, ""
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            right_dna_sequence += nucleotide

        vt_check = set_vt(dna_sequence=right_dna_sequence, vt_length=observed_length + 1, nucleotides=nucleotides)

        wrong_dna_sequence = introduce_errors(right_dna_sequence=right_dna_sequence, errors=errors,
                                              nucleotides=nucleotides)
        sequence_data.append((right_dna_sequence, start_index, vt_check, wrong_dna_sequence))
        monitor(current_repeat + 1, repeats)

    previous_time = time()
    for index in range(len(sequence_data)):
        right_dna_sequence, start_index, vt_check, wrong_dna_sequence = sequence_data[index]
        repaired_strings, params = repair_dna(dna_sequence=wrong_dna_sequence, accessor=accessor,
                                              start_index=start_index, observed_length=observed_length,
                                              vt_check=vt_check, has_indel=True)
        found_flag = right_dna_sequence in repaired_strings
        records.append([params[0], params[1], params[2], len(repaired_strings), found_flag])
        count += found_flag
        monitor(index + 1, len(sequence_data), extra={"count": count})

    used_times = time() - previous_time
    random.seed(None)

    return array(records), used_times


def calculate_minimum_reads(random_seed, accessor, vertices, observed_length, repeats, dna_length, errors):
    nucleotides = ["A", "C", "G", "T"]

    random.seed(random_seed)
    records, monitor, maximum_reads = [], Monitor(), 1
    for current_repeat in range(repeats):
        start_index = random.choice(vertices)

        # create ideal DNA string randomly.
        vertex_index, right_dna_sequence = start_index, ""
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            right_dna_sequence += nucleotide

        vt_check = set_vt(dna_sequence=right_dna_sequence, vt_length=observed_length + 1, nucleotides=nucleotides)

        reads, result_set = 1, {}
        while True:
            wrong_dna_sequence = introduce_errors(right_dna_sequence=right_dna_sequence, errors=errors,
                                                  nucleotides=nucleotides)

            repaired_strings, _ = repair_dna(dna_sequence=wrong_dna_sequence, accessor=accessor,
                                             start_index=start_index, observed_length=observed_length,
                                             vt_check=vt_check, has_indel=True)

            for repaired_dna_sequence in repaired_strings:
                if repaired_dna_sequence in result_set:
                    result_set[repaired_dna_sequence] += 1
                else:
                    result_set[repaired_dna_sequence] = 1

            strings, counts = [], []
            for string, count in result_set.items():
                strings.append(string)
                counts.append(count)
            counts = array(counts)

            if len(counts) > 0 and len(where(counts == max(counts))[0]) == 1:
                best_dna_sequence = strings[where(counts == max(counts))[0][0]]
                if best_dna_sequence == right_dna_sequence:
                    maximum_reads = max([maximum_reads, reads])
                    break

            reads += 1

        monitor(current_repeat + 1, repeats, extra={"current reads": reads, "maximum reads": maximum_reads})

        records.append(reads)

    random.seed(None)

    return records


def evaluate_established_depth(random_seed, accessor, vertices, observed_length, classes, dna_length,
                               errors, established_depth):
    nucleotides = ["A", "C", "G", "T"]

    random.seed(random_seed)
    right_set, obtained_set, monitor = [], {}, Monitor()
    for index in range(classes):
        while True:  # create ideal DNA string randomly.
            start_index = random.choice(vertices)
            vertex_index, right_dna_sequence = start_index, ""
            for location in range(dna_length):
                used_indices = where(accessor[vertex_index] >= 0)[0]
                used_index = random.choice(used_indices)
                nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
                right_dna_sequence += nucleotide

            if right_dna_sequence not in right_set:
                right_set.append(right_dna_sequence)
                break

        vt_check = set_vt(dna_sequence=right_dna_sequence, vt_length=observed_length + 1, nucleotides=nucleotides)

        for _ in range(established_depth):
            wrong_dna_sequence = introduce_errors(right_dna_sequence=right_dna_sequence, errors=errors)

            repaired_strings, _ = repair_dna(dna_sequence=wrong_dna_sequence, accessor=accessor,
                                             start_index=start_index, observed_length=observed_length,
                                             vt_check=vt_check, has_indel=True)

            for repaired_dna_sequence in repaired_strings:
                if repaired_dna_sequence in obtained_set:
                    obtained_set[repaired_dna_sequence] += 1
                else:
                    obtained_set[repaired_dna_sequence] = 1

        monitor(index + 1, classes)

    random.seed(None)

    return right_set, obtained_set


def screen_edges_for_repair(accessor, maximum_rounds=0):
    print("Obtain latter map from the accessor.")
    latter_map = accessor_to_latter_map(accessor=accessor, verbose=True)
    code_rate = approximate_capacity(accessor=accessor, repeats=4, maximum_iteration=100, verbose=False)
    print("Original code rate is %.5f.\n" % float(code_rate))

    current, results = 1, []
    while True:
        accessor, latter_map, arc, scores = remove_nasty_arc(accessor=accessor, latter_map=latter_map,
                                                             iteration=current, verbose=True,
                                                             has_insertion=True, has_deletion=True)
        new_code_rate = approximate_capacity(accessor=accessor, repeats=4, maximum_iteration=100, verbose=False)

        print("Current code rate is %.5f.\n" % float(new_code_rate))

        vertices = obtain_vertices(accessor=accessor)
        correction_data = {}
        for error_time in [1, 2, 3, 4]:
            records, runtimes = evaluate_repair_errors(random_seed=2021, error_times=error_time,
                                                       accessor=accessor, vertices=vertices,
                                                       observed_length=10, repeats=2000,
                                                       dna_length=100, check_iterations=error_time + 2)
            correction_data[(100, error_time)] = (records, runtimes)

        results.append((arc, new_code_rate, correction_data))

        # The observed length is 10, therefore, at least 10 nucleotides encode one bit.
        if new_code_rate < log(4.0) / log(len(accessor)) or 0 < maximum_rounds == current:
            break

        current += 1

    return results


def spiderweb_step(random_seed, accessor, vertices, dna_length, observed_length, repeats, errors):
    nucleotides = ["A", "C", "G", "T"]

    print("Check DNA sequence of length " + str(dna_length) + " imposing " + str(errors) + " error(s).")
    random.seed(random_seed)
    records, count, monitor = zeros(shape=(repeats,), dtype=int), 0, Monitor()
    for repeat in range(repeats):
        start_index = random.choice(vertices)

        # create ideal DNA sequence randomly.
        vertex_index, right_dna_sequence = start_index, ""
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            right_dna_sequence += nucleotide

        vt_check = set_vt(dna_sequence=right_dna_sequence, vt_length=observed_length + 1, nucleotides=nucleotides)

        wrong_dna_sequence = introduce_errors(right_dna_sequence=right_dna_sequence, errors=errors)

        repaired_sequences, params = repair_dna(dna_sequence=wrong_dna_sequence, accessor=accessor,
                                                start_index=start_index, observed_length=observed_length,
                                                vt_check=vt_check, has_indel=True)
        records[count] = params[-1]
        monitor(count, repeats)

    return records


def hedges_step(random_seed, bio_filter, dna_length, repeats, errors):
    nucleotides = ["A", "C", "G", "T"]

    def hash_function(source_value):
        target_value = array(source_value, dtype=longlong) * array(3935559000370003845, dtype=longlong)
        target_value += array(2691343689449507681, dtype=longlong)
        target_value ^= target_value >> 21
        target_value ^= target_value << 37
        target_value ^= target_value >> 4
        target_value *= array(4768777513237032717, dtype=longlong)
        target_value ^= target_value << 20
        target_value ^= target_value >> 41
        target_value ^= target_value << 5
        return target_value

    def hedges_encode(binary_message, strand_index, pattern, salt_number=46, previous_number=8, low_order_number=10):
        dna_sequence, available_nucleotides, bit_location, pattern_flag = "", nucleotides, 0, 0
        salt_index = strand_index % (2 ** salt_number)  # s bits of salt (strand ID).

        while True:
            try:
                while bit_location < len(binary_message):
                    bit_index = bit_location % (2 ** low_order_number)  # low-order q bits of the bit position index i.

                    if bit_location - previous_number >= 0:
                        previous_info = binary_message[bit_location - previous_number: bit_location]
                        previous_value = bit_to_number(previous_info, is_string=False)
                    else:
                        previous_value = 0

                    hash_value = hash_function(bit_index | previous_value | salt_index)
                    if len(available_nucleotides) > 0:
                        bit_number = pattern[pattern_flag]
                        message_bit = binary_message[bit_location: bit_location + bit_number]
                        bit_value = bit_to_number(message_bit, is_string=False) if len(message_bit) > 0 else 0
                        nucleotide = available_nucleotides[(hash_value + bit_value) % len(available_nucleotides)]
                        bit_location += bit_number
                        pattern_flag = (pattern_flag + 1) % len(pattern)
                    else:
                        raise ValueError("DNA sequence (index = " + str(strand_index) + ") "
                                         + "cannot be encoded because of the established constraints!")

                    dna_sequence += nucleotide

                    available_nucleotides = []
                    for potential_nucleotide in nucleotides:
                        if bio_filter.valid(dna_sequence + potential_nucleotide, only_last=True):
                            available_nucleotides.append(potential_nucleotide)
            except ValueError:
                dna_sequence = ""

            if len(dna_sequence) > 0:
                break

        return dna_sequence

    def hedges_decode(dna_sequence, strand_index, bit_length, pattern,
                      salt_number=46, previous_number=8, low_order_number=10, initial_score=0.0,
                      correct_penalty=-0.035, insert_penalty=1.0, delete_penalty=1.0, mutate_penalty=1.0):
        # s bits of salt (strand index).
        salt_value = strand_index % (2 ** salt_number)

        class HypothesisNode:

            def __init__(self, pattern_flag, message, sequence):
                self.pattern_flag, self.message, self.sequence = pattern_flag, message, sequence

            def next(self, nucleotide_index, current_score):
                follow_vertices, follow_scores, follow_indices = [], [], []

                # collect the available nucleotides in this location.
                available_nucleotides = []
                for potential_nucleotide in nucleotides:
                    if bio_filter.valid(self.sequence + potential_nucleotide, only_last=True):
                        available_nucleotides.append(potential_nucleotide)

                if len(available_nucleotides) == 0:  # this path is blocked, stop running.
                    return [], [], [], []

                # low-order q bits of the bit position index i.
                bit_index = len(self.message) % (2 ** low_order_number)

                if len(self.message) - previous_number >= 0:  # p previous concatenated bits.
                    previous_info = self.message[len(self.message) - previous_number:]
                    previous_value = bit_to_number(previous_info, is_string=False)
                else:
                    previous_value = 0

                for message_bit in product([0, 1], repeat=pattern[self.pattern_flag]):
                    hash_value = hash_function(bit_index | previous_value | salt_value)
                    bit_value = bit_to_number(list(message_bit), is_string=False) if len(message_bit) > 0 else 0
                    nucleotide = available_nucleotides[(hash_value + bit_value) % len(available_nucleotides)]
                    message, sequence = self.message + list(message_bit), self.sequence + nucleotide

                    node = HypothesisNode((self.pattern_flag + 1) % len(pattern), message, sequence)
                    if nucleotide == dna_sequence[nucleotide_index]:  # assume that current nucleotide is correct.
                        follow_vertices.append(node)
                        follow_scores.append(current_score + correct_penalty)
                        follow_indices.append(nucleotide_index + 1)
                    else:
                        # assume that current nucleotide is mutated.
                        follow_vertices.append(node)
                        follow_scores.append(current_score + mutate_penalty)
                        follow_indices.append(nucleotide_index + 1)

                        # assume that current nucleotide is inserted, the (i + 1)-th nucleotide is i-th nucleotide.
                        if nucleotide_index + 1 < len(dna_sequence) and \
                                nucleotide == dna_sequence[nucleotide_index + 1]:
                            follow_vertices.append(node)
                            follow_scores.append(current_score + insert_penalty)
                            follow_indices.append(nucleotide_index + 2)

                        # assume that current nucleotide is deleted.
                        follow_vertices.append(HypothesisNode((self.pattern_flag + 1) % len(pattern),
                                                              message, sequence))
                        follow_scores.append(current_score + delete_penalty)
                        follow_indices.append(nucleotide_index)

                # noinspection PyShadowingNames
                return follow_vertices, follow_scores, follow_indices, [len(v.message) for v in follow_vertices]

        terminal_indices, heap_size, search_step, monitor = None, 1, 0, Monitor()
        heap = {"v": [HypothesisNode(0, [], "")], "s": [initial_score], "i": [0], "l": [0]}  # priority heap

        while True:  # repair by A star search (score priority).
            chuck_indices, chuck_score = where(array(heap["s"]) == min(heap["s"]))[0], min(heap["s"])

            for chuck_index in chuck_indices:
                # noinspection PyTypeChecker
                heap["s"][chuck_index] = len(dna_sequence) + 1  # set the chuck vertex to inaccessible.
                heap_size -= 1

                if heap["i"][chuck_index] >= len(dna_sequence):
                    continue

                follow_info = heap["v"][chuck_index].next(heap["i"][chuck_index], chuck_score)
                heap_size, search_step = heap_size + len(follow_info[0]), search_step + len(follow_info)
                heap["v"], heap["l"] = heap["v"] + follow_info[0], heap["l"] + follow_info[3]
                heap["s"], heap["i"] = heap["s"] + follow_info[1], heap["i"] + follow_info[2]
                heap["v"][chuck_index], current_process = None, heap["i"][chuck_index] + 1
                monitor(current_process, len(dna_sequence),
                        extra={"size": heap_size, "score": "%.2f" % chuck_score, "search step": search_step})
                # the first chain of hypotheses to decode the required bytes of message wins.
                if bit_length == max(heap["l"]):
                    if current_process < len(dna_sequence):
                        print()  # not finish.

                    candidates = []
                    for terminal_index in where(array(heap["l"]) == bit_length)[0]:
                        candidates.append(array(heap["v"][terminal_index].message, dtype=int))

                    return candidates, heap_size, search_step, current_process

    random.seed(seed=random_seed)
    records = zeros(shape=(repeats,), dtype=int)
    for repeat in range(repeats):
        print("Run test " + str(repeat + 1) + "/" + str(repeats) + " with " + str(errors) + " errors.")
        used_message = random.randint(0, 2, size=(dna_length // 3 + 1,))
        right_dna_sequence = hedges_encode(binary_message=used_message, strand_index=repeat, pattern=[1, 0, 0])
        wrong_dna_sequence = introduce_errors(right_dna_sequence=right_dna_sequence, errors=errors)
        availables, size, step, process = hedges_decode(dna_sequence=wrong_dna_sequence, strand_index=repeat,
                                                        bit_length=len(used_message), pattern=[1, 0, 0],
                                                        correct_penalty=-0.324)
        for available in availables:
            if all(available[:-1] == used_message[:-1]):
                break
        records[repeat] = step

    return records
