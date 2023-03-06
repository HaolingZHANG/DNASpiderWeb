from Chamaeleo.methods.flowed import YinYangCode, DNAFountain
from Chamaeleo.utils.indexer import connect_all
from collections import defaultdict
from itertools import product
from numpy import random, array, sum, min, max, where, longlong
from warnings import filterwarnings

from dsw import encode, Monitor, number_to_bit, obtain_vertices, bit_to_number

from experiments import local_bio_filters, load_data

filterwarnings("ignore", category=RuntimeWarning)


def generate_database(random_seed, total_length, times):
    dataset = []
    random.seed(random_seed)
    for _ in range(times):
        dataset.append(random.randint(low=0, high=2, size=(total_length,), dtype=int))
    return dataset


def insert_index(data, index_length):
    return connect_all(bit_segments=data, index_binary_length=index_length, need_logs=False)[0]


def trans_fountain(dataset, shape, index_length, random_seed, param_number):
    class Code(DNAFountain):

        def trans(self, bit_arrays, screen, false_threshold):
            lfsr, used_seeds, dna_strings, monitor = DNAFountain.LFSR().lfsr_s_p(), dict(), [], Monitor()
            check_bit_arrays, done_bit_arrays, chunk_to_droplets = [None] * len(bit_arrays), set(), defaultdict(set)
            self.prng = DNAFountain.PRNG(number=len(bit_arrays), delta=self.delta, c=self.c_dist)
            false_times = 0

            while True:
                seed = next(lfsr)
                if seed in used_seeds:
                    continue

                droplet = DNAFountain.Droplet()
                dna_string = "".join(droplet.get_dna(seed, self.prng, bit_arrays, self.header_size))

                if screen.valid(dna_string=dna_string, only_last=False):  # calculate local biochemical constraints
                    false_times = 0
                    dna_strings.append(dna_string)

                    # decode synchronously.
                    droplet = DNAFountain.Droplet()
                    droplet.init_binaries(self.prng, dna_string, self.header_size)
                    for chunk_num in droplet.chuck_indices:
                        chunk_to_droplets[chunk_num].add(droplet)
                    self.update_droplets(droplet, check_bit_arrays, done_bit_arrays, chunk_to_droplets)
                else:
                    false_times += 1

                monitor(len(done_bit_arrays), len(bit_arrays))

                if len(done_bit_arrays) == len(bit_arrays):
                    return dna_strings

                if false_times >= false_threshold:
                    monitor(len(bit_arrays), len(bit_arrays))
                    return []

    print("Calculate DNA Fountain.")
    random.seed(random_seed)
    header_size = index_length // 8
    parameter_set = random.uniform(low=array([0.1, 0.01]), high=array([1, 0.1]), size=(param_number, 2))
    records, current, total = [], 0, 12 * len(dataset) * param_number

    for filter_index, bio_filter in local_bio_filters.items():
        for parameter_index, (delta, c) in enumerate(parameter_set):
            for data_index, data in enumerate(dataset):
                matrix = data.reshape(shape)
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                method = Code(c_dist=c, delta=delta, header_size=header_size)
                obtained_dna_strings = method.trans(bit_arrays=matrix.tolist(), screen=bio_filter, false_threshold=1e6)

                if len(obtained_dna_strings) > 0:
                    nucleotide_number = sum([len(dna_string) for dna_string in obtained_dna_strings])
                    records.append([1, int(filter_index), parameter_index, data_index,
                                    (8 * 1024 * 64) / float(nucleotide_number)])
                else:
                    records.append([1, int(filter_index), parameter_index, data_index,
                                    -1])

                current += 1

    return records


def trans_yinyang(dataset, shape, index_length, random_seed, param_number):
    class Code(YinYangCode):

        def trans(self, bit_arrays, screen):
            total_count, dna_strings, monitor = len(bit_arrays), [], Monitor()
            while len(bit_arrays) > 0:
                fixed_bit_segment, is_finish = bit_arrays.pop(), False
                for pair_time in range(self.max_iterations):
                    if len(bit_arrays) > 0:
                        selected_index = random.randint(0, len(bit_arrays))
                        selected_bit_array, dna_string = bit_arrays[selected_index], [[], []]
                        support_nucleotide_1, support_nucleotide_2 = self.virtual_nucleotide, self.virtual_nucleotide

                        for bit_1, bit_2 in zip(fixed_bit_segment, selected_bit_array):
                            current_nucleotide_1 = self._bits_to_nucleotide(bit_1, bit_2, support_nucleotide_1)
                            current_nucleotide_2 = self._bits_to_nucleotide(bit_2, bit_1, support_nucleotide_2)
                            dna_string[0].append(current_nucleotide_1)
                            dna_string[1].append(current_nucleotide_2)
                            support_nucleotide_1, support_nucleotide_2 = current_nucleotide_1, current_nucleotide_2
                        dna_string = ["".join(dna_string[0]), "".join(dna_string[1])]

                        if screen.valid(dna_string=dna_string[0], only_last=False):  # calculate local constraints.
                            is_finish = True
                            dna_strings.append(dna_string[0])
                            del bit_arrays[selected_index]
                            break

                        if screen.valid(dna_string=dna_string[1], only_last=False):  # calculate local constraints.
                            is_finish = True
                            dna_strings.append(dna_string[1])
                            del bit_arrays[selected_index]
                            break

                if not is_finish:
                    # ignore random strategy for faster transcoding (getting show).
                    dna_strings.append("A" * len(fixed_bit_segment))

                monitor(total_count - len(bit_arrays), total_count)

            return dna_strings

    # load Yin-Yang Code rule based on the rule index.
    def find_rule_parameters(index):
        ref_1 = [3, 5, 6, 9, 10, 12]
        ref_2 = [[5, 6, 9, 10], [3, 6, 9, 12], [3, 5, 10, 12], [3, 5, 10, 12], [3, 6, 9, 12], [5, 6, 9, 10]]
        ref_3 = ["A", "T", "C", "G"]
        yang = number_to_bit(decimal_number=ref_1[(index % 1536) // 256], bit_length=4)
        yin = [number_to_bit(decimal_number=ref_2[index % 1536 // 256][(index % 256) // 64], bit_length=4),
               number_to_bit(decimal_number=ref_2[index % 1536 // 256][(index % 64) // 16], bit_length=4),
               number_to_bit(decimal_number=ref_2[index % 1536 // 256][(index % 16) // 4], bit_length=4),
               number_to_bit(decimal_number=ref_2[index % 1536 // 256][(index % 4) // 1], bit_length=4)]
        virtual = ref_3[index // 1536]

        return yang, yin, virtual

    print("Calculate Yin-Yang Code.")
    random.seed(random_seed)
    total_indices = array(list(range(6144)))
    random.shuffle(total_indices)
    chosen_indices = total_indices[:param_number]
    records, current, total = [], 0, 12 * len(dataset) * len(chosen_indices)

    for filter_index, bio_filter in local_bio_filters.items():
        for param_index, rule_index in enumerate(chosen_indices):
            for data_index, data in enumerate(dataset):
                matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                random.seed(random_seed)
                yang_rule, yin_rule, virtual_nucleotide = find_rule_parameters(index=rule_index)
                method = Code(yang_rule=yang_rule, yin_rule=yin_rule, virtual_nucleotide=virtual_nucleotide)
                obtained_dna_strings = method.trans(bit_arrays=matrix.tolist(), screen=bio_filter)
                nucleotide_number = sum([len(dna_string) for dna_string in obtained_dna_strings])
                records.append([2, int(filter_index) - 1, param_index, data_index,
                                (8 * 1024 * 64) / float(nucleotide_number)])
                current += 1

    return records


def trans_hedges(dataset, shape, index_length):
    # a better implementation is shown in https://github.com/HaolingZHANG/pyHEDGES.
    salt_number, previous_number, low_order_number = 46, 8, 10
    patterns = [[2, 1], [2, 1, 1, 1, 1], [1], [1, 1, 0], [1, 0], [1, 0, 0]]
    correct_penalties = [-0.035, -0.082, -0.127, -0.229, -0.265, -0.324]

    def hash_function(source_value):
        """
        Obtain the target value from the source value based on the well-accepted hash function.

        :param source_value: source bit value.
        :type source_value: int

        :return: target value after the hash function.
        :rtype: int
        """
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

    class HypothesisNode:

        def __init__(self, pattern_flag, message, string):
            self.pattern_flag, self.message, self.string = pattern_flag, message, string

        def next(self, whole_string, nucleotide_index, current_score, pattern, used_filter, salt_value,
                 correct_penalty, mutate_penalty, insert_penalty, delete_penalty):
            follow_vertices, follow_scores, follow_indices = [], [], []

            # collect the available nucleotides in this location.
            available_nucleotides = []
            for potential_nucleotide in ["A", "C", "G", "T"]:
                if used_filter.valid(self.string + potential_nucleotide, only_last=True):
                    available_nucleotides.append(potential_nucleotide)

            if len(available_nucleotides) == 0:  # this path is blocked, stop running.
                return [], [], [], []

            # low-order q bits of the bit position index i.
            location_value = len(self.message) % (2 ** low_order_number)

            if len(self.message) - previous_number >= 0:  # p previous concatenated bits.
                previous_info = self.message[len(self.message) - previous_number:]
                previous_value = bit_to_number(previous_info, is_string=False)
            else:
                previous_value = 0

            for message_bit in product([0, 1], repeat=pattern[self.pattern_flag]):
                hash_value = hash_function(location_value | previous_value | salt_value)
                bit_value = bit_to_number(list(message_bit), is_string=False) if len(message_bit) > 0 else 0
                nucleotide = available_nucleotides[(hash_value + bit_value) % len(available_nucleotides)]
                message, string = self.message + list(message_bit), self.string + nucleotide
                node = HypothesisNode((self.pattern_flag + 1) % len(pattern), message, string)
                if nucleotide == whole_string[nucleotide_index]:  # assume that current nucleotide is correct.
                    follow_vertices.append(node)
                    follow_scores.append(current_score + correct_penalty)
                    follow_indices.append(nucleotide_index + 1)
                else:
                    # assume that current nucleotide is mutated.
                    follow_vertices.append(node)
                    follow_scores.append(current_score + mutate_penalty)
                    follow_indices.append(nucleotide_index + 1)

                    # assume that current nucleotide is inserted, the (i + 1)-th nucleotide is i-th nucleotide.
                    if nucleotide_index + 1 < len(whole_string) and nucleotide == whole_string[nucleotide_index + 1]:
                        follow_vertices.append(node)
                        follow_scores.append(current_score + insert_penalty)
                        follow_indices.append(nucleotide_index + 2)

                    # assume that current nucleotide is deleted.
                    follow_vertices.append(node)
                    follow_scores.append(current_score + delete_penalty)
                    follow_indices.append(nucleotide_index)

            return follow_vertices, follow_scores, follow_indices, [len(v.message) for v in follow_vertices]

    def hedges_encode(binary_message, strand_index, pattern, used_filter):
        dna_string, available_nucleotides, bit_location, pattern_flag = "", ["A", "C", "G", "T"], 0, 0
        salt_value = strand_index % (2 ** salt_number)  # s bits of salt (strand ID).
        while bit_location < len(binary_message):
            # low-order q bits of the bit position index i.
            location_value = bit_location % (2 ** low_order_number)

            if bit_location - previous_number >= 0:
                previous_info = binary_message[bit_location - previous_number: bit_location]
                previous_value = bit_to_number(previous_info, is_string=False)
            else:
                previous_value = 0

            hash_value = hash_function(location_value | previous_value | salt_value)
            if len(available_nucleotides) > 0:
                bit_number = pattern[pattern_flag]
                message_bit = binary_message[bit_location: bit_location + bit_number]
                bit_value = bit_to_number(message_bit, is_string=False) if len(message_bit) > 0 else 0
                nucleotide = available_nucleotides[(hash_value + bit_value) % len(available_nucleotides)]
                bit_location += bit_number
                pattern_flag = (pattern_flag + 1) % len(pattern)
            else:
                raise ValueError("DNA string (index = " + str(strand_index) + ") "
                                 + "cannot be encoded because of the established constraints!")

            dna_string += nucleotide

            available_nucleotides = []
            for potential_nucleotide in ["A", "C", "G", "T"]:
                if used_filter.valid(dna_string + potential_nucleotide, only_last=True):
                    available_nucleotides.append(potential_nucleotide)

        return dna_string

    def hedges_decode(dna_string, strand_index, bit_length, pattern, used_filter, initial_score, heap_limitation,
                      correct_penalty, mutate_penalty, insert_penalty, delete_penalty):
        monitor, terminal_indices, heap_size = Monitor(), None, 1
        heap = {"v": [HypothesisNode(0, [], "")], "s": [initial_score], "i": [0], "l": [0]}  # priority heap

        while True:  # repair by A star search (score priority).
            chuck_indices, chuck_score = where(array(heap["s"]) == min(heap["s"]))[0], min(heap["s"])

            for chuck_index in chuck_indices:
                # noinspection PyTypeChecker
                heap["s"][chuck_index] = len(dna_string) + 1  # set the chuck vertex to inaccessible.
                heap_size -= 1

                if heap["i"][chuck_index] >= len(dna_string):
                    continue

                # noinspection PyUnresolvedReferences
                follow_info = heap["v"][chuck_index].next(dna_string, heap["i"][chuck_index], chuck_score,
                                                          pattern, used_filter, strand_index % (2 ** salt_number),
                                                          correct_penalty, mutate_penalty, insert_penalty,
                                                          delete_penalty)
                heap_size += len(follow_info[0])

                heap["v"], heap["l"] = heap["v"] + follow_info[0], heap["l"] + follow_info[3]
                heap["s"], heap["i"] = heap["s"] + follow_info[1], heap["i"] + follow_info[2]

                current_process = heap["i"][chuck_index] + 1
                monitor(current_process, len(dna_string), extra={"size": len(heap["v"]), "score": "%.2f" % chuck_score})

                # the first chain of hypotheses to decode the required bytes of message wins.
                if bit_length == max(heap["l"]) or heap_size >= heap_limitation:
                    if current_process < len(dna_string):
                        print()  # not finish.

                    candidates = []
                    for terminal_index in where(array(heap["l"]) == bit_length)[0]:
                        candidates.append((heap["v"][terminal_index].message, heap["v"][terminal_index].string))

                    return candidates, heap_size, current_process

    records, current, total = [], 0, 12 * len(dataset) * 6
    for filter_index, bio_filter in local_bio_filters.items():
        for param_index, (p, c_penalty) in enumerate(zip(patterns, correct_penalties)):
            for data_index, data in enumerate(dataset):
                matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                nucleotide_number, miss_flag = 0, False
                for index, bit_array in enumerate(matrix):
                    encoded_string = hedges_encode(binary_message=bit_array, strand_index=index,
                                                   pattern=p, used_filter=bio_filter)
                    nucleotide_number += len(encoded_string)
                    decoded_info = hedges_decode(dna_string=encoded_string, strand_index=index,
                                                 bit_length=len(bit_array), pattern=p, used_filter=bio_filter,
                                                 initial_score=0, heap_limitation=1e5, correct_penalty=c_penalty,
                                                 mutate_penalty=1.0, insert_penalty=1.0, delete_penalty=1.0)
                    if encoded_string not in decoded_info[0]:
                        miss_flag = True
                        break

                if miss_flag:
                    records.append([3, int(filter_index) - 1, param_index, data_index,
                                    (8 * 1024 * 64) / float(nucleotide_number)])
                else:
                    records.append([3, int(filter_index) - 1, param_index, data_index,
                                    -1])

                current += 1

    return records


def trans_spiderweb(dataset, shape, index_length, random_seed, param_number):
    print("Calculate generated algorithms.")
    records, current, total = [], 0, 12 * len(dataset) * param_number  # number of vertex.
    coding_graphs = load_data(load_path="./data/graph_coding.pkl")

    for filter_index in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
        accessor = coding_graphs[filter_index]
        vertices = obtain_vertices(accessor=accessor)
        random.seed(random_seed)
        random.shuffle(vertices)
        chosen_indices = vertices[:param_number]
        for param_index, start_index in enumerate(chosen_indices):
            for data_index, data in enumerate(dataset):
                matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                nucleotide_number, monitor = 0, Monitor()

                for index, binary_message in enumerate(matrix):
                    dna_string = encode(binary_message=matrix, accessor=accessor, start_index=start_index)
                    nucleotide_number += len(dna_string)
                    monitor(index + 1, len(matrix))

                records.append([4, int(filter_index) - 1, param_index, data_index,
                                (8 * 1024 * 64) / float(nucleotide_number)])
                current += 1

    return records


def conversion_step(accessor, dna_length, start_index, repeats, random_seed):
    nucleotides = "ACGT"

    random.seed(random_seed)

    def multiplication(number, base):
        if base == "0":
            return "0", 1

        if base == "1":
            return number, 1

        number = [int(item) for item in number]

        remainder, times = 0, 0
        for index in range(len(number))[::-1]:
            current, times = number[index] * int(base) + remainder, times + 1
            if current >= 10:
                number[index] = current % 10
                remainder = current // 10
            else:
                number[index] = current
                remainder = 0

        while remainder > 0:
            number.insert(0, remainder % 10)
            remainder //= 10
            times += 1

        result = "".join(list(map(str, number)))

        return result, times

    def addition(number, base):
        number, base = list(number), list(base.zfill(len(number)))

        result, times = [0 for _ in range(len(number) + 1)], 0
        for index in range(len(number) - 1, -1, -1):
            sum_value = int(number[index]) + int(base[index]) + int(result[index + 1])
            if sum_value < 10:
                result[index + 1] = sum_value
                times += 1
            else:
                flag = 0
                while sum_value > 0:
                    result[index + 1 - flag] = sum_value % 10
                    sum_value //= 10
                    flag += 1
                    times += 1

        result = "".join(list(map(str, result)))

        return result if result[0] != "0" else result[1:], times

    def division(number, base):
        if base == "0":
            return "0", "0", 1

        if base == "1":
            return number, "0", 1

        if len(number) == 1 and number[0] < base:
            return "0", number[0], 1

        number, new_number, remainder, times = [int(item) for item in number], [], 0, 0
        for index, quotient in enumerate(number):
            current = quotient + remainder * 10
            if current >= int(base):
                new_number.append(current // int(base))
                remainder = current - new_number[-1] * int(base)
            else:
                new_number.append(0)
                remainder = current
            times += 1

        quotient = "".join(list(map(str, new_number)))

        for index in range(len(quotient)):
            times += 1
            if quotient[index] != "0":
                return quotient[index:], str(remainder), times

        return "0", str(remainder), times

    def to_number(dna_sequence):
        nucleotide_values, times = list(map(nucleotides.index, dna_sequence)), 0

        decimal_number = "0"
        for nucleotide_value in nucleotide_values:
            # multiply by length of usage of nucleotides.
            decimal_number, sub_times = multiplication(number=decimal_number, base=str(len(nucleotides)))
            times += sub_times
            # add current nucleotide value.
            decimal_number, sub_times = addition(number=decimal_number, base=str(nucleotide_value))
            times += sub_times

        return decimal_number, times

    def to_binary(decimal_number):
        one_array, times = [], 0
        while decimal_number != "0":  # decimal number > 0
            decimal_number, remainder, sub_times = division(number=decimal_number, base="2")
            one_array.insert(0, int(remainder))
            times += sub_times

        return one_array, times

    results, monitor = [], Monitor()
    for repeat in range(repeats):
        # create ideal DNA string randomly.
        vertex_index, right_dna_string = start_index, ""
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            right_dna_string += nucleotide

        intermediate_value, t_1 = to_number(dna_sequence=right_dna_string)
        _, t_2 = to_binary(decimal_number=intermediate_value)
        results.append([t_1, t_2])
        monitor(repeat + 1, repeats, extra={"step 1": t_1, "step 2": t_2})

    random.seed(None)

    return array(results)
