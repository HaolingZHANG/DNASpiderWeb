# noinspection PyPackageRequirements
from Chamaeleo.methods.flowed import YinYangCode, DNAFountain
# noinspection PyPackageRequirements
from Chamaeleo.utils.indexer import connect_all
from collections import defaultdict
from numpy import random, sum

from dsw import encode, Monitor, number_to_bit


def generated(binary_messages, accessor, start_index):
    nucleotide_number, monitor = 0, Monitor()

    for current, binary_message in enumerate(binary_messages):
        dna_string = encode(binary_message=binary_message, accessor=accessor, start_index=start_index)
        nucleotide_number += len(dna_string)
        monitor.output(current + 1, len(binary_messages))

    return nucleotide_number


def hedges(binary_messages, mapping, bio_filter):
    nucleotide_number, monitor = 0, Monitor()

    # consider the out-degree 3 is 2 and ignore the hash function for rough implementation.
    for current, binary_message in enumerate(binary_messages):
        dna_string, available_nucleotides, bit_location = "", mapping, 0
        while bit_location < len(binary_message):
            if len(available_nucleotides) == 1:
                nucleotide = available_nucleotides[0]
            elif len(available_nucleotides) == 2 or len(available_nucleotides) == 3:
                value = binary_message[bit_location]
                nucleotide = available_nucleotides[value]
                bit_location += 1
            else:
                if bit_location + 2 < len(binary_message):
                    value = binary_message[bit_location] * 2 + binary_message[bit_location + 1]
                else:
                    value = binary_message[bit_location]

                nucleotide = available_nucleotides[value]
                bit_location += 2

            dna_string += nucleotide

            available_nucleotides = []
            for potential_nucleotide in mapping:
                if bio_filter.valid(dna_string + potential_nucleotide, only_last=True):
                    available_nucleotides.append(potential_nucleotide)

            if len(available_nucleotides) == 0:
                return -1

        nucleotide_number += len(dna_string)
        monitor.output(current + 1, len(binary_messages))

    return nucleotide_number


def yinyang(binary_messages, bio_filter, rule_index):

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
                    # ignore random strategy for faster transcoding (getting results).
                    dna_strings.append("A" * len(fixed_bit_segment))

                monitor.output(total_count - len(bit_arrays), total_count)

            return dna_strings

    # load Yin-Yang Code rule based on the rule index.
    def find_rule_param(index):
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

    yang_rule, yin_rule, virtual_nucleotide = find_rule_param(index=rule_index)
    method = Code(yang_rule=yang_rule, yin_rule=yin_rule, virtual_nucleotide=virtual_nucleotide)
    obtained_dna_strings = method.trans(bit_arrays=binary_messages, screen=bio_filter)

    return sum([len(dna_string) for dna_string in obtained_dna_strings])


def fountain(binary_messages, bio_filter, header_size, c, delta):

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

                monitor.output(len(done_bit_arrays), len(bit_arrays))

                if len(done_bit_arrays) == len(bit_arrays):
                    return dna_strings

                if false_times >= false_threshold:
                    monitor.output(len(bit_arrays), len(bit_arrays))
                    return []

    method = Code(c_dist=c, delta=delta, header_size=header_size)
    obtained_dna_strings = method.trans(bit_arrays=binary_messages, screen=bio_filter, false_threshold=1e6)

    if len(obtained_dna_strings) > 0:
        nucleotide_number = sum([len(dna_string) for dna_string in obtained_dna_strings])
    else:
        nucleotide_number = -1

    return nucleotide_number


def generate(random_seed, total_length, times):
    dataset = []
    random.seed(random_seed)
    for _ in range(times):
        dataset.append(random.randint(low=0, high=2, size=(total_length,), dtype=int))
    return dataset


def insert_index(data, index_length):
    return connect_all(bit_segments=data, index_binary_length=index_length, need_logs=False)[0]
