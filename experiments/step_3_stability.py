__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


# noinspection PyPackageRequirements
from Chamaeleo.utils.indexer import connect_all
# noinspection PyPackageRequirements
from Chamaeleo.methods.flowed import YinYangCode, DNAFountain
from collections import defaultdict
from copy import deepcopy
from itertools import product, permutations
# noinspection PyPackageRequirements
from matplotlib import pyplot, patches, MatplotlibDeprecationWarning
from numpy import random, array, load, save, sum, median, min, max, where
from os import path
from warnings import filterwarnings

from dsw import encode, Monitor, bit_to_number, calculus_division, LocalBioFilter

from experiments import colors, create_folders, obtain_filters
from experiments import balanced_gc_range, accepted_gc_bias, low_gc_bias, high_gc_bias, cut_segments, nanopore_segments

filterwarnings(action="ignore", category=MatplotlibDeprecationWarning)


def proposed(dataset, shape, index_length, random_seed):
    print("Calculate generated algorithms.")
    records, monitor, current, total = [], Monitor(), 0, 12 * 10000

    for filter_index in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
        encoded_graph = load(file="./results/data/a" + filter_index + "[g].npy")
        start_indices = where(load(file="./results/data/a" + filter_index + "[v].npy") == 1)[0]
        random.seed(random_seed)
        random.shuffle(start_indices)
        chosen_indices = start_indices[:100]
        for data_index, data in enumerate(dataset):
            bit_number = len(data)
            matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
            for start_index in chosen_indices:
                nucleotide_number = 0
                for binary_message in matrix:
                    dna_string = encode(binary_message=binary_message, accessor=encoded_graph, start_index=start_index)
                    nucleotide_number += len(dna_string)
                records.append([0, int(filter_index), data_index + 1, bit_number, nucleotide_number])
                current += 1
                monitor.output(current, total)

    return records


def yin_yang(dataset, shape, index_length, random_seed):

    class Code(YinYangCode):
        def trans(self, bit_segments, screen):
            total_count, dna_sequences = len(bit_segments), []
            while len(bit_segments) > 0:
                fixed_bit_segment, is_finish = bit_segments.pop(), False
                for pair_time in range(self.max_iterations):
                    if len(bit_segments) > 0:
                        selected_index = random.randint(0, len(bit_segments))
                        selected_bit_segment, dna_sequence = bit_segments[selected_index], [[], []]
                        support_nucleotide_1, support_nucleotide_2 = self.virtual_nucleotide, self.virtual_nucleotide
                        for bit_1, bit_2 in zip(fixed_bit_segment, selected_bit_segment):
                            current_nucleotide_1 = self._bits_to_nucleotide(bit_1, bit_2, support_nucleotide_1)
                            current_nucleotide_2 = self._bits_to_nucleotide(bit_2, bit_1, support_nucleotide_2)
                            dna_sequence[0].append(current_nucleotide_1)
                            dna_sequence[1].append(current_nucleotide_2)
                            support_nucleotide_1, support_nucleotide_2 = current_nucleotide_1, current_nucleotide_2
                        dna_sequence = ["".join(dna_sequence[0]), "".join(dna_sequence[1])]
                        access = True  # calculate local biochemical constraints
                        for sub_position in range(len(dna_sequence[0]) - 9):
                            if not screen.valid(dna_sequence[0][sub_position: sub_position + 10]):
                                access = False
                                break
                        if access:
                            is_finish = True
                            dna_sequences.append(dna_sequence[0])
                            del bit_segments[selected_index]
                            break
                        access = True  # calculate local biochemical constraints
                        for sub_position in range(len(dna_sequence[1]) - 9):
                            if not screen.valid(dna_sequence[1][sub_position: sub_position + 10]):
                                access = False
                                break
                        if access:
                            is_finish = True
                            dna_sequences.append(dna_sequence[1])
                            del bit_segments[selected_index]
                            break
                if not is_finish:
                    dna_sequences.append("A" * len(fixed_bit_segment))

            return dna_sequences

    def check(rule1, rule2):
        for index in range(len(rule1)):
            if rule1[index] != 0 and rule1[index] != 1:
                return False
        if sum(rule1) != 2:
            return False
        if rule1[0] == rule1[1]:
            same = [0, 1, 2, 3]
        elif rule1[0] == rule1[2]:
            same = [0, 2, 1, 3]
        else:
            same = [0, 3, 1, 2]
        for row in range(len(rule2)):
            if rule2[row][same[0]] + rule2[row][same[1]] != 1 or rule2[row][same[0]] * rule2[row][same[1]] != 0:
                return False
            if rule2[row][same[2]] + rule2[row][same[3]] != 1 or rule2[row][same[2]] * rule2[row][same[3]] != 0:
                return False

        return True

    def get_yyc_rules():
        rules = []
        temp_rule1 = ["".join(x) for x in product("01", repeat=4)]
        temp_rule2 = ["".join(x) for x in product("01", repeat=16)]
        for base in ["A", "T", "C", "G"]:
            for rule1index in range(len(temp_rule1)):
                for rule2index in range(len(temp_rule2)):
                    rule1 = list(map(int, list(temp_rule1[rule1index])))
                    rule2 = array(list(map(int, list(temp_rule2[rule2index])))).reshape(4, 4).tolist()
                    if check(rule1, rule2):
                        rules.append((rule1, rule2, base))

        return rules

    random.seed(random_seed)

    print("Calculate Yin-Yang Code.")
    total_rules = get_yyc_rules()
    total_indices = array(list(range(len(total_rules))))
    random.shuffle(total_indices)
    chosen_indices = total_indices[:100]
    records, monitor, current, total = [], Monitor(), 0, 12 * 10000

    for filter_index, bio_filter in obtain_filters().items():
        for data_index, data in enumerate(dataset):
            bit_number = len(data)
            matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
            for rule_index in chosen_indices:
                rule = total_rules[rule_index]
                random.seed(random_seed)
                method = Code(yang_rule=rule[0], yin_rule=rule[1], virtual_nucleotide=rule[2])
                dna_strings = method.trans(bit_segments=matrix.tolist(), screen=bio_filter)
                nucleotide_number = sum([len(dna_string) for dna_string in dna_strings])
                records.append([1, int(filter_index), data_index + 1, bit_number, nucleotide_number])
                current += 1
                monitor.output(current, total)

    return records


def fountain(dataset, shape, index_length):

    class Code(DNAFountain):
        def trans(self, bit_segments, screen):
            lfsr, used_seeds, dna_sequences = DNAFountain.LFSR().lfsr_s_p(), dict(), []
            check_bit_segments, done_segments, chunk_to_droplets = [None] * len(bit_segments), set(), defaultdict(set)
            self.prng = DNAFountain.PRNG(number=len(bit_segments), delta=self.delta, c=self.c_dist)
            false_times = 0
            while True:
                seed = next(lfsr)
                if seed in used_seeds:
                    continue
                droplet = DNAFountain.Droplet()
                dna_sequence = "".join(droplet.get_dna(seed, self.prng, bit_segments, self.header_size))
                access = True  # calculate local biochemical constraints
                for sub_position in range(len(dna_sequence) - 9):
                    if not screen.valid(dna_sequence[sub_position: sub_position + 10]):
                        access = False
                        break
                if access:
                    dna_sequences.append(dna_sequence)

                    # decode synchronously.
                    droplet = DNAFountain.Droplet()
                    droplet.init_binaries(self.prng, dna_sequence, self.header_size)
                    for chunk_num in droplet.chuck_indices:
                        chunk_to_droplets[chunk_num].add(droplet)
                    self.update_droplets(droplet, check_bit_segments, done_segments, chunk_to_droplets)
                    false_times = 0
                else:
                    false_times += 1

                if len(done_segments) == len(bit_segments):
                    return dna_sequences

                if false_times >= 10000:
                    return []

    header_size = index_length // 8
    parameter_set = random.uniform(low=array([0.1, 0.01]), high=array([1, 0.1]), size=(100, 2))

    print("Calculate DNA Fountain.")
    records, monitor, current, total = [], Monitor(), 0, 12 * 10000

    for filter_index, bio_filter in obtain_filters().items():
        for data_index, data in enumerate(dataset):
            for delta, c in parameter_set:
                bit_number = len(data)
                matrix = data.reshape(shape).tolist()
                method = Code(c_dist=c, delta=delta, header_size=header_size)
                dna_strings = method.trans(bit_segments=matrix, screen=bio_filter)
                if len(dna_strings) > 0:
                    nucleotide_number = sum([len(dna_string) for dna_string in dna_strings])
                else:
                    nucleotide_number = -1
                records.append([2, int(filter_index), data_index + 1, bit_number, nucleotide_number])
                current += 1
                monitor.output(current, total)

    return records


def hedges(dataset, shape, index_length):

    class ReplacedFilter(LocalBioFilter):

        def __init__(self, observed_length, max_homopolymer_runs=None, gc_range=None, undesired_motifs=None):
            super().__init__(max_homopolymer_runs=max_homopolymer_runs,
                             gc_range=gc_range, undesired_motifs=undesired_motifs)
            self._observed_length = observed_length

        def valid(self, dna_string):
            observed_dna = dna_string[-self._observed_length:]

            for nucleotide in observed_dna:
                if nucleotide not in "ACGT":
                    return False

            if self._max_homopolymer_runs is not None:
                for nucleotide in "ACGT":
                    if nucleotide * (1 + self._max_homopolymer_runs) in observed_dna:
                        return False

            if self._gc_range is not None:
                gc_count = float(observed_dna.count("C") + observed_dna.count("G"))
                if gc_count > self._gc_range[1] * self._observed_length:
                    return False
                at_count = float(observed_dna.count("A") + observed_dna.count("T"))
                if at_count > (1 - self._gc_range[0]) * self._observed_length:
                    return False

            if self._undesired_motifs is not None:
                for special in self._undesired_motifs:
                    if special in observed_dna:
                        return False
                    reverse_complement = special.replace("A", "t").replace("C", "g")
                    reverse_complement = reverse_complement.replace("G", "c").replace("T", "a")
                    reverse_complement = reverse_complement[::-1].upper()
                    if reverse_complement in observed_dna:
                        return False

            return True

    def transcode(bit_segments, screen, mapping):
        dna_sequences, m = [], Monitor()
        for bit_index, bit_segment in enumerate(bit_segments):
            decimal_number, dna_sequence = bit_to_number(bit_segment, is_string=False), ""
            available_flag, available_nucleotides = True, deepcopy(mapping)
            # print(decimal_number)
            while decimal_number > 0:
                if len(available_nucleotides) > 1:  # current vertex contains information
                    decimal_number, remainder = calculus_division(number=decimal_number,
                                                                  base=str(len(available_nucleotides)))
                    nucleotide = available_nucleotides[int(remainder)]
                else:  # current vertex does not contain information.
                    nucleotide = available_nucleotides[0]

                dna_sequence += nucleotide

                available_nucleotides = []
                for potential_nucleotide in mapping:
                    if screen.valid(dna_sequence + potential_nucleotide):
                        available_nucleotides.append(potential_nucleotide)

                if len(available_nucleotides) == 0:
                    available_flag = False
                    break

            if not available_flag:
                dna_sequence = ""

            dna_sequences.append(dna_sequence)
            m.output(bit_index + 1, len(bit_segments))

        return dna_sequences

    filters = {
        "01": ReplacedFilter(max_homopolymer_runs=2, gc_range=balanced_gc_range, undesired_motifs=cut_segments,
                             observed_length=10),
        "02": ReplacedFilter(max_homopolymer_runs=1,
                             observed_length=10),  # Goldman method
        "03": ReplacedFilter(gc_range=low_gc_bias,
                             observed_length=10),
        "04": ReplacedFilter(max_homopolymer_runs=2, gc_range=accepted_gc_bias, undesired_motifs=nanopore_segments,
                             observed_length=10),
        "05": ReplacedFilter(max_homopolymer_runs=2, gc_range=accepted_gc_bias,
                             observed_length=10),
        "06": ReplacedFilter(gc_range=high_gc_bias,
                             observed_length=10),
        "07": ReplacedFilter(max_homopolymer_runs=3, gc_range=accepted_gc_bias,
                             observed_length=10),
        "08": ReplacedFilter(max_homopolymer_runs=4, gc_range=accepted_gc_bias,
                             observed_length=10),  # DNA Fountain and Yin-Yang Code
        "09": ReplacedFilter(max_homopolymer_runs=3,
                             observed_length=10),  # Church method
        "10": ReplacedFilter(max_homopolymer_runs=4,
                             observed_length=10),
        "11": ReplacedFilter(max_homopolymer_runs=5,
                             observed_length=10),
        "12": ReplacedFilter(max_homopolymer_runs=6,
                             observed_length=10)
    }

    print("Calculate HEDGES Code.")
    records, monitor, current, total = [], Monitor(), 0, 12 * 24 * 100

    for filter_index, bio_filter in filters.items():
        for data_index, data in enumerate(dataset):
            for value in permutations("ACGT", 4):
                print(current, total)
                bit_number = len(data)
                matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
                dna_strings = transcode(bit_segments=matrix, screen=bio_filter, mapping=list(value))
                if len(dna_strings) > 0:
                    nucleotide_number = sum([len(dna_string) for dna_string in dna_strings])
                else:
                    nucleotide_number = -1
                print([3, int(filter_index), data_index + 1, bit_number, nucleotide_number])
                records.append([3, int(filter_index), data_index + 1, bit_number, nucleotide_number])
                current += 1
                # monitor.output(current, total)

    return records


def generate(random_seed, total_length, times):
    dataset = []
    random.seed(random_seed)
    for _ in range(times):
        dataset.append(random.randint(low=0, high=2, size=(total_length,), dtype=int))
    return dataset


def insert_index(data, index_length):
    return connect_all(bit_segments=data, index_binary_length=index_length, need_logs=False)[0]


def evaluate(dataset):
    if not path.exists("./results/data/step_3_stability_evaluation.npy"):
        records = []
        records += proposed(dataset=dataset, shape=(4096, 128), index_length=16, random_seed=2021)
        records += hedges(dataset=dataset, shape=(4096, 128), index_length=16)
        records += yin_yang(dataset=dataset, shape=(4096, 128), index_length=16, random_seed=2021)
        records += fountain(dataset=dataset, shape=(2048, 256), index_length=32)
        save(file="./results/data/step_3_stability_evaluation.npy", arr=array(records))


def draw():
    display_data = [[[] for _ in range(12)] for _ in range(3)]
    for sample in load(file="./results/data/step_3_stability_evaluation.npy"):
        if sample[3] > 0:
            display_data[int(sample[0])][int(sample[1]) - 1].append((8 * 1024 * 64) / sample[3])
        else:
            display_data[int(sample[0])][int(sample[1]) - 1].append(-1)

    pyplot.figure(figsize=(10, 5))
    bias = [0.3, 0, -0.3]
    used = [[colors["algo1"], colors["algo2"]], [colors["yyco1"], colors["yyco2"]], [colors["foun1"], colors["foun2"]]]
    labels = ["SPIDER-WEB", "Yin-Yang Code", "DNA Fountain"]
    for bias_index, method_data in enumerate(display_data):
        for location, data in enumerate(method_data):
            if sum(data) == -len(data):
                pyplot.scatter([location + bias[bias_index]], [0], color=used[bias_index][0], marker="x", s=10)
            elif max(data) - min(data) < 0.01:
                result = median(data)
                pyplot.hlines(result, location + bias[bias_index] - 0.1, location + bias[bias_index] + 0.1,
                              linewidths=1, edgecolors=used[bias_index][0], zorder=3)
                pyplot.scatter([location + bias[bias_index]], result,
                               color="white", edgecolor=used[bias_index][0], linewidth=1, s=5, zorder=4)
            else:
                violin = pyplot.violinplot(dataset=data, positions=[location + bias[bias_index]],
                                           bw_method=0.5, showextrema=False, widths=0.2)
                for patch in violin["bodies"]:
                    patch.set_edgecolor(used[bias_index][0])
                    patch.set_facecolor(used[bias_index][1])
                    patch.set_linewidth(1)
                    patch.set_alpha(1)
                pyplot.scatter([location + bias[bias_index]], median(data),
                               color="white", edgecolor=used[bias_index][0], linewidth=1, s=5, zorder=4)

            if location % 2 != 0:
                pyplot.fill_between([location - 0.5, location + 0.5], [-0.08, -0.08], [2.08, 2.08],
                                    color=colors["diffs"], zorder=0)
    legends = [patches.Patch(facecolor=used[index][1], edgecolor=used[index][0], linewidth=1, label=labels[index])
               for index in range(3)][::-1]
    pyplot.legend(handles=legends, loc="upper left", fontsize=8)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylim(-0.08, 2.08)
    pyplot.xticks(range(12), range(1, 13), fontsize=8)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
                  ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0"], fontsize=8)
    pyplot.xlabel("constraint set", fontsize=8)
    pyplot.ylabel("code rate", fontsize=8)

    pyplot.savefig("./results/figures/[3-1] stability evaluation.png", format="png", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()
    # whole_data = generate(random_seed=2021, total_length=8 * 1024 * 64, times=100)
    # evaluate(dataset=whole_data)
    draw()
