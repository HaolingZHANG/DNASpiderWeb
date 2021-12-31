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

from dsw import encode, Monitor, bit_to_number, calculus_division

from experiments import colors, create_folders, obtain_filters

filterwarnings(action="ignore", category=MatplotlibDeprecationWarning)


def proposed(dataset, shape, index_length, random_seed):
    print("Calculate generated algorithms.")
    records, monitor, current, total = [], Monitor(), 0, 12 * len(dataset) * 100  # number of vertex.

    for filter_index in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
        accessor = load(file="./results/data/a" + filter_index + "[g].npy")
        vertices = where(load(file="./results/data/a" + filter_index + "[v].npy") == 1)[0]
        random.seed(random_seed)
        random.shuffle(vertices)
        chosen_indices = vertices[:100]
        for data_index, data in enumerate(dataset):
            matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
            for start_index in chosen_indices:
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                nucleotide_number, monitor = 0, Monitor()
                for index, binary_message in enumerate(matrix):
                    dna_string = encode(binary_message=binary_message, accessor=accessor, start_index=start_index)
                    nucleotide_number += len(dna_string)
                    monitor.output(index + 1, len(matrix))
                records.append([0, int(filter_index), nucleotide_number])
                current += 1

    return records


def hedges(dataset, shape, index_length):

    def trans(binary_messages, screen, mapping):
        # as the Section "Imposing DNA Output Sequence Constraints" in HEDGES work.
        dna_strings, monitor = [], Monitor()
        for bit_index, bit_segment in enumerate(binary_messages):
            quotient, dna_string = bit_to_number(bit_segment, is_string=True), ""
            available_nucleotides = deepcopy(mapping)
            while quotient != "0":
                if len(available_nucleotides) > 1:  # current vertex contains information
                    quotient, remainder = calculus_division(number=quotient, base=str(len(available_nucleotides)))
                    nucleotide = available_nucleotides[int(remainder)]
                else:  # current vertex does not contain information.
                    nucleotide = available_nucleotides[0]

                dna_string += nucleotide

                available_nucleotides = []
                for potential_nucleotide in mapping:
                    if screen.valid(dna_string + potential_nucleotide, only_last=True):  # check the choice of C#.
                        available_nucleotides.append(potential_nucleotide)

                if len(available_nucleotides) == 0:
                    return []

            dna_strings.append(dna_string)
            monitor.output(bit_index + 1, len(binary_messages))

        return dna_strings

    print("Calculate HEDGES Code.")
    records, current, total = [], 0, 12 * 24 * 100

    for filter_index, bio_filter in obtain_filters().items():
        for data_index, data in enumerate(dataset):
            for value in permutations("ACGT", 4):
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
                obtain_dna_strings = trans(binary_messages=matrix, screen=bio_filter, mapping=list(value))
                if len(obtain_dna_strings) > 0:
                    nucleotide_number = sum([len(dna_string) for dna_string in obtain_dna_strings])
                else:
                    nucleotide_number = -1
                records.append([1, int(filter_index), nucleotide_number])
                current += 1

    return records


def fountain(dataset, shape, index_length):

    class Code(DNAFountain):

        def trans(self, binary_messages, screen, false_threshold):
            lfsr, used_seeds, dna_strings, monitor = DNAFountain.LFSR().lfsr_s_p(), dict(), [], Monitor()
            check_binary_messages, done_binary_messages = [None] * len(binary_messages), set()
            chunk_to_droplets = defaultdict(set)
            self.prng = DNAFountain.PRNG(number=len(binary_messages), delta=self.delta, c=self.c_dist)
            false_times = 0
            while True:
                seed = next(lfsr)
                if seed in used_seeds:
                    continue
                droplet = DNAFountain.Droplet()
                dna_string = "".join(droplet.get_dna(seed, self.prng, binary_messages, self.header_size))
                if screen.valid(dna_string=dna_string, only_last=False):  # calculate local biochemical constraints
                    dna_strings.append(dna_string)

                    # decode synchronously.
                    droplet = DNAFountain.Droplet()
                    droplet.init_binaries(self.prng, dna_string, self.header_size)
                    for chunk_num in droplet.chuck_indices:
                        chunk_to_droplets[chunk_num].add(droplet)
                    self.update_droplets(droplet, check_binary_messages, done_binary_messages, chunk_to_droplets)
                    false_times = 0
                else:
                    false_times += 1

                monitor.output(len(done_binary_messages), len(binary_messages))

                if len(done_binary_messages) == len(binary_messages):
                    return dna_strings

                if false_times >= false_threshold:
                    monitor.output(len(binary_messages), len(binary_messages))
                    return []

    header_size = index_length // 8
    parameter_set = random.uniform(low=array([0.1, 0.01]), high=array([1, 0.1]), size=(100, 2))

    print("Calculate DNA Fountain.")
    records, current, total = [], 0, 12 * len(dataset) * len(parameter_set)

    for filter_index, bio_filter in obtain_filters().items():
        for data_index, data in enumerate(dataset):
            for delta, c in parameter_set:
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                matrix = data.reshape(shape).tolist()
                method = Code(c_dist=c, delta=delta, header_size=header_size)
                obtained_dna_strings = method.trans(binary_messages=matrix, screen=bio_filter, false_threshold=1e6)
                if len(obtained_dna_strings) > 0:
                    nucleotide_number = sum([len(dna_string) for dna_string in obtained_dna_strings])
                else:
                    nucleotide_number = -1
                records.append([2, int(filter_index), nucleotide_number])
                current += 1

    return records


def yin_yang(dataset, shape, index_length, random_seed):

    class Code(YinYangCode):

        def trans(self, bit_segments, screen):
            total_count, dna_sequences, monitor = len(bit_segments), [], Monitor()
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

                        if screen.valid(dna_string=dna_sequence[0], only_last=False):  # calculate local constraints.
                            is_finish = True
                            dna_sequences.append(dna_sequence[0])
                            del bit_segments[selected_index]
                            break

                        if screen.valid(dna_string=dna_sequence[1], only_last=False):  # calculate local constraints.
                            is_finish = True
                            dna_sequences.append(dna_sequence[1])
                            del bit_segments[selected_index]
                            break

                if not is_finish:
                    dna_sequences.append("A" * len(fixed_bit_segment))  # for faster transcoding (getting results).

                monitor.output(total_count - len(bit_segments), total_count)

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
        print("Get Yin-Yang Code rules.")
        rules = []
        temp_rule1 = ["".join(x) for x in product("01", repeat=4)]
        temp_rule2 = ["".join(x) for x in product("01", repeat=16)]
        x, y, monitor = 1, 4 * len(temp_rule1) * len(temp_rule2), Monitor()
        for base in ["A", "T", "C", "G"]:
            for rule1index in range(len(temp_rule1)):
                for rule2index in range(len(temp_rule2)):
                    rule1 = list(map(int, list(temp_rule1[rule1index])))
                    rule2 = array(list(map(int, list(temp_rule2[rule2index])))).reshape(4, 4).tolist()
                    if check(rule1, rule2):
                        rules.append((rule1, rule2, base))
                    monitor.output(x, y)
                    x += 1

        return rules

    total_rules = get_yyc_rules()

    random.seed(random_seed)

    print("Calculate Yin-Yang Code.")
    total_indices = array(list(range(len(total_rules))))
    random.shuffle(total_indices)
    chosen_indices = total_indices[:100]
    records, current, total = [], 0, 12 * len(dataset) * len(chosen_indices)

    for filter_index, bio_filter in obtain_filters().items():
        for data_index, data in enumerate(dataset):
            matrix = array(insert_index(data=data.reshape(shape).tolist(), index_length=index_length))
            for rule_index in chosen_indices:
                print("task (" + str(current + 1) + " / " + str(total) + ").")
                rule = total_rules[rule_index]
                random.seed(random_seed)
                method = Code(yang_rule=rule[0], yin_rule=rule[1], virtual_nucleotide=rule[2])
                obtained_dna_strings = method.trans(bit_segments=matrix.tolist(), screen=bio_filter)
                nucleotide_number = sum([len(dna_string) for dna_string in obtained_dna_strings])
                records.append([3, int(filter_index), nucleotide_number])
                current += 1

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
        records += fountain(dataset=dataset, shape=(2048, 256), index_length=32)
        records += yin_yang(dataset=dataset, shape=(4096, 128), index_length=16, random_seed=2021)
        save(file="./results/data/step_3_stability_evaluation.npy", arr=array(records))


def draw():
    display_data = [[[] for _ in range(12)] for _ in range(4)]
    for sample in load(file="./results/data/step_3_stability_evaluation.npy"):
        if sample[2] > 0:
            display_data[int(sample[0])][int(sample[1]) - 1].append((8 * 1024 * 64) / sample[2])
        else:
            display_data[int(sample[0])][int(sample[1]) - 1].append(-1)

    pyplot.figure(figsize=(10, 5), tight_layout=True)
    bias = [0.375, 0.125, -0.375, -0.125]
    used_colors = [[colors["algo1"], colors["algo2"]], [colors["hedc1"], colors["hedc2"]],
                   [colors["foun1"], colors["foun2"]], [colors["yyco1"], colors["yyco2"]]]
    labels = ["SPIDER-WEB", "HEDGES", "DNA Fountain", "Yin-Yang Code"]
    for bias_index, method_data in enumerate(display_data):
        for location, data in enumerate(method_data):
            if sum(data) == -len(data):
                pyplot.scatter([location + bias[bias_index]], [0], color=used_colors[bias_index][0], marker="x", s=15)
            else:
                if max(data) - min(data) < 0.02:
                    result = median(data)
                    pyplot.hlines(result, location + bias[bias_index] - 0.08, location + bias[bias_index] + 0.08,
                                  linewidths=1, edgecolors=used_colors[bias_index][0], zorder=3)
                    pyplot.scatter([location + bias[bias_index]], result,
                                   color="white", edgecolor=used_colors[bias_index][0], linewidth=1, s=5, zorder=4)
                else:
                    violin = pyplot.violinplot(dataset=data, positions=[location + bias[bias_index]],
                                               bw_method=0.5, showextrema=False, widths=0.16)
                    for patch in violin["bodies"]:
                        patch.set_edgecolor(used_colors[bias_index][0])
                        patch.set_facecolor(used_colors[bias_index][1])
                        patch.set_linewidth(1)
                        patch.set_alpha(1)
                    pyplot.scatter([location + bias[bias_index]], median(data),
                                   color="white", edgecolor=used_colors[bias_index][0], linewidth=1, s=5, zorder=4)

            if location % 2 != 0:
                pyplot.fill_between([location - 0.5, location + 0.5], [-0.08, -0.08], [2.08, 2.08],
                                    color=colors["diffs"], zorder=0)

    legends = [patches.Patch(facecolor=used_colors[index][1], edgecolor=used_colors[index][0],
                             linewidth=1, label=labels[index]) for index in [0, 1, 3, 2]]
    pyplot.legend(handles=legends, loc="upper left", fontsize=8)
    pyplot.xlim(-0.5, 11.5)
    pyplot.ylim(-0.08, 2.08)
    pyplot.xticks(range(12), range(1, 13), fontsize=8)
    pyplot.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
                  ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0"], fontsize=8)
    pyplot.xlabel("constraint set", fontsize=8)
    pyplot.ylabel("code rate", fontsize=8)
    pyplot.savefig("./results/figures/[3-1] stability evaluation.svg", format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()
    whole_data = generate(random_seed=2021, total_length=8 * 1024 * 64, times=100)
    evaluate(dataset=whole_data)
    draw()
