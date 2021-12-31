__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from collections import Counter
from itertools import combinations
# noinspection PyPackageRequirements
from matplotlib import pyplot
from numpy import random, load, save, zeros, ones, ones_like, array, linspace, concatenate
from numpy import union1d, mean, sum, cumsum, log, log10, where
from pickle import dump as psave
from pickle import load as pload
from os import path, remove

from dsw import Monitor, obtain_leaf_vertices, number_to_dna, accessor_to_latter_map, encode
from dsw import remove_useless, approximate_capacity, repair_dna

from experiments import colors, create_folders


def repair_examples():

    def substitution_example(accessor, vertices):
        start_index = where(vertices == 1)[0][0]
        binary_message = [0, 1] * 20
        wrong_location = 14
        right_dna_string = encode(binary_message=binary_message, accessor=accessor, start_index=start_index)
        for nucleotide in list(filter(lambda n: n != right_dna_string[wrong_location], ["A", "C", "G", "T"])):
            wrong_dna_string = list(right_dna_string)
            wrong_dna_string[wrong_location] = nucleotide
            wrong_dna_string = "".join(wrong_dna_string)
            print("Right: " + right_dna_string)
            print("Wrong: " + wrong_dna_string)

            repaired_dna_strings = repair_dna(dna_string=wrong_dna_string, accessor=accessor, start_index=start_index,
                                              has_insertion=True, has_deletion=True, verbose=True)
            if repaired_dna_strings is not None:
                print()
                print("Repair strategy:")
                for repair_type, repair_info in repaired_dna_strings.items():
                    if len(repair_info) > 0:
                        for index, data in enumerate(repair_info):
                            occur_location, used_nucleotide, obtained_dna_string = data
                            result = str(obtained_dna_string == right_dna_string).rjust(5) + " | " + obtained_dna_string
                            if index > 0:
                                print("                       = " + result)
                            else:
                                if repair_type == "S":
                                    print("repair by substitution = " + result)
                                elif repair_type == "I":
                                    print("repair by insertion    = " + result)
                                else:
                                    print("repair by deletion     = " + result)
            else:
                print("Cannot find any repair strategy.")
            print()
            print()

    def insertion_example(accessor, vertices):
        start_index = where(vertices == 1)[0][0]
        binary_message = [0, 1] * 20
        wrong_location = 14
        right_dna_string = encode(binary_message=binary_message, accessor=accessor, start_index=start_index)
        for nucleotide in list(["A", "C", "G", "T"]):
            wrong_dna_string = list(right_dna_string)
            wrong_dna_string.insert(wrong_location, nucleotide)
            wrong_dna_string = "".join(wrong_dna_string)
            print("Right: " + right_dna_string)
            print("Wrong: " + wrong_dna_string)

            repaired_dna_strings = repair_dna(dna_string=wrong_dna_string, accessor=accessor, start_index=start_index,
                                              has_insertion=True, has_deletion=True, verbose=True)
            if repaired_dna_strings is not None:
                print()
                print("Repair strategy:")
                for repair_type, repair_info in repaired_dna_strings.items():
                    if len(repair_info) > 0:
                        for index, data in enumerate(repair_info):
                            occur_location, used_nucleotide, obtained_dna_string = data
                            result = str(obtained_dna_string == right_dna_string).rjust(5) + " | " + obtained_dna_string
                            if index > 0:
                                print("                       = " + result)
                            else:
                                if repair_type == "S":
                                    print("repair by substitution = " + result)
                                elif repair_type == "I":
                                    print("repair by insertion    = " + result)
                                else:
                                    print("repair by deletion     = " + result)
            else:
                print("Cannot find any repair strategy.")
            print()
            print()

    def deletion_example(accessor, vertices):
        start_index = where(vertices == 1)[0][0]
        binary_message = [0, 1] * 20
        wrong_location = 14
        right_dna_string = encode(binary_message=binary_message, accessor=accessor, start_index=start_index)
        wrong_dna_string = list(right_dna_string)
        del wrong_dna_string[wrong_location]
        wrong_dna_string = "".join(wrong_dna_string)
        print("Right: " + right_dna_string)
        print("Wrong: " + wrong_dna_string)

        repaired_dna_strings = repair_dna(dna_string=wrong_dna_string, accessor=accessor, start_index=start_index,
                                          has_insertion=True, has_deletion=True, verbose=True)
        if repaired_dna_strings is not None:
            print()
            print("Repair strategy:")
            for repair_type, repair_info in repaired_dna_strings.items():
                if len(repair_info) > 0:
                    for index, data in enumerate(repair_info):
                        occur_location, used_nucleotide, obtained_dna_string = data
                        result = str(obtained_dna_string == right_dna_string).rjust(5) + " | " + obtained_dna_string
                        if index > 0:
                            print("                       = " + result)
                        else:
                            if repair_type == "S":
                                print("Repair by substitution = " + result)
                            elif repair_type == "I":
                                print("Repair by deletion     = " + result)
                            else:
                                print("Repair by insertion    = " + result)
        else:
            print("Cannot find any repair strategy.")
        print()
        print()

    a, v = load(file="./results/data/a01[g].npy"), load(file="./results/data/a01[v].npy")
    substitution_example(a, v)
    insertion_example(a, v)
    deletion_example(a, v)


def do_repair(right_dna_string, wrong_dna_string, used_accessor, used_index):
    repaired_dna_strings = repair_dna(dna_string=wrong_dna_string, accessor=used_accessor, start_index=used_index,
                                      has_insertion=True, has_deletion=True, verbose=False)
    if repaired_dna_strings is not None:
        repaired = [len(repaired_dna_strings["S"]), len(repaired_dna_strings["I"]), len(repaired_dna_strings["D"])]
        repair_flag = False
        for repair_type, repair_info in repaired_dna_strings.items():
            if len(repair_info) > 0:
                for occur_location, used_nucleotide, obtained_dna_string in repair_info:
                    if obtained_dna_string == right_dna_string:
                        repair_flag = True
                        break
            if repair_flag:
                break
        return repaired + [repair_flag]
    else:
        return [0, 0, 0, 0]


def evaluate_repair_mechanism(accessor, vertices, wrong_location, repeats):
    print("Evaluate the repair mechanism.")
    nucleotides = ["A", "C", "G", "T"]
    monitor, current, total, records = Monitor(), 1, repeats * len(vertices), []
    for _ in range(repeats):
        for start_index in vertices:
            # obtain ideal DNA string.
            dna_string, vertex_index = "", start_index
            wrong_flag = False
            for _ in range(wrong_location * 2):
                try:
                    used_indices = where(accessor[vertex_index] >= 0)[0]
                    used_index = random.choice(used_indices)
                    nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
                    dna_string += nucleotide
                except ValueError:
                    wrong_flag = True
                    break
                except IndexError:
                    wrong_flag = True
                    break

            if not wrong_flag:
                # occur substitution (3 types)
                for nucleotide in list(filter(lambda n: n != dna_string[wrong_location], ["A", "C", "G", "T"])):
                    s_dna_string = list(dna_string)
                    s_dna_string[wrong_location] = nucleotide
                    s_dna_string = "".join(s_dna_string)
                    situations = do_repair(wrong_dna_string=s_dna_string, right_dna_string=dna_string,
                                           used_accessor=accessor, used_index=start_index)
                    records.append([current, start_index, 0] + situations)

                # occur insertion (4 types)
                for nucleotide in list(["A", "C", "G", "T"]):
                    i_dna_string = list(dna_string)
                    i_dna_string.insert(wrong_location, nucleotide)
                    i_dna_string = "".join(i_dna_string)
                    situations = do_repair(wrong_dna_string=i_dna_string, right_dna_string=dna_string,
                                           used_accessor=accessor, used_index=start_index)
                    records.append([current, start_index, 1] + situations)

                # occur deletion (1 type)
                d_dna_string = list(dna_string)
                del d_dna_string[wrong_location]
                d_dna_string = "".join(d_dna_string)
                situations = do_repair(wrong_dna_string=d_dna_string, right_dna_string=dna_string,
                                       used_accessor=accessor, used_index=start_index)
                records.append([current, start_index, 2] + situations)

            monitor.output(current, total)
            current += 1

    return array(records).astype(int)


def calculate_performance(latter_map, observed_length, evaluation_repeats, vertex_count,
                          capacity_repeats=4, maximum_iteration=100, nucleotides=None):
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    print("Convert latter map to accessor.")
    monitor = Monitor()
    repaired_accessor = -ones(shape=(len(nucleotides) ** observed_length, len(nucleotides)), dtype=int)
    for index, (former_vertex, latter_vertices) in enumerate(latter_map.items()):
        for latter_vertex in latter_vertices:
            repaired_accessor[former_vertex, latter_vertex % len(nucleotides)] = latter_vertex
        monitor.output(index + 1, len(latter_map))

    capacity = approximate_capacity(accessor=repaired_accessor, repeats=capacity_repeats,
                                    maximum_iteration=maximum_iteration, verbose=True)
    print("Finally, the capacity is approximated as " + str(capacity) + ".")
    # The observed length is 10, therefore, at least 10 nucleotides encode one bit.
    # Otherwise, the Perron–Frobenius theorem is not suitable because the multi-connected branch.
    if capacity < 0.1:
        return None

    counts = [0 for _ in range(len(nucleotides))]
    for index, latters in latter_map.items():
        counts[len(latters) - 1] += 1
    counts = array(counts)

    repaired_vertices = zeros(shape=(len(nucleotides) ** 10,), dtype=bool)
    for vertex_index, vertex in enumerate(repaired_accessor):
        if len(where(vertex >= 0)[0]) >= 2:
            repaired_vertices[vertex_index] = 1

    repaired_vertices = where(repaired_vertices > 0)[0]
    random.shuffle(repaired_vertices)
    repaired_vertices = repaired_vertices[:vertex_count]
    info = evaluate_repair_mechanism(accessor=repaired_accessor, vertices=repaired_vertices,
                                     wrong_location=observed_length + 5, repeats=evaluation_repeats)

    strategy_numbers = sum(info[:, 4:7], axis=1)
    strategy_numbers[strategy_numbers == 0] = 1
    print("Average heap level is " + str(mean(strategy_numbers)) + ".")

    return capacity, counts, info


def calculate_maximum_score(latter_map, round_number, nucleotides=None, observed_length=10, repair_indel=True):
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    print("Calculate the score for each directed edge in round " + str(round_number) + ".")
    currents, depth, monitor = list(latter_map.keys()), observed_length - 1, Monitor()
    scores = zeros(shape=(len(nucleotides) ** observed_length, len(nucleotides)), dtype=int)
    for current, current_index in enumerate(currents):
        mutate_branches = []  # substitution
        for latter_index in latter_map[current_index]:
            mutate_branches.append(obtain_leaf_vertices(latter_index, depth, latter_map=latter_map))
        for one, two in combinations(range(len(mutate_branches)), 2):
            score = len(union1d(mutate_branches[one], mutate_branches[two]))
            scores[current_index, latter_map[current_index][one] % len(nucleotides)] += score
            scores[current_index, latter_map[current_index][two] % len(nucleotides)] += score

        if repair_indel:
            delete_branch = [obtain_leaf_vertices(current_index, depth, latter_map=latter_map)]  # deletion
            for index in range(len(mutate_branches)):
                score = len(union1d(mutate_branches[index], delete_branch))
                scores[current_index, latter_map[current_index][index] % len(nucleotides)] += score
            del delete_branch

            for index, former_index in enumerate(latter_map[current_index]):  # insertion
                if former_index in latter_map:
                    for latter_index in latter_map[former_index]:
                        insert_branch = obtain_leaf_vertices(latter_index, depth, latter_map=latter_map)
                        score = len(union1d(mutate_branches[index], insert_branch))
                        scores[current_index, latter_map[current_index][index] % len(nucleotides)] += score

        monitor.output(current + 1, len(currents))

    return scores


def analyze_repair_match(task_seed, repeats, vertex_number):
    if not path.exists("./results/data/step_4_repairability_match.npy"):
        filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08"]
        whole_data, observed_length = [], 10
        random.seed(task_seed)
        for filter_index in filter_indices:
            if not path.exists("./results/temp/match" + filter_index + ".npy"):
                print("Analyze the generated algorithm under constraint set " + filter_index)
                accessor = load(file="./results/data/a" + filter_index + "[g].npy")
                vertices = where(load(file="./results/data/a" + filter_index + "[v].npy") == 1)[0]
                random.shuffle(vertices)
                chosen_vertices = vertices[:vertex_number]
                info = evaluate_repair_mechanism(accessor=accessor, vertices=chosen_vertices,
                                                 wrong_location=observed_length + 5, repeats=repeats)

                records = concatenate((ones(shape=(1, len(info))).T * int(filter_index), info), axis=1).astype(int)
                save(file="./results/temp/match" + filter_index + ".npy", arr=records)
            else:
                records = load(file="./results/temp/match" + filter_index + ".npy")

            whole_data += records.tolist()

        save(file="./results/data/step_4_repairability_match.npy", arr=array(whole_data))

        for filter_index in filter_indices:
            remove(path="./results/temp/match" + filter_index + ".npy")


def screen_edges_for_repair(evaluation_repeats, vertex_count):
    if not path.exists("./results/data/step_4_repairability_screen.pkl"):
        accessor = load(file="./results/data/a01[g].npy")
        print("Obtain latter map from accessor with the constraint set 1.")
        latter_map = accessor_to_latter_map(accessor=accessor, verbose=True)
        nucleotides, observed_length, repair_indel = ["A", "C", "G", "T"], 10, True
        print()

        print("Screen directed edges based on the maximum score.")
        round_number, records, depth = 1, [], observed_length - 1
        while True:  # try to screen for the determine (rather then probabilistic) repair model.
            scores = calculate_maximum_score(latter_map=latter_map, round_number=round_number, nucleotides=nucleotides,
                                             observed_length=observed_length, repair_indel=repair_indel)

            remove_edge, maximum_score, out = None, 0, len(nucleotides)
            locations = where(scores > 0)
            for vertex_index, latter_value in zip(locations[0], locations[1]):
                out_value, flag = len(latter_map[vertex_index]), False
                if maximum_score < scores[vertex_index, latter_value]:
                    flag = True
                if maximum_score == scores[vertex_index, latter_value] and out_value < out:
                    flag = True

                if flag:
                    out = out_value
                    current = vertex_index
                    latter = int((current * len(nucleotides) + latter_value) % (len(nucleotides) ** observed_length))
                    remove_edge = (current, latter)
                    maximum_score = scores[vertex_index, latter_value]

            scores = scores.reshape(-1)
            scores = scores[scores > 0].tolist()
            counter, score_record, value = Counter(scores), [], maximum_score
            while value > 0:
                if value in counter:
                    score_record.append([value, counter[value]])
                value -= 1

            score_record = array(score_record).T

            performance_record = calculate_performance(latter_map=latter_map, observed_length=observed_length,
                                                       evaluation_repeats=evaluation_repeats, vertex_count=vertex_count,
                                                       nucleotides=nucleotides)

            print()
            if performance_record is None:
                break

            capacity, counts, repair_record = performance_record

            with open("./results/temp/screen" + str(round_number - 1).zfill(4) + ".pkl", "wb") as file:
                psave((capacity, counts, score_record, repair_record), file)

            if remove_edge is None:
                break
            else:
                former, latter = remove_edge
                location = latter_map[former].index(latter)
                del latter_map[former][location]
                if len(latter_map[former]) == 0:
                    del latter_map[former]

                former_dna_string = number_to_dna(decimal_number=int(former), dna_length=10, nucleotides=nucleotides)
                latter_dna_string = number_to_dna(decimal_number=int(latter), dna_length=10, nucleotides=nucleotides)
                edge_info = former_dna_string + " -> " + latter_dna_string
                print("Remove the directed edge: " + edge_info + " with the score " + str(maximum_score) + ".")
                print("Remain scores are: \n" + str(score_record) + ".")
                print("The out-degree information of this graph is as follow:")
                print("Out degree | " + " " * 5 + "1" + " " * 5 + "2" + " " * 5 + "3" + " " * 5 + "4")
                count_info = ""
                for count in counts:
                    count_info += " " * (6 - len(str(count))) + str(count)
                print("Count      | " + count_info)

                previous_count = len(latter_map)
                latter_map = remove_useless(latter_map, threshold=1, verbose=False)
                current_count = len(latter_map)

                transcode_flag = False
                print("Remove " + str(previous_count - current_count) + " useless vertices that have no out-degree.")
                for current_index, latter_indices in latter_map.items():
                    if len(latter_indices) >= 2:
                        transcode_flag = True
                if not transcode_flag:
                    print("There are no nodes available for transcoding, exit!")
                    break

            round_number += 1

        whole_data = []
        for data_index in range(round_number):
            with open("./results/temp/screen" + str(data_index).zfill(4) + ".pkl", "rb") as file:
                whole_data.append(pload(file))

        with open("./results/data/step_4_repairability_screen.pkl", "wb") as file:
            psave(whole_data, file)

        for data_index in range(round_number):
            remove(path="./results/temp/screen" + str(data_index).zfill(4) + ".pkl")


def draw_average_length(heap_preset_size, reference_length_1, reference_length_2):
    error_rates = linspace(0, 0.1, 101)

    figure = pyplot.figure(figsize=(10, 5))
    pyplot.subplots_adjust(wspace=0.3)
    subgraph_1 = pyplot.subplot(1, 2, 1)

    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08"]
    data = load(file="./results/data/step_4_repairability_match.npy")
    strategy_numbers = [[] for _ in range(len(filter_indices))]
    for sample in data:
        strategy_number = sum(sample[4:7])
        if strategy_number == 0:
            strategy_number = 1
        strategy_numbers[sample[0] - 1].append(strategy_number)

    original_heap_sizes = array([mean(out_degrees) for out_degrees in strategy_numbers])

    with open("./results/data/step_1_compatibility_capacities.pkl", "rb") as file:
        original_capacities, _ = pload(file)

    original_matrix = []
    for error_rate in error_rates:
        average_lengths = log(heap_preset_size) / log(original_heap_sizes)
        if error_rate > 0:
            values = log10(average_lengths / error_rate)
        else:
            values = ones_like(average_lengths) * 5
        original_matrix.append(values - 1)

    original_matrix = array(original_matrix).T

    locations_1, locations_2 = [], []
    for values in original_matrix:
        collect_value = None
        for location, value in enumerate(values):
            if value < log10(float(reference_length_1)) - 1:
                collect_value = location
                break
        if collect_value is None:
            collect_value = len(values)

        locations_1.append(collect_value / 1000.0)

    locations_1 = array(locations_1)
    locations_1 = locations_1[locations_1 > 0]

    for values in original_matrix:
        collect_value = None
        for location, value in enumerate(values):
            if value < log10(float(reference_length_2)) - 1:
                collect_value = location
                break
        if collect_value is None:
            collect_value = len(values)

        locations_2.append(collect_value / 1000.0)

    locations_2 = array(locations_2)
    locations_2 = locations_2[locations_2 > 0]

    info_1 = pyplot.pcolormesh(error_rates, range(9), original_matrix,
                               cmap=pyplot.get_cmap(name="RdYlGn"), vmin=0, vmax=4)

    pyplot.scatter(locations_1, linspace(1, 8, 8) - 0.5, s=25, marker="^",
                   color=colors["trad2"], edgecolor="black", label="200bp", zorder=2)
    pyplot.scatter(locations_2, linspace(1, 8, 8) - 0.5, s=25, marker="p",
                   color=colors["trad2"], edgecolor="black", label="1000bp", zorder=2)

    pyplot.legend(loc="upper right", fontsize=8)
    pyplot.xlabel("total DNA error rate", fontsize=8)
    pyplot.xlim(0, 0.1)
    pyplot.xticks([0, 0.02, 0.04, 0.06, 0.08, 0.1], ["0%", "2%", "4%", "6%", "8%", "10%"], fontsize=8)

    pyplot.ylabel("constraint set (capacity)", fontsize=8)
    y_tick_labels = [filter_indices[index] + "(%.4f)" % original_capacities[index]
                     for index in range(len(filter_indices))]
    pyplot.yticks(linspace(1, 8, 8) - 0.5, y_tick_labels, fontsize=8)
    pyplot.ylim(0, 8)

    subgraph_2 = pyplot.subplot(1, 2, 2)

    with open("./results/data/step_4_repairability_screen.pkl", "rb") as file:
        screened_data_group = pload(file)

    screened_capacities, heap_sizes = [], []
    for capacity, counts, score_record, repair_record in screened_data_group:
        screened_capacities.append(capacity)
        repair_data = sum(repair_record[:, 3:6], axis=1)
        repair_data[repair_data == 0] = 1
        heap_sizes.append(mean(repair_data))

    screened_capacities = array(screened_capacities)
    heap_sizes = array(heap_sizes)

    screened_capacity_labels = linspace(0.2, 1.0, 81)
    temp = [[] for _ in range(len(screened_capacity_labels))]
    for capacity, level in zip(screened_capacities, log(heap_preset_size) / log(heap_sizes)):
        temp[int(round(capacity, 2) * 100) - 20].append(level)

    screened_capacity_levels = [mean(values) if len(values) > 0 else -1 for values in temp]
    for index in range(len(screened_capacity_levels)):
        if screened_capacity_levels[index] == -1:
            estimated_value = (screened_capacity_levels[index - 1] + screened_capacity_levels[index + 1]) / 2.0
            screened_capacity_levels[index] = estimated_value
    screened_capacity_levels = array(screened_capacity_levels)

    matrix = []
    for error_rate in error_rates:
        if error_rate > 0:
            values = log10(screened_capacity_levels / error_rate)
        else:
            values = ones_like(screened_capacity_levels) * 5
        matrix.append(values - 1)

    matrix = array(matrix).T

    thresholds_1 = []
    for values in matrix.T:
        collect_value = None
        for location, value in enumerate(values):
            if value < log10(float(reference_length_1)) - 1:
                collect_value = (location + 20 - 1) / 100.0
                break
        if collect_value is None:
            collect_value = (len(values) + 20 - 1) / 100.0

        thresholds_1.append(collect_value)

    thresholds_1 = array(thresholds_1)

    thresholds_2 = []
    for values in matrix.T:
        collect_value = None
        for location, value in enumerate(values):
            if value < log10(float(reference_length_2)) - 1:
                collect_value = (location + 20 - 1) / 100.0
                break
        if collect_value is None:
            collect_value = (len(values) + 20 - 1) / 100.0

        thresholds_2.append(collect_value)

    pyplot.pcolormesh(error_rates, screened_capacity_labels, matrix,
                      cmap=pyplot.get_cmap(name="RdYlGn"), vmin=0, vmax=4)

    ignore = len(thresholds_1[thresholds_1 == 1.0]) - 1
    pyplot.plot(error_rates[ignore: len(thresholds_1)], thresholds_1[ignore:],
                color="black", linestyle=":", linewidth=1, label="200bp", zorder=2)
    pyplot.plot(error_rates[:len(thresholds_2)], thresholds_2,
                color="black", linewidth=1, label="1000bp", zorder=2)

    pyplot.legend(loc="upper right", fontsize=8)

    pyplot.xlabel("total DNA error rate", fontsize=8)
    pyplot.xticks([0, 0.02, 0.04, 0.06, 0.08, 0.1], ["0%", "2%", "4%", "6%", "8%", "10%"], fontsize=8)
    pyplot.xlim(0, 0.1)
    pyplot.ylabel("screened capacity from constraint set 01", fontsize=8)
    pyplot.yticks([0.2, 0.4, 0.6, 0.8, 1], ["0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=8)
    pyplot.ylim(0.2, 1)

    bar = figure.colorbar(info_1, ax=[subgraph_1, subgraph_2])
    bar.ax.set_title("average DNA length", fontsize=8)
    bar.set_ticks([0, 1, 2, 3, 4])
    bar.set_ticklabels([r"$\leq$" + "1e1", r"$=$" + "1e2", r"$=$" + "1e3", r"$=$" + "1e4", r"$\geq$" + "1e5"])
    bar.update_ticks()
    bar.ax.tick_params(labelsize=8)

    figure.text(0.043, 0.9, "A", va="center", ha="center")
    figure.text(0.436, 0.9, "B", va="center", ha="center")
    pyplot.savefig("./results/figures/[4-1] repairability totally.svg", format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


def draw_screen_record(heap_preset_size, error_rate):
    with open("./results/data/step_4_repairability_screen.pkl", "rb") as file:
        screened_data_group = pload(file)

    capacity_group, out_group, length_group, score_group, data_number = [], [], [], [], 0
    for capacity, counts, score_record, repair_record in screened_data_group:
        capacity_group.append(capacity)

        repair_data = sum(repair_record[:, 3:6], axis=1)
        repair_data[repair_data == 0] = 1
        heap_size = mean(repair_data)
        length_group.append(log(heap_preset_size) / log(float(heap_size)) / error_rate)

        out_group.append(counts.tolist())

        level_values = zeros(shape=(41,))
        for score, count in zip(score_record[0][::-1], score_record[1][::-1]):
            level = int(score / 4096.0 * 40 + 0.5)
            level_values[level] += count
        score_group.append(level_values.tolist())

        data_number += 1

    iterations = linspace(0, data_number - 1, data_number)
    bottoms = zeros(shape=(data_number,))

    capacity_group, out_group = array(capacity_group), array(out_group).T
    length_group, score_group = array(length_group), array(score_group).T

    figure = pyplot.figure(figsize=(10, 5), tight_layout=True)
    pyplot.subplots_adjust(hspace=0.3, wspace=0.3)

    pyplot.subplot(2, 2, 1)
    pyplot.fill_between(iterations, bottoms, capacity_group, color=colors["diffs"], linewidth=0, zorder=1)
    pyplot.plot(iterations, capacity_group, color="silver", linewidth=1, zorder=2)

    pyplot.xticks([0, 393, 786, 1179, 1572], ["", "", "", "", ""], fontsize=8)
    pyplot.yticks([0, 0.3, 0.6, 0.9, 1.2], ["0.0", "0.3", "0.6", "0.9", "1.2"], fontsize=8)
    pyplot.ylabel("screened capacity", fontsize=8)
    pyplot.xlim(0, 1572)
    pyplot.ylim(0, 1.2)

    pyplot.subplot(2, 2, 2)
    pyplot.fill_between(iterations, bottoms, length_group, color=colors["diffs"], linewidth=0, zorder=1)
    pyplot.plot(list(range(len(length_group))), length_group, color="silver", linewidth=1, zorder=2)

    pyplot.xticks([0, 393, 786, 1179, 1572], ["", "", "", "", ""], fontsize=8)
    pyplot.yticks([0, 150, 300, 450, 600], [0, 150, 300, 450, 600], fontsize=8)
    pyplot.ylabel("average DNA length", fontsize=8)
    pyplot.xlim(0, 1572)
    pyplot.ylim(0, 600)

    pyplot.subplot(2, 2, 3)
    gradient_colors = pyplot.get_cmap(name="Blues")(linspace(0, 1, 5))
    area_info = cumsum(out_group, axis=0)[::-1]
    for index, area in enumerate(area_info):
        pyplot.fill_between(iterations, bottoms, area, color=gradient_colors[4 - index], linewidth=0)

    pyplot.fill_between([1250, 1550], [200 - 140, 200 - 140], [1450 - 140, 1450 - 140],
                        color="white", alpha=0.7)
    pyplot.vlines(1250, 200 - 140, 1450 - 140, color=colors["diffs"], linewidth=0.75, zorder=5)
    pyplot.vlines(1550, 200 - 140, 1450 - 140, color=colors["diffs"], linewidth=0.75, zorder=5)
    pyplot.hlines(200 - 140, 1250, 1550, color=colors["diffs"], linewidth=0.75, zorder=5)
    pyplot.hlines(1450 - 140, 1250, 1550, color=colors["diffs"], linewidth=0.75, zorder=5)
    pyplot.text(x=1400, y=1250 - 140, s="out-degree", va="center", ha="center", fontsize=8, zorder=5)
    pyplot.fill_between([1350, 1400], [300 - 140, 300 - 140], [450 - 140, 450 - 140],
                        color=gradient_colors[0], zorder=5)
    pyplot.fill_between([1350, 1400], [450 - 140, 450 - 140], [600 - 140, 600 - 140],
                        color=gradient_colors[1], zorder=5)
    pyplot.fill_between([1350, 1400], [600 - 140, 600 - 140], [750 - 140, 750 - 140],
                        color=gradient_colors[2], zorder=5)
    pyplot.fill_between([1350, 1400], [750 - 140, 750 - 140], [900 - 140, 900 - 140],
                        color=gradient_colors[3], zorder=5)
    pyplot.fill_between([1350, 1400], [900 - 140, 900 - 140], [1050 - 140, 1050 - 140],
                        color=gradient_colors[3], zorder=5)
    pyplot.vlines(1350, 300 - 140, 1050 - 140, color="black", linewidth=0.75, zorder=6)
    pyplot.vlines(1400, 300 - 140, 1050 - 140, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines(300 - 140, 1350, 1400, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines(1050 - 140, 1350, 1400, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((300 + 450) / 2 - 140, 1400, 1420, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((450 + 600) / 2 - 140, 1400, 1420, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((600 + 750) / 2 - 140, 1400, 1420, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((750 + 900) / 2 - 140, 1400, 1420, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((900 + 1050) / 2 - 140, 1400, 1420, color="black", linewidth=0.75, zorder=6)
    pyplot.text(1440, (300 + 450) / 2 - 10 - 140, 0, va="center", ha="center", fontsize=8, zorder=6)
    pyplot.text(1440, (450 + 600) / 2 - 10 - 140, 1, va="center", ha="center", fontsize=8, zorder=6)
    pyplot.text(1440, (600 + 750) / 2 - 10 - 140, 2, va="center", ha="center", fontsize=8, zorder=6)
    pyplot.text(1440, (750 + 900) / 2 - 10 - 140, 3, va="center", ha="center", fontsize=8, zorder=6)
    pyplot.text(1440, (900 + 1050) / 2 - 10 - 140, 4, va="center", ha="center", fontsize=8, zorder=6)

    pyplot.xlabel("screening iteration based on maximum overlap score", fontsize=8)
    pyplot.xticks([0, 393, 786, 1179, 1572], [0, 393, 786, 1179, 1572], fontsize=8)
    pyplot.yticks([0, 600, 1200, 1800, 2400], [0, 600, 1200, 1800, 2400], fontsize=8)
    pyplot.ylabel("vertex number", fontsize=8)
    pyplot.xlim(0, 1572)
    pyplot.ylim(0, 2400)

    pyplot.subplot(2, 2, 4)
    gradient_colors = pyplot.get_cmap(name="YlOrRd")(linspace(0, 1, 101))
    area_info = cumsum(score_group, axis=0)[::-1]
    for index, area in enumerate(area_info):
        pyplot.fill_between(iterations, bottoms, area, color=gradient_colors[100 - index], linewidth=0)

    pyplot.fill_between([1250, 1550], [(200 - 140) * 2, (200 - 140) * 2], [(1450 - 140) * 2, (1450 - 140) * 2],
                        color="white", alpha=0.7)
    pyplot.vlines(1250, (200 - 140) * 2, (1450 - 140) * 2, color=colors["diffs"], linewidth=0.75, zorder=5)
    pyplot.vlines(1550, (200 - 140) * 2, (1450 - 140) * 2, color=colors["diffs"], linewidth=0.75, zorder=5)
    pyplot.hlines((200 - 140) * 2, 1250, 1550, color=colors["diffs"], linewidth=0.75, zorder=5)
    pyplot.hlines((1450 - 140) * 2, 1250, 1550, color=colors["diffs"], linewidth=0.75, zorder=5)
    pyplot.text(x=1400, y=(1250 - 140) * 2, s="score", va="center", ha="center", fontsize=8, zorder=5)
    interval = (1050 - 300) / len(gradient_colors) * 2
    for index in range(len(gradient_colors)):
        lower, upper = 300 + interval * index, 300 + interval * (index + 1)
        pyplot.fill_between([1310, 1360], [lower, lower], [upper, upper], color=gradient_colors[index], zorder=5)
    pyplot.vlines(1310, (300 - 140) * 2, (1050 - 140) * 2, color="black", linewidth=0.75, zorder=6)
    pyplot.vlines(1360, (300 - 140) * 2, (1050 - 140) * 2, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((300 - 140) * 2, 1310, 1360, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((1050 - 140) * 2, 1310, 1360, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((300 - 140) * 2, 1360, 1380, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((487.5 - 140) * 2, 1360, 1380, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((675 - 140) * 2, 1360, 1380, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((862.5 - 140) * 2, 1360, 1380, color="black", linewidth=0.75, zorder=6)
    pyplot.hlines((1050 - 140) * 2, 1360, 1380, color="black", linewidth=0.75, zorder=6)
    pyplot.text(1405, (300 - 140) * 2 - 15, 0, va="center", ha="center", fontsize=8, zorder=6)
    pyplot.text(1445, (487.5 - 140) * 2 - 15, 1024, va="center", ha="center", fontsize=8, zorder=6)
    pyplot.text(1443, (675 - 140) * 2 - 15, 2048, va="center", ha="center", fontsize=8, zorder=6)
    pyplot.text(1445, (862.5 - 140) * 2 - 15, 3072, va="center", ha="center", fontsize=8, zorder=6)
    pyplot.text(1445, (1050 - 140) * 2 - 15, 4096, va="center", ha="center", fontsize=8, zorder=6)

    pyplot.xlabel("screening iteration based on maximum overlap score", fontsize=8)
    pyplot.ylabel("edge number", fontsize=8)
    pyplot.xticks([0, 393, 786, 1179, 1572], [0, 393, 786, 1179, 1572], fontsize=8)
    pyplot.yticks([0, 1200, 2400, 3600, 4800], [0, 1200, 2400, 3600, 4800], fontsize=8)
    pyplot.xlim(0, 1572)
    pyplot.ylim(0, 4800)

    figure.align_labels()

    figure.text(0.020, 0.97, "A", va="center", ha="center")
    figure.text(0.513, 0.97, "B", va="center", ha="center")
    figure.text(0.020, 0.52, "C", va="center", ha="center")
    figure.text(0.513, 0.52, "D", va="center", ha="center")
    pyplot.savefig("./results/figures/[4-2] repairability record.svg", format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


def draw_heap_situation():
    figure = pyplot.figure(figsize=(10, 5), tight_layout=True)

    pyplot.subplot(1, 2, 1)
    filter_indices = ["01", "02", "03", "04", "05", "06", "07", "08"]
    data = load(file="./results/data/step_4_repairability_match.npy")
    strategy_numbers = [[] for _ in range(len(filter_indices))]
    for sample in data:
        strategy_number = sum(sample[4:7])
        if strategy_number == 0:
            strategy_number = 1
        strategy_numbers[sample[0] - 1].append(strategy_number)

    with open("./results/data/step_1_compatibility_capacities.pkl", "rb") as file:
        original_capacities, _ = pload(file)

    violin = pyplot.violinplot(strategy_numbers, bw_method=0.5, showextrema=False)
    for patch in violin["bodies"]:
        patch.set_edgecolor(colors["algo1"])
        patch.set_facecolor(colors["algo2"])
        patch.set_linewidth(1)
        patch.set_alpha(1)

    averages = [mean(values) for values in strategy_numbers]
    pyplot.scatter(linspace(1, 8, 8), averages, color="white", edgecolor=colors["algo1"], linewidth=1, s=15, zorder=3)

    pyplot.xlabel("original capacity (constraint set)", fontsize=8)
    pyplot.xlim(0.5, 8.5)
    x_tick_labels = ["%.4f\n" % original_capacities[index] + "(" + filter_indices[index] + ")"
                     for index in range(len(filter_indices))]
    pyplot.xticks(linspace(1, 8, 8), x_tick_labels, fontsize=8)
    pyplot.ylabel("number of repair strategy", fontsize=8)
    pyplot.ylim(0, 64)
    pyplot.yticks([0, 8, 16, 24, 32, 40, 48, 56, 64], [0, 8, 16, 24, 32, 40, 48, 56, 64], fontsize=8)

    pyplot.subplot(1, 2, 2)

    range_labels = ["[0.2,0.3)", "[0.3,0.4)", "[0.4,0.5)", "[0.5,0.6)",
                    "[0.6,0.7)", "[0.7,0.8)", "[0.8,0.9)", "[0.9,1.0)"]

    with open("./results/data/step_4_repairability_screen.pkl", "rb") as file:
        screened_data_group = pload(file)

    strategy_numbers = [[] for _ in range(len(range_labels))]
    for capacity, _, _, repair_record in screened_data_group:
        if capacity < 1:
            repair_data = sum(repair_record[:, 3:6], axis=1)
            repair_data[repair_data == 0] = 1
            strategy_numbers[int((capacity - 0.2) * 10)] += repair_data.tolist()

    numbers = [len(values) for values in strategy_numbers]

    for index in range(len(strategy_numbers)):
        violin = pyplot.violinplot(strategy_numbers[index], positions=[index + 1], bw_method=0.5, showextrema=False)
        for patch in violin["bodies"]:
            patch.set_edgecolor(colors["algo1"])
            patch.set_facecolor(colors["algo2"])
            patch.set_linewidth(1)
            patch.set_alpha(1)

    averages = [mean(values) for values in strategy_numbers]
    pyplot.scatter(linspace(1, 8, 8), averages, color="white", edgecolor=colors["algo1"], linewidth=1, s=15, zorder=3)

    pyplot.xlabel("screened capacity range (sample number)", fontsize=8)
    pyplot.xlim(0.5, 8.5)
    x_tick_labels = [range_label + "\n(" + ("%.0e" % number).replace("+0", "") + ")"
                     for range_label, number in zip(range_labels, numbers)]
    pyplot.xticks(linspace(1, 8, 8), x_tick_labels, fontsize=8)
    pyplot.ylabel("  ", fontsize=8)
    pyplot.ylim(0, 16)
    pyplot.yticks([0, 2, 4, 6, 8, 10, 12, 14, 16], [0, 2, 4, 6, 8, 10, 12, 14, 16], fontsize=8)

    figure.text(0.02, 0.98, "A", va="center", ha="center")
    figure.text(0.51, 0.98, "B", va="center", ha="center")

    pyplot.savefig("./results/figures/[4-3] repairability heaps.svg", format="svg", bbox_inches="tight", dpi=600)
    pyplot.close()


if __name__ == "__main__":
    create_folders()

    repair_examples()

    analyze_repair_match(task_seed=2021, repeats=100, vertex_number=100)
    screen_edges_for_repair(evaluation_repeats=50, vertex_count=40)

    draw_average_length(heap_preset_size=1e6, reference_length_1=200, reference_length_2=1000)
    draw_screen_record(heap_preset_size=1e6, error_rate=0.1)
    draw_heap_situation()
