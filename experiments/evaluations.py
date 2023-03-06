from hashlib import md5
from itertools import permutations
from numpy import random, array, zeros, abs, sum, min, max, linalg, real, log2, log10, ceil, where
from os.path import exists, getsize

from dsw import Monitor, obtain_vertices, find_vertices, connect_valid_graph, connect_coding_graph
from dsw import encode, bit_to_number, calculus_division, approximate_capacity
from dsw import get_complete_accessor, accessor_to_adjacency_matrix, adjacency_matrix_to_accessor

from experiments import local_bio_filters, special_filter, load_data, save_data
from experiments.code_encode import generate_database, conversion_step
from experiments.code_encode import trans_fountain, trans_yinyang, trans_hedges, trans_spiderweb
from experiments.code_repair import evaluate_repair_errors, calculate_minimum_reads, evaluate_established_depth
from experiments.code_repair import spiderweb_step, hedges_step


def evaluate_performance():
    if (not exists(path="./raw/graph_primal.pkl")) or (not exists(path="./raw/graph_coding.pkl")):
        # parameters during the generation process are recorded by artificially.
        primal_graphs, coding_graphs = {}, {}
        for index, bio_filter in local_bio_filters.items():
            print("Calculate graph " + str(index) + ".")
            vertices = find_vertices(observed_length=10, bio_filter=bio_filter, verbose=True)
            print()
            primal_accessor = connect_valid_graph(observed_length=10, vertices=vertices, verbose=True)
            print()
            _, coding_accessor = connect_coding_graph(observed_length=10, vertices=vertices, threshold=2, verbose=True)
            print()
            primal_graphs[index], coding_graphs[index] = primal_accessor, coding_accessor
        save_data(save_path="./raw/graph_primal.pkl", information=primal_graphs)
        save_data(save_path="./raw/graph_coding.pkl", information=coding_graphs)

    if not exists(path="./raw/capacity_reference.pkl"):
        primal_graphs = load_data(load_path="./raw/graph_primal.pkl")
        capacities, data = [], []
        for index in range(12):
            accessor = primal_graphs[str(index + 1).zfill(2)]
            capacity, processes = approximate_capacity(accessor=accessor, repeats=10, process=True)
            capacities.append(capacity)
            data.append((capacity, processes))
        save_data(save_path="./raw/capacity_reference.pkl", information=(capacities, data))

    if not exists(path="./raw/coding_evaluation_1.pkl"):
        dataset, records = generate_database(random_seed=2021, total_length=8 * 1024 * 64, times=100), []
        records += trans_fountain(dataset, shape=(2048, 256), index_length=32, random_seed=2021, param_number=100)
        records += trans_yinyang(dataset, shape=(4096, 128), index_length=16, random_seed=2021, param_number=100)
        records += trans_hedges(dataset, shape=(4096, 128), index_length=16)
        records += trans_spiderweb(dataset, shape=(4096, 128), index_length=16, random_seed=2021, param_number=100)
        save_data(save_path="./raw/coding_evaluation_1.pkl", information=array(records))

    if not exists(path="./raw/coding_evaluation_2.pkl"):
        monitor, task_seed, bit_length, repeats = Monitor(), 2021, 200, 100
        coding_graphs = load_data(load_path="./raw/graph_coding.pkl")
        results_1, results_2 = [], {}
        for index in range(12):
            filter_index = str(index + 1).zfill(2)
            print("Evaluate coding algorithm " + filter_index + ".")
            accessor = coding_graphs[filter_index]
            vertices = obtain_vertices(accessor=accessor)
            random.seed(task_seed)
            random.shuffle(vertices)
            encoded_matrix, nucleotide_numbers = random.randint(0, 2, size=(repeats, bit_length)), []
            for current, start_index in enumerate(vertices[:repeats]):
                for binary_message in encoded_matrix:
                    dna_string = encode(binary_message=binary_message, accessor=accessor, start_index=start_index)
                    nucleotide_numbers.append(len(dna_string))
                monitor(current + 1, repeats)
            results_1.append(bit_length / array(nucleotide_numbers))
        primal_graphs = load_data(load_path="./raw/graph_primal.pkl")
        for filter_index, bio_filter in local_bio_filters.items():
            accessor = primal_graphs[filter_index]
            out_degrees = sum(where(accessor >= 0, 1, 0), axis=1)
            if len(out_degrees[out_degrees == 1]) == 0:
                continue
            print("Evaluate coding algorithm " + filter_index + " [out-degree 1].")
            random.seed(task_seed)
            encoded_matrix, nucleotide_numbers, current = random.randint(0, 2, size=(repeats, bit_length)), [], 0
            for mapping in permutations(["A", "C", "G", "T"]):
                for binary_message in encoded_matrix:
                    available, quotient, dna_string, monitor = mapping, bit_to_number(binary_message), "", Monitor()
                    while quotient != "0":
                        if len(available) > 1:  # current vertex contains information.
                            quotient, remainder = calculus_division(number=quotient, base=str(len(available)))
                            dna_string += available[int(remainder)]
                        elif len(available) == 1:  # current vertex does not contain information.
                            dna_string += available[0]
                        available = []
                        for potential_nucleotide in mapping:
                            if bio_filter.valid(dna_string + potential_nucleotide, only_last=True):
                                available.append(potential_nucleotide)
                        if len(available) == 0:
                            dna_string = ""
                            break
                    if len(dna_string) > 0:
                        nucleotide_numbers.append(len(dna_string))
                    else:
                        nucleotide_numbers.append(-1)

                    monitor(current + 1, len(encoded_matrix) * 24)
                    current += 1
            results_2[filter_index] = bit_length / array(nucleotide_numbers)
            save_data(save_path="./raw/coding_evaluation_2.pkl", information=(results_1, results_2))


def evaluate_correction():
    random_seed = 2021

    if not exists(path="./raw/special_graph.npy"):
        print("Generate the special coding digraph.")
        vertices = find_vertices(observed_length=10, bio_filter=special_filter, verbose=True)
        _, accessor = connect_coding_graph(observed_length=10, vertices=vertices, threshold=1, verbose=True)
        save_data(save_path="./raw/special_graph.npy", information=accessor)

    if not exists(path="./raw/correction_evaluation_1.pkl"):
        print("Evaluate correction rate of the special coding digraph.")
        records = {}
        accessor = load_data(load_path="./raw/special_graph.npy")
        vertices = obtain_vertices(accessor)
        for errors in range(1, 10, 1):
            results = evaluate_repair_errors(random_seed=random_seed, dna_length=200, repeats=10000, errors=errors,
                                             accessor=accessor, vertices=vertices, observed_length=10)

            records[errors / 200.0] = (results[0], results[1])
        save_data(save_path="./raw/correction_evaluation_1.pkl", information=records)

    if not exists(path="./raw/correction_evaluation_2.pkl"):
        print("Evaluate minimum pretreatment-free depth of the special coding digraph.")
        records = {}
        accessor = load_data(load_path="./raw/special_graph.npy")
        vertices = obtain_vertices(accessor)
        for errors in [1, 2, 3, 4, 5, 6, 7, 8]:
            results = calculate_minimum_reads(random_seed=random_seed, dna_length=200, repeats=10000, errors=errors,
                                              accessor=accessor, vertices=vertices, observed_length=10)
            records[errors / 200.0] = results
        save_data(save_path="./raw/correction_evaluation_2.pkl", information=records)

    if not exists(path="./raw/correction_evaluation_3.pkl"):
        print("Evaluate correction performance of the special coding digraph under 1E+6 sequence diversity.")
        accessor = load_data(load_path="./raw/special_graph.npy")
        vertices = obtain_vertices(accessor)
        repeats = 100

        records, sequence_diversity = {"loss": zeros(shape=(8, 10)), "gap": zeros(shape=(8, 10))}, 1000000
        for index_1, error_number in enumerate([1, 2, 3, 4, 5, 6, 7, 8]):
            for index_2, reads_number in enumerate([5, 10, 15, 20, 25, 30, 35, 40, 45, 50]):
                maximum_loss, smallest_gap = 0, 50
                for repeat in range(repeats):
                    results = evaluate_established_depth(random_seed=random_seed + repeat,
                                                         dna_length=200, observed_length=10,
                                                         accessor=accessor, vertices=vertices, errors=error_number,
                                                         classes=sequence_diversity, established_depth=reads_number)
                    right_strings, obtained_set = results
                    right_strings = set(right_strings)
                    final_results = sorted(obtained_set.items(), key=lambda sample: sample[1], reverse=True)
                    loss_count, maximum_incorrect, minimum_correct = 0, 0, 50

                    for rank_index, (obtained_string, count) in enumerate(final_results):
                        if obtained_string not in right_strings:
                            if rank_index < sequence_diversity:
                                loss_count += 1
                            maximum_incorrect = max([maximum_incorrect, count])
                        else:
                            minimum_correct = min([minimum_correct, count])

                    maximum_loss = max([loss_count, maximum_loss])
                    smallest_gap = min([smallest_gap, minimum_correct - maximum_incorrect])

                records["loss"][index_1, index_2] = maximum_loss
                records["gap"][index_1, index_2] = smallest_gap

            save_data(save_path="./raw/correction_evaluation_3.pkl", information=records)

    if not exists(path="./raw/correction_evaluation_4.pkl") or not exists(path="./raw/correction_evaluation_5.pkl"):
        print("Evaluate minimum reads number of the special coding digraph under different sequence diversities.")
        accessor = load_data(load_path="./raw/special_graph.npy")
        vertices = obtain_vertices(accessor)
        repeats = 100

        records = {}
        for error_number in [1, 8]:
            records[(error_number / 200, 0.990)] = zeros(shape=(6,))
            records[(error_number / 200, 0.999)] = zeros(shape=(6,))
            records[(error_number / 200, 1.000)] = zeros(shape=(6,))
            for index, sequence_diversity in [10, 100, 1000, 10000, 100000, 1000000]:
                reads_number, flags = 1, 0
                while True:
                    count_group, maximum_incorrect_group = [], []
                    for repeat in range(repeats):
                        results = evaluate_established_depth(random_seed=random_seed + repeat,
                                                             dna_length=200, observed_length=10,
                                                             accessor=accessor, vertices=vertices, errors=error_number,
                                                             classes=sequence_diversity, established_depth=reads_number)
                        right_strings, obtained_set = results
                        right_strings = set(right_strings)
                        final_results = sorted(obtained_set.items(), key=lambda sample: sample[1], reverse=True)
                        success_count = sequence_diversity

                        for rank_index, (obtained_string, count) in enumerate(final_results):
                            if obtained_string not in right_strings:
                                if rank_index < sequence_diversity:
                                    success_count -= 1

                        count_group.append(success_count)

                    count_group = array(count_group)

                    if min(count_group) >= sequence_diversity * 0.99 and flags == 0:
                        records[(error_number / 200, 0.990)][index] = reads_number
                        flags += 1

                    if min(count_group) >= sequence_diversity * 0.999 and flags == 1:
                        records[(error_number / 200, 0.999)][index] = reads_number
                        flags += 1

                    if min(count_group) == sequence_diversity and flags == 2:
                        records[(error_number / 200, 0.999)][index] = reads_number
                        break

                    reads_number += 1

            save_data(save_path="./raw/correction_evaluation_4.pkl", information=records)
            # correction_evaluation_5.pkl is the curve fitting data.
            # save_data(save_path="./raw/correction_evaluation_5.pkl", information=fitted_data)

    if not exists(path="./raw/correction_evaluation_6.pkl"):
        print("Evaluate non-blocking reads of the special coding digraph under different sequence diversities.")
        accessor = load_data(load_path="./raw/special_graph.npy")
        vertices = obtain_vertices(accessor)
        repeats = 100

        records = {0.005: zeros(shape=(5, 6)), 0.040: zeros(shape=(5, 6))}
        for error_number in [1, 8]:
            for reads_number in [10, 20, 30, 40, 50]:
                for sequence_diversity in [10, 100, 1000, 10000, 100000, 1000000]:
                    maximum_incorrect_group = []
                    for repeat in range(repeats):
                        results = evaluate_established_depth(random_seed=random_seed + repeat,
                                                             dna_length=200, observed_length=10,
                                                             accessor=accessor, vertices=vertices, errors=error_number,
                                                             classes=sequence_diversity, established_depth=reads_number)
                        right_strings, obtained_set = results
                        right_strings = set(right_strings)
                        final_results = sorted(obtained_set.items(), key=lambda sample: sample[1], reverse=True)

                        maximum_incorrect = 0
                        for obtained_string, count in final_results:
                            if obtained_string not in right_strings:
                                maximum_incorrect = max([maximum_incorrect, count])
                                break
                        maximum_incorrect_group.append(maximum_incorrect)

                    maximum_incorrect_group = array(maximum_incorrect_group)
                    index_1, index_2 = reads_number // 10 - 1, int(log10(sequence_diversity)) - 1
                    records[error_number / 200.0][index_1, index_2] = max(maximum_incorrect_group)

            save_data(save_path="./raw/correction_evaluation_6.pkl", information=records)

    if not exists(path="./raw/correction_evaluation_7.pkl"):
        records = {"spiderweb": {}, "hedges": {}}
        accessor = load_data(load_path="./raw/special_graph.npy")
        vertices = obtain_vertices(accessor)
        repeats = 1000

        for errors in range(9):
            results = spiderweb_step(random_seed=random_seed, accessor=accessor, vertices=vertices,
                                     dna_length=200, observed_length=10, repeats=repeats, errors=errors)
            records["spiderweb"][errors / 200.0] = results
        for errors in range(9):
            results = hedges_step(random_seed=random_seed, bio_filter=special_filter,
                                  dna_length=200, observed_length=10, repeats=repeats, errors=errors)
            records["hedges"][errors / 200.0] = results
        save_data(save_path="./raw/correction_evaluation_7.pkl", information=records)


def evaluate_reliability():
    if not exists("./raw/capacity_evaluation.pkl"):
        test_times, terminal_length, random_seed, records = 100, 6, 2021, {}

        random.seed(random_seed)
        feed_matrix = accessor_to_adjacency_matrix(get_complete_accessor(observed_length=2))
        data, errors, monitor = [], [], Monitor()
        for time in range(test_times):
            matrix = feed_matrix.copy()
            for position in array(list(where(feed_matrix == 1)[:2])).T:
                if random.random(1)[0] <= 0.5:
                    matrix[position[0], position[1]] = 0
            monitor(current_state=time + 1, total_state=test_times)
            calculated_eigenvalue = real(max(linalg.eig(matrix)[0]))
            calculated_value = log2(calculated_eigenvalue) if calculated_eigenvalue > 0.0 else 0.0
            approximate_value = approximate_capacity(accessor=adjacency_matrix_to_accessor(matrix), repeats=10)
            data.append(approximate_value)
            errors.append(abs(approximate_value - calculated_value))

        records["detail"] = (data, errors)

        whole_data = {}
        for length in range(2, terminal_length + 1):
            print("calculate length " + str(length))
            feed_matrix = accessor_to_adjacency_matrix(get_complete_accessor(observed_length=length))
            values, errors, monitor = [], [], Monitor()
            for time in range(test_times):
                matrix = feed_matrix.copy()
                for position in array(list(where(feed_matrix == 1)[:2])).T:
                    if random.random(1)[0] <= 0.5:
                        matrix[position[0], position[1]] = 0
                monitor(current_state=time + 1, total_state=test_times)
                calculated_eigenvalue = real(max(linalg.eig(matrix)[0]))
                calculated_value = log2(calculated_eigenvalue) if calculated_eigenvalue > 0.0 else 0.0
                approximate_value = approximate_capacity(accessor=adjacency_matrix_to_accessor(matrix), repeats=10)
                values.append(approximate_value)
                errors.append(abs(approximate_value - calculated_value))
            whole_data[length] = (values, errors)

        records["extend"] = whole_data

        save_data(save_path="./raw/capacity_evaluation.pkl", information=records)


def evaluate_conversion():
    if not exists("./raw/conversion_evaluation.pkl"):
        accessor = load_data(load_path="./raw/special_graph.npy")
        vertices = obtain_vertices(accessor)
        repeats = 100

        record = {}
        for index, start_vertex in enumerate(vertices):
            results = conversion_step(accessor=accessor, dna_length=200, start_index=start_vertex,
                                      repeats=repeats, random_seed=2021)
            record[start_vertex] = results

        save_data(save_path="./raw/conversion_evaluation.pkl", information=record)


if __name__ == "__main__":
    evaluate_performance()
    evaluate_correction()
    evaluate_reliability()
    evaluate_conversion()

    child_paths = [
        "generation_evaluation.pkl",
        "special_graph.npy",
        "graph_primal.pkl",
        "graph_coding.pkl",
        "coding_evaluation_1.pkl",
        "coding_evaluation_2.pkl",
        "conversion_evaluation.pkl",
        "correction_evaluation_1.pkl",
        "correction_evaluation_2.pkl",
        "correction_evaluation_3.pkl",
        "correction_evaluation_4.pkl",
        "correction_evaluation_5.pkl",
        "correction_evaluation_6.pkl",
        "correction_evaluation_7.pkl",
        "capacity_reference.pkl",
        "capacity_evaluation.pkl"
    ]
    for child_path in child_paths:
        if "pkl" in child_path or "npy" in child_path:
            md5_hash = md5()
            with open("./raw/" + child_path, "rb") as f:
                for byte_block in iter(lambda: f.read(4096), b""):
                    md5_hash.update(byte_block)
            file_size = getsize("./raw/" + child_path)
            print(child_path, (md5_hash.hexdigest()).upper(), ceil(file_size / 1024).astype(int))
