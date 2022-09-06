from itertools import permutations
from numpy import random, array, arange, zeros, abs, sum, min, max, linalg, real, log2, where, any, save
from os.path import exists
from pickle import load, dump

from dsw import Monitor, obtain_vertices, find_vertices, connect_valid_graph, connect_coding_graph
from dsw import encode, decode, bit_to_number, calculus_division, approximate_capacity
from dsw import get_complete_accessor, accessor_to_adjacency_matrix, adjacency_matrix_to_accessor

from experiments import local_bio_filters
from experiments.code_encode import generate_database, trans_fountain, trans_yinyang, trans_hedges, trans_spiderweb
from experiments.code_repair import evaluate_repair_errors, calculate_minimum_reads, evaluate_established_depth
from experiments.code_repair import spiderweb_step, hedges_step


def generate_algorithms():
    if exists("./data/primal_graphs.pkl") and exists("./data/coding_graphs.pkl"):
        return

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

    with open("./data/primal_graphs.pkl", "wb") as file:
        dump(file=file, obj=primal_graphs)
    with open("./data/coding_graphs.pkl", "wb") as file:
        dump(file=file, obj=coding_graphs)


def evaluate_performance():
    if not exists("./data/capacity_reference.pkl"):
        with open("./data/primal_graphs.pkl", "rb") as file:
            primal_graphs = load(file=file)

        capacities, data = [], []
        for index in range(12):
            accessor = primal_graphs[str(index + 1).zfill(2)]
            capacity, processes = approximate_capacity(accessor=accessor, repeats=10, need_process=True)
            capacities.append(capacity)
            data.append((capacity, processes))

        with open("./data/capacity_reference.pkl", "wb") as file:
            dump((capacities, data), file)

    if not exists("./data/performance_evaluation_1.pkl"):
        dataset = generate_database(random_seed=2021, total_length=8 * 1024 * 64, times=100)

        records = []
        records += trans_fountain(dataset, shape=(2048, 256), index_length=32, random_seed=2021, param_number=100)
        records += trans_yinyang(dataset, shape=(4096, 128), index_length=16, random_seed=2021, param_number=100)
        records += trans_hedges(dataset, shape=(4096, 128), index_length=16)
        records += trans_spiderweb(dataset, shape=(4096, 128), index_length=16, random_seed=2021, param_number=100)

        with open("./data/performance_evaluation_1.pkl", "wb") as file:
            dump(obj=array(records), file=file)

    if not exists("./data/performance_evaluation_2.pkl"):
        monitor, task_seed, bit_length, repeats = Monitor(), 2021, 100, 100

        with open("./data/coding_graphs.pkl", "rb") as file:
            coding_graphs = load(file=file)

        transcode_results_1 = []
        for index in range(12):
            filter_index = str(index + 1).zfill(2)
            print("Evaluate coding algorithm " + filter_index + ".")
            accessor = coding_graphs[filter_index]
            vertices = obtain_vertices(accessor=accessor)
            random.seed(task_seed)
            encoded_matrix, nucleotide_numbers = random.randint(0, 2, size=(repeats, bit_length)), []
            for current, start_index in enumerate(vertices):
                dna_strings, nucleotide_number = [], 0
                for binary_message in encoded_matrix:
                    dna_string = encode(binary_message=binary_message, accessor=accessor, start_index=start_index)
                    nucleotide_number += len(dna_string)
                    dna_strings.append(dna_string)
                nucleotide_numbers.append(nucleotide_number)
                decoded_matrix = []  # temp the correctness of transcoding.
                for dna_string in dna_strings:
                    binary_message = decode(dna_string=dna_string, bit_length=bit_length,
                                            accessor=accessor, start_index=start_index)
                    decoded_matrix.append(binary_message)
                decoded_matrix = array(decoded_matrix)
                if any(encoded_matrix != decoded_matrix):
                    raise ValueError("Encoded matrix is not equal to decoded matrix!")

                monitor.output(current + 1, len(vertices))

            transcode_results_1.append(bit_length / array(nucleotide_numbers))

        with open("./data/primal_graphs.pkl", "rb") as file:
            primal_graphs = load(file=file)

        transcode_results_2 = []
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
                    available, nucleotide_number = mapping, 0
                    quotient, dna_string, monitor = bit_to_number(binary_message), "", Monitor()
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

                    monitor.output(current + 1, len(encoded_matrix) * 24)
                    current += 1

            transcode_results_2[filter_index] = bit_length / array(nucleotide_numbers)

            with open("./data/performance_evaluation_2.pkl", "wb") as file:
                dump(obj=(transcode_results_1, transcode_results_2), file=file)


def evaluate_correction():
    random_seed = 2021

    if not exists("./data/correction_evaluation_1.pkl"):
        print("Evaluate correction rate of the coding graph of constraint set 01.")
        records = {}
        with open("./data/coding_graphs.pkl", "rb") as file:
            accessor = load(file)["01"]
        vertices = obtain_vertices(accessor)
        for errors in [1, 2, 4, 6]:
            results, time_cost = evaluate_repair_errors(random_seed=random_seed, dna_length=200,
                                                        repeats=10000, errors=errors,
                                                        accessor=accessor, vertices=vertices, observed_length=10)
            records[errors / 200.0] = (results, time_cost / 10000.0)
        with open("./data/correction_evaluation_1.pkl", "wb") as file:
            dump(obj=records, file=file)

    if not exists("./data/correction_evaluation_2.pkl"):
        print("Evaluate minimum pretreatment-free depth of the coding graph of constraint set 01.")
        records = {}
        with open("./data/coding_graphs.pkl", "rb") as file:
            accessor = load(file)["01"]
        vertices = obtain_vertices(accessor)
        for errors in [1, 2, 4, 6]:
            results = calculate_minimum_reads(random_seed=random_seed, dna_length=200, repeats=10000,
                                              errors=errors, accessor=accessor, vertices=vertices, observed_length=10)
            records[errors / 200.0] = results
        with open("./data/correction_evaluation_2.pkl", "wb") as file:
            dump(obj=records, file=file)

    if not exists("./data/correction_evaluation_3.pkl"):
        print("Evaluate success rate of the coding graph of constraint set 01 with established depths under 2% errors.")
        records = []
        with open("./data/coding_graphs.pkl", "rb") as file:
            accessor = load(file)["01"]
        vertices = obtain_vertices(accessor)
        results = evaluate_established_depth(random_seed=random_seed, dna_length=200,
                                             accessor=accessor, vertices=vertices, observed_length=10,
                                             classes=10000, errors=4, established_depth=5)
        right_strings, obtained_set = results
        right_strings = set(right_strings)
        final_results = sorted(obtained_set.items(), key=lambda sample: sample[1], reverse=True)
        for obtained_string, count in final_results:
            records.append([obtained_string in right_strings, count])
        with open("./data/correction_evaluation_3.pkl", "wb") as file:
            dump(obj=array(records), file=file)

    if not exists("./data/correction_evaluation_4.pkl"):
        print("Evaluate success rate of the coding graph of constraint set 01 with established depths under 2% errors.")
        records = {0.005: zeros(shape=(20, 100, 3)), 0.010: zeros(shape=(20, 100, 3)),
                   0.020: zeros(shape=(20, 100, 3)), 0.030: zeros(shape=(20, 100, 3))}
        with open("./data/coding_graphs.pkl", "rb") as file:
            accessor = load(file)["01"]
        vertices = obtain_vertices(accessor)
        for error_index, errors in enumerate([1, 2, 4, 6]):
            for depth_index, depth in enumerate(arange(20) + 1):
                for repeat in range(100):
                    results = evaluate_established_depth(random_seed=random_seed + repeat, dna_length=200,
                                                         accessor=accessor, vertices=vertices, observed_length=10,
                                                         classes=10000, errors=errors, established_depth=depth)
                    right_strings, obtained_set = results
                    right_strings = set(right_strings)
                    final_results = sorted(obtained_set.items(), key=lambda sample: sample[1], reverse=True)
                    success, right_count, wrong_count = 0, None, None
                    for obtained_string, count in final_results[:10000]:
                        if obtained_string in right_strings:
                            success += 1
                            right_count = min([right_count, count]) if right_count is not None else count
                        else:
                            wrong_count = max([wrong_count, count]) if wrong_count is not None else count
                    if wrong_count is None:
                        for obtained_string, count in final_results[10000:]:
                            if obtained_string not in right_strings:
                                wrong_count = max([wrong_count, count]) if wrong_count is not None else count
                    records[errors / 200.0][depth_index, repeat] = [success / 10000, right_count, wrong_count]
        with open("./data/correction_evaluation_4.pkl", "wb") as file:
            dump(obj=records, file=file)

    if not exists("./data/correction_evaluation_5.pkl"):
        print("Evaluate the influence of constraint capacity on correction rate.")
        records = {}
        with open("./data/coding_graphs.pkl", "rb") as file:
            accessors = load(file)
        for set_index, accessor in accessors.items():
            print("Calculate constraint set " + set_index + ".")
            vertices = obtain_vertices(accessor)
            results = evaluate_repair_errors(random_seed=random_seed, dna_length=200, repeats=10000,
                                             errors=2, accessor=accessor, vertices=vertices, observed_length=10)
            records[set_index] = results
        with open("./data/correction_evaluation_5.pkl", "wb") as file:
            dump(obj=records, file=file)

    if not exists("./data/correction_evaluation_6.pkl"):
        records = {"spiderweb": {}, "hedges": {}}
        with open("./data/coding_graphs.pkl", "rb") as file:
            accessor = load(file)["01"]
        vertices = obtain_vertices(accessor)
        for errors in [0, 1, 2, 4, 6]:
            results = spiderweb_step(random_seed=random_seed, accessor=accessor, vertices=vertices,
                                     dna_length=200, observed_length=10, repeats=10000, errors=errors)
            records["spiderweb"][errors / 200.0] = results
        for errors in [0, 1, 2, 4, 6]:
            results = hedges_step(random_seed=2021, bio_filter=local_bio_filters["01"],
                                  dna_length=200, observed_length=10, repeats=10000, errors=errors)
            records["hedges"][errors / 200.0] = results
        with open("./data/correction_evaluation_6.pkl", "wb") as file:
            dump(obj=records, file=file)


def evaluate_reconstruction():
    if not exists("./data/reconstruction_evaluation.npy"):
        random_seed, sequencing_length, repeats, threshold = 2021, 100, 100, 10000
        nucleotides, reconstructions = ["A", "C", "G", "T"], []
        with open("./data/coding_graphs.pkl", "rb") as file:
            accessors = load(file)

        for index in range(12):
            filter_index = str(index + 1).zfill(2)
            accessor = accessors[filter_index]
            vertices = obtain_vertices(accessor)
            records = zeros(shape=(repeats, threshold + 1))
            random.seed(random_seed)
            for repeat in range(repeats):
                print("Do the reconstruction test for the constraint set " + filter_index
                      + " with sequencing random " + str(sequencing_length) + "bp as repeat " + str(repeat + 1))

                monitor = Monitor()
                visited, number, total, used = zeros(shape=(4 ** 10,), dtype=bool), 0, len(vertices), []
                while sum(visited) < total:
                    vertex_index = random.choice(vertices)
                    dna_string = ""
                    number += 1
                    for _ in range(sequencing_length):
                        used_index = random.choice(where(accessor[vertex_index] >= 0)[0])
                        nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
                        visited[vertex_index] = True
                        dna_string += nucleotide
                    used.append(dna_string)

                    records[repeat, number] = sum(visited)
                    monitor.output(current_state=sum(visited), total_state=total,
                                   extra={"number of DNA string": number})

                    if number == threshold:
                        monitor.output(current_state=total, total_state=total,
                                       extra={"number of DNA string": "> " + str(threshold)})
                        break

            random.seed(None)

        reconstructions = array(reconstructions)

        save(file="./data/reconstruction_evaluation.npy", arr=reconstructions)


def evaluate_reliability():
    if not exists("./data/capacity_evaluation.pkl"):
        test_times, terminal_length, random_seed, records = 100, 6, 2021, {}

        random.seed(random_seed)
        feed_matrix = accessor_to_adjacency_matrix(get_complete_accessor(observed_length=2))
        data, errors, monitor = [], [], Monitor()
        for time in range(test_times):
            matrix = feed_matrix.copy()
            for position in array(list(where(feed_matrix == 1)[:2])).T:
                if random.random(1)[0] <= 0.5:
                    matrix[position[0], position[1]] = 0
            monitor.output(current_state=time + 1, total_state=test_times)
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
                monitor.output(current_state=time + 1, total_state=test_times)
                calculated_eigenvalue = real(max(linalg.eig(matrix)[0]))
                calculated_value = log2(calculated_eigenvalue) if calculated_eigenvalue > 0.0 else 0.0
                approximate_value = approximate_capacity(accessor=adjacency_matrix_to_accessor(matrix), repeats=10)
                values.append(approximate_value)
                errors.append(abs(approximate_value - calculated_value))
            whole_data[length] = (values, errors)

        records["extend"] = whole_data

        with open("./data/capacity_evaluation.pkl", "wb") as file:
            dump(obj=records, file=file)


if __name__ == "__main__":
    # generate_algorithms()
    # evaluate_performance()
    evaluate_correction()
    # evaluate_reconstruction()
    # evaluate_reliability()
