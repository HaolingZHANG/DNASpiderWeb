# noinspection PyPackageRequirements
from Bio.pairwise2 import align
from datetime import datetime
from numpy import random, array, load, sum, std, max, mean, log, where
from warnings import filterwarnings

from dsw import Monitor, encode, set_vt, repair_dna, obtain_vertices
from dsw import accessor_to_latter_map, approximate_capacity, remove_nasty_arc

filterwarnings("ignore", category=RuntimeWarning)


def show_single_examples():
    with open("./data/coding_graphs.pkl", "rb") as file:
        accessor = load(file, allow_pickle=True)["01"]
    vertices = obtain_vertices(accessor)

    start_index = vertices[0]
    binary_message = [0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1]
    wrong_location = 14
    right_dna_string, vt_check = encode(binary_message=binary_message, accessor=accessor, start_index=start_index,
                                        vt_length=10)
    for nucleotide in list(filter(lambda n: n != right_dna_string[wrong_location], ["A", "C", "G", "T"])):
        wrong_dna_string = list(right_dna_string)
        wrong_dna_string[wrong_location] = nucleotide
        wrong_dna_string = "".join(wrong_dna_string)
        print("Right: " + right_dna_string)
        print("Wrong: " + wrong_dna_string)

        repaired_dna_strings, additions = repair_dna(dna_string=wrong_dna_string,
                                                     accessor=accessor, start_index=start_index,
                                                     observed_length=10, vt_check=vt_check,
                                                     has_insertion=True, has_deletion=True, verbose=True)
        if len(repaired_dna_strings) > 0:
            reference_length = max([len(repaired_dna_string) for repaired_dna_string in repaired_dna_strings])
            print()
            print("Repair strategy:")
            print("right DNA string    | " + right_dna_string.ljust(reference_length) + " | ")
            for repaired_dna_string in repaired_dna_strings:
                flag = right_dna_string == repaired_dna_string
                print("repaired DNA string | " + repaired_dna_string.ljust(reference_length) + " | " + str(flag))
        else:
            print("Cannot find any repair strategy.")
        print()
        print()

    start_index = vertices[0]
    binary_message = [0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1]
    wrong_location = 14
    right_dna_string, vt_check = encode(binary_message=binary_message, accessor=accessor, start_index=start_index,
                                        vt_length=10)
    for nucleotide in list(["A", "C", "G", "T"]):
        wrong_dna_string = list(right_dna_string)
        wrong_dna_string.insert(wrong_location, nucleotide)
        wrong_dna_string = "".join(wrong_dna_string)
        print("Right: " + right_dna_string)
        print("Wrong: " + wrong_dna_string)

        repaired_dna_strings, additions = repair_dna(dna_string=wrong_dna_string,
                                                     accessor=accessor, start_index=start_index,
                                                     observed_length=10, vt_check=vt_check,
                                                     has_insertion=True, has_deletion=True, verbose=True)
        if len(repaired_dna_strings) > 0:
            reference_length = max([len(repaired_dna_string) for repaired_dna_string in repaired_dna_strings])
            print()
            print("Repair strategy:")
            print("right DNA string    | " + right_dna_string.ljust(reference_length) + " | ")
            for repaired_dna_string in repaired_dna_strings:
                flag = right_dna_string == repaired_dna_string
                print("repaired DNA string | " + repaired_dna_string.ljust(reference_length) + " | " + str(flag))
        else:
            print("Cannot find any repair strategy.")
        print()
        print()

    start_index = vertices[0]
    binary_message = [0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1]
    wrong_location = 14
    right_dna_string, vt_check = encode(binary_message=binary_message, accessor=accessor, start_index=start_index,
                                        vt_length=10)
    wrong_dna_string = list(right_dna_string)
    del wrong_dna_string[wrong_location]
    wrong_dna_string = "".join(wrong_dna_string)
    print("Right: " + right_dna_string)
    print("Wrong: " + wrong_dna_string)

    repaired_dna_strings, additions = repair_dna(dna_string=wrong_dna_string,
                                                 accessor=accessor, start_index=start_index,
                                                 observed_length=10, vt_check=vt_check,
                                                 has_insertion=True, has_deletion=True, verbose=True)
    if len(repaired_dna_strings) > 0:
        reference_length = max([len(repaired_dna_string) for repaired_dna_string in repaired_dna_strings])
        print()
        print("Repair strategy:")
        print("right DNA string    | " + right_dna_string.ljust(reference_length) + " | ")
        for repaired_dna_string in repaired_dna_strings:
            flag = right_dna_string == repaired_dna_string
            print("repaired DNA string | " + repaired_dna_string.ljust(reference_length) + " | " + str(flag))
    else:
        print("Cannot find any repair strategy.")
    print()
    print()


def show_multiple_examples():
    with open("./data/coding_graphs.pkl", "rb") as file:
        accessor = load(file, allow_pickle=True)["01"]
    vertices = obtain_vertices(accessor)
    dna_length, observed_length, nucleotides, error_times = 50, 10, ["A", "C", "G", "T"], 2
    random.seed(2021)

    for current_repeat, start_index in enumerate(vertices[:4]):
        # create ideal DNA string randomly.
        vertex_index, right_dna_string = start_index, ""
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            right_dna_string += nucleotide

        right_vt_check = set_vt(dna_string=right_dna_string, vt_length=observed_length + 1, nucleotides=nucleotides)

        # introduce error under most dense.
        wrong_dna_string = list(right_dna_string)
        for _ in range(error_times):
            wrong_type = random.randint(0, 3)
            error_location = random.randint(observed_length, len(wrong_dna_string) - observed_length - 1)
            if wrong_type == 0:  # substitution
                nucleotide = random.choice(list(filter(lambda n: n != wrong_dna_string[error_location], nucleotides)))
                wrong_dna_string[error_location] = nucleotide
            elif wrong_type == 1:  # insertion
                nucleotide = random.choice(nucleotides)
                wrong_dna_string.insert(error_location, nucleotide)
            else:  # deletion
                del wrong_dna_string[error_location]
        wrong_dna_string = "".join(wrong_dna_string)

        max_score, shown_right, shown_wrong = 0, "", ""
        for alignment in align.localxx(right_dna_string, wrong_dna_string):
            if alignment.score > max_score:
                max_score, shown_right, shown_wrong = alignment.score, alignment.seqA, alignment.seqB
        difference = sum(array([bool(right != wrong) for right, wrong in zip(shown_right, shown_wrong)]))
        error_rate = difference / len(shown_right)

        print("task (" + str(current_repeat + 1) + " / 4) " + "with error rate %.2f" % (error_rate * 100) + "%")
        print("right DNA string = " + shown_right)
        print("wrong DNA string = " + shown_wrong)
        print()

        repaired_dna_strings, additions = repair_dna(dna_string=wrong_dna_string,
                                                     accessor=accessor, start_index=start_index,
                                                     observed_length=observed_length,
                                                     vt_check=right_vt_check,
                                                     check_iterations=error_times + 1,
                                                     has_insertion=True, has_deletion=True)

        print("find error(s)    = " + str(additions[0] or additions[1]))
        if len(repaired_dna_strings) > 0:
            print("find repair(s)   = True")
            reference_length = max([len(dna_string) for dna_string in repaired_dna_strings])
            for repaired_dna_string in repaired_dna_strings:
                current_vt_check = set_vt(dna_string=repaired_dna_string, vt_length=observed_length + 1,
                                          nucleotides=nucleotides)
                correction_flag = (right_vt_check == current_vt_check)
                shown_check = repaired_dna_string.ljust(reference_length + 1)
                shown_flag = str(correction_flag).ljust(5)
                situation = repaired_dna_string == right_dna_string
                print("check DNA string = " + shown_check + " | " + shown_flag + " | " + str(situation))
        else:
            print("find repair(s)   = False")

        print()

    random.seed(None)


def introduce_errors(right_dna_string, error_times, observed_length, nucleotides=None):
    if error_times == 0:
        return right_dna_string

    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    while True:
        wrong_dna_string = list(right_dna_string)  # introduce error under most dense.

        for _ in range(error_times):
            wrong_type = random.randint(0, 3)  # the probability of three types of errors is the same.
            error_location = random.randint(observed_length + 1, len(wrong_dna_string) - observed_length - 1)
            if wrong_type == 0:  # substitution
                nucleotide = random.choice(list(filter(lambda n: n != wrong_dna_string[error_location],
                                                       nucleotides)))
                wrong_dna_string[error_location] = nucleotide
            elif wrong_type == 1:  # insertion
                nucleotide = random.choice(nucleotides)
                wrong_dna_string.insert(error_location, nucleotide)
            else:  # deletion
                del wrong_dna_string[error_location]
        wrong_dna_string = "".join(wrong_dna_string)

        if wrong_dna_string != right_dna_string:
            return wrong_dna_string


def evaluate_repair_multiple_errors(random_seed, accessor, vertices, observed_length,
                                    repeats, dna_length, error_times, check_iterations, nucleotides=None):
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    print("Check DNA string of length " + str(dna_length) + " imposing " + str(error_times) + " error(s).")
    random.seed(random_seed)
    records, monitor = [], Monitor()
    for current_repeat in range(repeats):
        start_index = random.choice(vertices)

        # create ideal DNA string randomly.
        vertex_index, right_dna_string = start_index, ""
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            right_dna_string += nucleotide

        vt_check = set_vt(dna_string=right_dna_string, vt_length=observed_length + 1, nucleotides=nucleotides)

        wrong_dna_string = introduce_errors(right_dna_string=right_dna_string,
                                            error_times=error_times, observed_length=observed_length,
                                            nucleotides=nucleotides)

        max_score, shown_right, shown_wrong = 0, "", ""
        for alignment in align.localxx(right_dna_string, wrong_dna_string):
            if alignment.score > max_score:
                max_score, shown_right, shown_wrong = alignment.score, alignment.seqA, alignment.seqB

        differences = where(array([bool(right != wrong) for right, wrong in zip(shown_right, shown_wrong)]))[0]
        actual_error_rate = len(differences) / len(shown_right)
        differences = differences / 1000.0
        variation_coefficient = std(differences) / mean(differences)  # variation coefficient

        previous_time = datetime.now()
        # repair with check list of Varshamov-Tenengolts code.
        repaired_dna_strings, additions = repair_dna(dna_string=wrong_dna_string,
                                                     accessor=accessor, start_index=start_index,
                                                     observed_length=observed_length, vt_check=vt_check,
                                                     check_iterations=check_iterations,
                                                     has_insertion=True, has_deletion=True)
        time_cost = (datetime.now() - previous_time).total_seconds()

        records.append([actual_error_rate, variation_coefficient,
                        additions[0], additions[2], additions[1], additions[0] or additions[1],
                        len(repaired_dna_strings), right_dna_string in repaired_dna_strings, time_cost])
        monitor.output(current_repeat + 1, repeats)

    random.seed(None)

    return array(records)


def calculate_spiderweb_minimum_reads(random_seed, accessor, vertices, observed_length, repeats, dna_length,
                                      error_times, check_iterations, nucleotides=None):
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    print("Check DNA string of length " + str(dna_length) + " imposing " + str(error_times) + " error(s).")
    random.seed(random_seed)
    records, monitor, maximum_reads = [], Monitor(), 1
    for current_repeat in range(repeats):
        start_index = random.choice(vertices)

        # create ideal DNA string randomly.
        vertex_index, right_dna_string = start_index, ""
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            right_dna_string += nucleotide

        vt_check = set_vt(dna_string=right_dna_string, vt_length=observed_length + 1, nucleotides=nucleotides)

        reads, result_set = 1, {}  # zeros(shape=(4, dna_length), dtype=int)
        while True:
            wrong_dna_string = introduce_errors(right_dna_string=right_dna_string,
                                                error_times=error_times, observed_length=observed_length,
                                                nucleotides=nucleotides)

            # repair with check list of Varshamov-Tenengolts code.
            repaired_dna_strings, additions = repair_dna(dna_string=wrong_dna_string,
                                                         accessor=accessor, start_index=start_index,
                                                         observed_length=observed_length, vt_check=vt_check,
                                                         check_iterations=check_iterations,
                                                         has_insertion=True, has_deletion=True)
            for repaired_dna_string in repaired_dna_strings:
                if repaired_dna_string in result_set:
                    result_set[repaired_dna_string] += 1
                else:
                    result_set[repaired_dna_string] = 1

            strings, counts = [], []
            for string, count in result_set.items():
                strings.append(string)
                counts.append(count)
            counts = array(counts)

            if where(counts == max(counts))[0] == 1:
                best_dna_string = strings[where(counts == max(counts))[0][0]]
                if best_dna_string == right_dna_string:
                    maximum_reads = max([maximum_reads, reads])
                    break

            reads += 1

        monitor.output(current_repeat + 1, repeats, extra={"current reads": reads, "maximum reads": maximum_reads})

        records.append(reads)

    random.seed(None)

    return records


def evaluate_spiderweb_established_reads(random_seed, accessor, vertices, observed_length, classes, dna_length,
                                         error_times, check_iterations, established_reads, nucleotides=None):
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    print("Check DNA string of length " + str(dna_length) + " imposing " + str(error_times) + " error(s).")
    random.seed(random_seed)
    right_set, obtained_set, monitor = [], {}, Monitor()
    while True:
        start_index = random.choice(vertices)

        # create ideal DNA string randomly.
        vertex_index, right_dna_string = start_index, ""
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            right_dna_string += nucleotide

        if right_dna_string not in right_set:
            right_set.append(right_dna_string)
        else:
            continue

        vt_check = set_vt(dna_string=right_dna_string, vt_length=observed_length + 1, nucleotides=nucleotides)

        for read in range(established_reads):
            wrong_dna_string = introduce_errors(right_dna_string=right_dna_string,
                                                error_times=error_times, observed_length=observed_length,
                                                nucleotides=nucleotides)

            # repair with check list of Varshamov-Tenengolts code.
            repaired_dna_strings, additions = repair_dna(dna_string=wrong_dna_string,
                                                         accessor=accessor, start_index=start_index,
                                                         observed_length=observed_length, vt_check=vt_check,
                                                         check_iterations=check_iterations,
                                                         has_insertion=True, has_deletion=True)
            for repaired_dna_string in repaired_dna_strings:
                if repaired_dna_string in obtained_set:
                    obtained_set[repaired_dna_string] += 1
                else:
                    obtained_set[repaired_dna_string] = 1

        monitor.output(len(right_set), classes)

        if len(right_set) == classes:
            break

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
                                                             repair_insertion=True, repair_deletion=True)
        new_code_rate = approximate_capacity(accessor=accessor, repeats=4, maximum_iteration=100, verbose=False)

        print("Current code rate is %.5f.\n" % float(new_code_rate))

        vertices = obtain_vertices(accessor=accessor)
        correction_data = {}
        for error_time in [1, 2, 3, 4]:
            records, runtimes = evaluate_repair_multiple_errors(random_seed=2021, error_times=error_time,
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
