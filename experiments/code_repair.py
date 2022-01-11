# noinspection PyPackageRequirements
from Bio.pairwise2 import align
from numpy import random, array, load, sum, where

from dsw import Monitor, encode, repair_dna, number_to_dna


def show_single_examples():
    accessor, vertices = load(file="./results/data/a01[g].npy"), load(file="./results/data/a01[v].npy")

    start_index = where(vertices == 1)[0][0]
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

        detect_flag, repaired_dna_strings = repair_dna(dna_string=wrong_dna_string,
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

    start_index = where(vertices == 1)[0][0]
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

        detect_flag, repaired_dna_strings = repair_dna(dna_string=wrong_dna_string,
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

    start_index = where(vertices == 1)[0][0]
    binary_message = [0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1]
    wrong_location = 14
    right_dna_string, vt_check = encode(binary_message=binary_message, accessor=accessor, start_index=start_index,
                                        vt_length=10)
    wrong_dna_string = list(right_dna_string)
    del wrong_dna_string[wrong_location]
    wrong_dna_string = "".join(wrong_dna_string)
    print("Right: " + right_dna_string)
    print("Wrong: " + wrong_dna_string)

    detect_flag, repaired_dna_strings = repair_dna(dna_string=wrong_dna_string,
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
    accessor, vertices = load("./results/data/a01[g].npy"), where(load("./results/data/a01[v].npy") == 1)[0]
    dna_length, vt_length, observed_length, nucleotides, error_times = 50, 10, 10, ["A", "C", "G", "T"], 2
    random.seed(2021)

    for current_repeat, start_index in enumerate(vertices[:4]):
        # create ideal DNA string randomly.
        vertex_index, right_dna_string, vt_value = start_index, "", 0
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            vt_value += (used_index * location) % (4 ** vt_length)
            right_dna_string += nucleotide
        right_vt_check = number_to_dna(decimal_number=int(vt_value), dna_length=vt_length)

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

        detect_flag, repaired_dna_strings = repair_dna(dna_string=wrong_dna_string,
                                                       accessor=accessor, start_index=start_index,
                                                       observed_length=observed_length,
                                                       check_iterations=2, has_insertion=True, has_deletion=True)

        print("find error(s)    = " + str(detect_flag))
        if len(repaired_dna_strings) > 0:
            print("find repair(s)   = True")
            reference_length = max([len(dna_string) for dna_string in repaired_dna_strings])
            for repaired_dna_string in repaired_dna_strings:
                vertex_index, vt_value = start_index, 0
                for location, nucleotide in enumerate(repaired_dna_string):
                    vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                    vt_value += (nucleotides.index(nucleotide) * location) % (4 ** len(right_vt_check))
                vt_check = number_to_dna(decimal_number=vt_value, dna_length=vt_length)

                correction_flag = (right_vt_check == vt_check)
                shown_check = repaired_dna_string.ljust(reference_length + 1)
                shown_flag = str(correction_flag).ljust(5)
                situation = repaired_dna_string == right_dna_string
                print("check DNA string = " + shown_check + " | " + shown_flag + " | " + str(situation))
        else:
            print("find repair(s)   = False")

        print()

    random.seed(None)


def evaluate_single_error(accessor, start_indices, repeats):
    monitor, current, total, records = Monitor(), 1, repeats * len(start_indices), []
    nucleotides, observed_length, vt_length, error_position, dna_length = ["A", "C", "G", "T"], 10, 10, 12, 30
    for _ in range(repeats):
        for start_index in start_indices:
            # create ideal DNA string randomly.
            vertex_index, right_dna_string, vt_value = start_index, "", 0
            for location in range(30):
                used_indices = where(accessor[vertex_index] >= 0)[0]
                used_index = random.choice(used_indices)
                nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
                vt_value += used_index * (location + 1)
                right_dna_string += nucleotide
            vt_value %= 4 ** vt_length
            vt_check = number_to_dna(decimal_number=int(vt_value), dna_length=vt_length)

            # occur substitution (3 types)
            for nucleotide in list(filter(lambda n: n != right_dna_string[error_position], nucleotides)):
                wrong_dna_string = list(right_dna_string)
                wrong_dna_string[error_position] = nucleotide
                wrong_dna_string = "".join(wrong_dna_string)

                # repair without check list of Varshamov-Tenengolts code.
                detect_flag_1, repaired_dna_strings_1 = repair_dna(dna_string=wrong_dna_string,
                                                                   accessor=accessor, start_index=start_index,
                                                                   observed_length=observed_length,
                                                                   has_insertion=True, has_deletion=True, verbose=False)
                # repair with check list of Varshamov-Tenengolts code.
                detect_flag_2, repaired_dna_strings_2 = repair_dna(dna_string=wrong_dna_string,
                                                                   accessor=accessor, start_index=start_index,
                                                                   observed_length=observed_length, vt_check=vt_check,
                                                                   has_insertion=True, has_deletion=True, verbose=False)
                correct_flag = len(repaired_dna_strings_2) == 1 and repaired_dna_strings_2[0] == right_dna_string
                records.append([current, 0, detect_flag_1, len(repaired_dna_strings_1), detect_flag_2, correct_flag])

            # occur insertion (4 types)
            for nucleotide in list(["A", "C", "G", "T"]):
                wrong_dna_string = list(right_dna_string)
                wrong_dna_string.insert(error_position, nucleotide)
                wrong_dna_string = "".join(wrong_dna_string)
                # repair without check list of Varshamov-Tenengolts code.
                detect_flag_1, repaired_dna_strings_1 = repair_dna(dna_string=wrong_dna_string,
                                                                   accessor=accessor, start_index=start_index,
                                                                   observed_length=observed_length,
                                                                   has_insertion=True, has_deletion=True, verbose=False)
                # repair with check list of Varshamov-Tenengolts code.
                detect_flag_2, repaired_dna_strings_2 = repair_dna(dna_string=wrong_dna_string,
                                                                   accessor=accessor, start_index=start_index,
                                                                   observed_length=observed_length, vt_check=vt_check,
                                                                   has_insertion=True, has_deletion=True, verbose=False)
                correct_flag = len(repaired_dna_strings_2) == 1 and repaired_dna_strings_2[0] == right_dna_string
                records.append([current, 1, detect_flag_1, len(repaired_dna_strings_1), detect_flag_2, correct_flag])

            # occur deletion (1 type)
            wrong_dna_string = list(right_dna_string)
            del wrong_dna_string[error_position]
            wrong_dna_string = "".join(wrong_dna_string)
            # repair without check list of Varshamov-Tenengolts code.
            detect_flag_1, repaired_dna_strings_1 = repair_dna(dna_string=wrong_dna_string,
                                                               accessor=accessor, start_index=start_index,
                                                               observed_length=observed_length,
                                                               has_insertion=True, has_deletion=True, verbose=False)
            # repair with check list of Varshamov-Tenengolts code.
            detect_flag_2, repaired_dna_strings_2 = repair_dna(dna_string=wrong_dna_string,
                                                               accessor=accessor, start_index=start_index,
                                                               observed_length=observed_length, vt_check=vt_check,
                                                               has_insertion=True, has_deletion=True, verbose=False)
            correct_flag = len(repaired_dna_strings_2) == 1 and repaired_dna_strings_2[0] == right_dna_string
            records.append([current, 2, detect_flag_1, len(repaired_dna_strings_1), detect_flag_2, correct_flag])

            monitor.output(current, total)
            current += 1

    return array(records).astype(int)


def evaluate_repair_multiple_errors(random_seed, accessor, vertices, observed_length, vt_length,
                                    repeats, dna_length, error_times, check_iterations, nucleotides=None):
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    print("Check DNA string of length " + str(dna_length) + " imposing " + str(error_times) + " error(s).")
    random.seed(random_seed)
    records, monitor = [], Monitor()
    for current_repeat in range(repeats):
        start_index = random.choice(vertices)

        # create ideal DNA string randomly.
        vertex_index, right_dna_string, vt_value = start_index, "", 0
        for location in range(dna_length):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            used_index = random.choice(used_indices)
            nucleotide, vertex_index = nucleotides[used_index], accessor[vertex_index][used_index]
            vt_value += used_index * (location + 1)
            right_dna_string += nucleotide
        vt_value %= 4 ** vt_length
        right_vt_check = number_to_dna(decimal_number=int(vt_value), dna_length=vt_length)

        # introduce error under most dense.
        wrong_dna_string = list(right_dna_string)
        for _ in range(error_times):
            wrong_type = random.randint(0, 3)
            error_location = random.randint(observed_length + 1, len(wrong_dna_string) - observed_length - 1)
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

        # repair with check list of Varshamov-Tenengolts code.
        detect_flag, repaired_dna_strings = repair_dna(dna_string=wrong_dna_string,
                                                       accessor=accessor, start_index=start_index,
                                                       observed_length=observed_length, vt_check=right_vt_check,
                                                       check_iterations=check_iterations,
                                                       has_insertion=True, has_deletion=True)

        correct_flag, false_positive_count = False, 0
        for repaired_dna_string in repaired_dna_strings:
            if repaired_dna_string == right_dna_string:
                correct_flag = True
            else:
                false_positive_count += 1

        records.append((error_rate, shown_right, shown_wrong, detect_flag, correct_flag, false_positive_count))
        monitor.output(current_repeat + 1, repeats)

    random.seed(None)

    return records
