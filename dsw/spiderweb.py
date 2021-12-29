__author__ = "Zhang, Haoling [hlzchn@gmail.com]"

from numpy import zeros, ones, array, random, sum, log, argsort, where

from dsw.operation import Monitor, calculus_addition, calculus_multiplication, calculus_division
from dsw.operation import bit_to_number, number_to_bit, number_to_dna
from dsw.graphized import obtain_latters, repair_by_match


def encode(binary_message, accessor, start_index, nucleotides=None, shuffles=None, removed_edges=None):
    """
    Encode a bit array by the specific accessor.

    :param binary_message: binary message.
    :type binary_message: numpy.ndarray

    :param accessor: accessor of the coding algorithm.
    :type accessor: numpy.ndarray

    :param start_index: virtual vertex to start encoding.
    :type start_index: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :param shuffles: shuffle relationships for bit-nucleotide mapping.
    :type: numpy.ndarray

    :param removed_edges: removed directed edges in accessor.
    :type removed_edges: list

    :return: DNA string encoded by this graph.
    :rtype: str

    ..example::
        >>> from numpy import array
        >>> from dsw import encode
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> binary_message = array([0, 1, 0, 1, 0, 1, 0, 1])
        >>> encode(accessor=accessor, binary_message=binary_message, start_index=1)
        'TCTCTCT'
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    used_accessor = accessor.copy()
    if removed_edges is not None:  # remove some directed edges if required
        for vertex_index, location in removed_edges:
            used_accessor[vertex_index, location] = -1

    quotient, dna_string, vertex_index = bit_to_number(binary_message), "", start_index

    while quotient != "0":
        used_indices = where(used_accessor[vertex_index] >= 0)[0]

        if len(used_indices) > 1:  # current vertex contains information.
            quotient, remainder = calculus_division(number=quotient, base=str(len(used_indices)))
            remainder = int(remainder)

            if shuffles is not None:  # shuffle remainder based on the inputted shuffles.
                remainder = argsort(shuffles[vertex_index, used_indices])[remainder]

            value = used_indices[remainder]

        elif len(used_indices) == 1:  # current vertex does not contain information.
            value = used_indices[0]

        else:  # current vertex is wrong.
            raise ValueError("Current vertex doesn't have an out-degree, the accessor or the start vertex is wrong!")

        nucleotide, vertex_index = nucleotides[value], used_accessor[vertex_index][value]

        dna_string += nucleotide

    return dna_string


def decode(dna_string, bit_length, accessor, start_index, nucleotides=None, shuffles=None, removed_edges=None):
    """
    Decode a DNA string by the specific accessor.

    :param dna_string: DNA string encoded bt this graph.
    :type dna_string: str

    :param bit_length: length of the bit array.
    :type bit_length: int

    :param accessor: accessor of the coding algorithm (consistent with the encoding process).
    :type accessor: numpy.ndarray

    :param start_index: virtual vertex to start decoding (consistent with the encoding process).
    :type start_index: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :param shuffles: shuffle relationships for bit-nucleotide mapping.
    :type: numpy.ndarray

    :param removed_edges: removed directed edges in accessor.
    :type removed_edges: list

    :return: binary message decoded by this graph.
    :rtype: numpy.ndarray

    ..example::
        >>> from numpy import array
        >>> from dsw import decode
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> dna_string = "TCTCTCT"
        >>> decode(accessor=accessor, dna_string=dna_string, start_index=1, bit_length=8)
        array([0, 1, 0, 1, 0, 1, 0, 1])
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    used_accessor = accessor.copy()
    if removed_edges is not None:  # remove some directed edges if required
        for vertex_index, location in removed_edges:
            used_accessor[vertex_index, location] = -1

    quotient, saved_values, vertex_index = "0", [], start_index

    for index, nucleotide in enumerate(dna_string):
        used_indices = where(used_accessor[vertex_index] >= 0)[0]

        if len(used_indices) > 1:  # current vertex contains information.
            used_nucleotides = [nucleotides[used_index] for used_index in used_indices]

            if nucleotide in used_nucleotides:  # check whether the DNA string is right currently.
                remainder = used_nucleotides.index(nucleotide)
            else:
                raise ValueError("At least one error is found in this DNA string!")

            if shuffles is not None:  # shuffle remainder based on the inputted shuffles.
                remainder = where(argsort(shuffles[vertex_index, used_indices]) == remainder)[0][0]

            saved_values.append((len(used_indices), remainder))
            vertex_index = used_accessor[vertex_index][nucleotides.index(nucleotide)]

        elif len(used_indices) == 1:  # current vertex does not contain information.
            used_nucleotide = nucleotides[used_indices[0]]
            if nucleotide == used_nucleotide:
                vertex_index = used_accessor[vertex_index][nucleotides.index(nucleotide)]
            else:
                raise ValueError("At least one error is found in this DNA string!")

        else:  # current vertex is wrong.
            raise ValueError("Current vertex doesn't have an out-degree, the accessor or the start vertex is wrong!")

    for index, (out_degree, number) in enumerate(saved_values[::-1]):
        quotient = calculus_multiplication(number=quotient, base=str(out_degree))
        quotient = calculus_addition(number=quotient, base=str(number))

    return array(number_to_bit(decimal_number=quotient, bit_length=bit_length), dtype=int)


def repair_dna(dna_string, accessor, start_index, has_insertion=False, has_deletion=False,
               nucleotides=None, verbose=False):
    """
    Repair the DNA string containing one (or more) errors.

    :param dna_string: DNA string waiting for recovery.
    :type dna_string: str

    :param accessor: accessor of the coding algorithm.
    :type accessor: numpy.ndarray

    :param start_index: virtual vertex to start encoding.
    :type start_index: int

    :param has_insertion: consider insertion error.
    :type has_insertion: bool

    :param has_deletion: consider deletion error.
    :type has_deletion: bool

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :param verbose: need to print log.
    :type verbose: bool

    :return: repaired DNA strings (may contain multiple repair results).
    :rtype: dict

    ..example::
        >>> from numpy import array
        >>> from dsw import decode
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> dna_string = "TCTCTATCTCTC"  # "TCTCTCTCTCTC" is original DNA string
        >>> repaired_dna_strings = repair_dna(dna_string=dna_string, accessor=accessor, start_index=1, \
                                              has_insertion=True, has_deletion=True, verbose=False)
        >>> repaired_dna_strings["S"]  # repair by substitution
        [(5, 'C', 'TCTCTCTCTCTC'), (5, 'G', 'TCTCTGTCTCTC'), (3, 'G', 'TCTGTATCTCTC'), (2, 'A', 'TCACTATCTCTC')]
        >>> repaired_dna_strings["I"]  # repair by insertion
        []
        >>> repaired_dna_strings["D"]  # repair by deletion
        [(4, 'T', 'TCTCATCTCTC')]
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    observed_length = int(log(len(accessor)) / log(len(nucleotides)))

    vertex_index, index_queue = start_index, [start_index]
    for index, nucleotide in enumerate(dna_string):
        used_indices = where(accessor[vertex_index] >= 0)[0]
        used_nucleotides = [nucleotides[used_index] for used_index in used_indices]

        if nucleotide in used_nucleotides:  # temp whether the dna_string is right currently.
            vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
            index_queue.append(vertex_index)
        else:
            if verbose:
                used_indices = where(accessor[index_queue[-1]] >= 0)[0]
                o_options = [nucleotides[used_index] for used_index in used_indices]
                b_options = accessor[index_queue[-1]][~(accessor[index_queue[-1]] == -1)]
                print("error occur in location " + str(index) + ", found " + dna_string[index]
                      + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))

            repaired_dna_strings = {"S": [], "I": [], "D": []}

            if verbose:
                used_indices = where(accessor[index_queue[index]] >= 0)[0]
                o_options = [nucleotides[used_index] for used_index in used_indices]
                b_options = accessor[index_queue[index]][~(accessor[index_queue[index]] == -1)]
                print("-" * 50)
                print("Search error in current location.")
                print("Check location = " + dna_string[index] + " | " + str(index) + ", found " + dna_string[index]
                      + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))

            # probe search
            probe_dna_strings = repair_by_match(dna_string=dna_string, accessor=accessor, index_queue=index_queue,
                                                occur_location=index, nucleotides=nucleotides,
                                                has_insertion=has_insertion, has_deletion=has_deletion, verbose=verbose)

            for repair_type, location, repaired_nucleotide, repaired_dna_string in probe_dna_strings:
                if (location, repaired_nucleotide, repaired_dna_string) not in repaired_dna_strings[repair_type]:
                    repaired_dna_strings[repair_type].append((location, repaired_nucleotide, repaired_dna_string))

                    if len(repaired_dna_strings["D"]) is not None:
                        print(index_queue)

            # recall search
            recall_queue = index_queue[::-1][:observed_length + 1]
            for recall_index, original_vertex_index in enumerate(recall_queue):
                error = index - recall_index - 1
                if verbose:
                    print("-" * 50)
                    print("Search the error in the previous (" + str((recall_index + 1)) + ") location.")
                    print("Check location = " + dna_string[error] + " | " + str(error))

                recall_dna_strings = repair_by_match(dna_string=dna_string, accessor=accessor, index_queue=index_queue,
                                                     occur_location=error, nucleotides=nucleotides,
                                                     has_insertion=has_insertion, has_deletion=has_deletion,
                                                     verbose=verbose)

                for repair_type, location, repaired_nucleotide, repaired_dna_string in recall_dna_strings:
                    if (location, repaired_nucleotide, repaired_dna_string) not in repaired_dna_strings[repair_type]:
                        repaired_dna_strings[repair_type].append((location, repaired_nucleotide, repaired_dna_string))

                        if len(repaired_dna_strings["D"]) is not None:
                            print(index_queue)

            return repaired_dna_strings


def find_vertices(observed_length, bio_filter, nucleotides=None, verbose=False):
    """
    Find valid vertices based on the given constraints.

    :param observed_length: length of the DNA string in a vertex.
    :type observed_length: int

    :param bio_filter: screening operation for identifying the valid DNA string (required the given constraints).
    :type bio_filter: dsw.biofilter.DefaultBioFilter

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :param verbose: need to print log.
    :type verbose: bool

    :return: available vertices.
    :rtype: numpy.ndarray

    ..example::
        >>> from dsw import LocalBioFilter, find_vertices
        >>> bio_filter = LocalBioFilter(max_homopolymer_runs=2, gc_range=[0.5, 0.5])
        >>> # "True" refers to the available index of vertex.
        >>> find_vertices(observed_length=2, bio_filter=bio_filter)
        array([False,  True,  True, False,  True, False, False,  True,  True,
               False, False,  True, False,  True,  True, False])
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    vertices, monitor = zeros(shape=(int(len(nucleotides) ** observed_length),), dtype=bool), Monitor()

    if verbose:
        print("Find valid vertices in this observed length of DNA string.")

    for vertex_index in range(len(vertices)):
        dna_string = number_to_dna(decimal_number=vertex_index, dna_length=observed_length)
        vertices[vertex_index] = bio_filter.valid(dna_string=dna_string)

        if verbose:
            monitor.output(vertex_index + 1, len(vertices))

    valid_rate = sum(vertices) / len(vertices)

    if valid_rate == 0:
        raise ValueError("No vertex is collected!")

    if verbose:
        print(str(round(valid_rate * 100, 2)) + "% (" + str(sum(vertices)) + ") valid vertices are collected.")

    return vertices


def connect_valid_graph(observed_length, vertices, nucleotides=None, verbose=False):
    """
    Connect a valid graph by valid vertices.

    :param observed_length: length of the DNA string in a vertex.
    :type observed_length: int

    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :param verbose: need to print log.
    :type verbose: bool

    :return: accessor of the valid graph.
    :rtype: numpy.ndarray

    ..example::
        >>> from numpy import array
        >>> from dsw import connect_valid_graph
        >>> vertices = array([False,  True,  True, False,  True, False, False,  True,  True, \
                              False, False,  True, False,  True,  True, False])
        >>> connect_valid_graph(observed_length=2, vertices=vertices)
        array([[-1, -1, -1, -1],
               [ 4, -1, -1,  7],
               [ 8, -1, -1, 11],
               [-1, -1, -1, -1],
               [-1,  1,  2, -1],
               [-1, -1, -1, -1],
               [-1, -1, -1, -1],
               [-1, 13, 14, -1],
               [-1,  1,  2, -1],
               [-1, -1, -1, -1],
               [-1, -1, -1, -1],
               [-1, 13, 14, -1],
               [-1, -1, -1, -1],
               [ 4, -1, -1,  7],
               [ 8, -1, -1, 11],
               [-1, -1, -1, -1]])
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    if vertices is None:
        raise ValueError("No collected vertex!")

    if verbose:
        print("Connect valid graph with valid vertices.")

    valid_rate, monitor = sum(vertices) / len(vertices), Monitor()

    if valid_rate > 0:
        accessor = -ones(shape=(int(len(nucleotides) ** observed_length), len(nucleotides)), dtype=int)

        for vertex_index in range(int(len(nucleotides) ** observed_length)):
            if vertices[vertex_index]:
                latters = obtain_latters(current=vertex_index, observed_length=observed_length)
                for position, latter_vertex_index in enumerate(latters):
                    if vertices[latter_vertex_index]:
                        accessor[vertex_index][position] = latter_vertex_index

            if verbose:
                monitor.output(vertex_index + 1, len(vertices))

        if verbose:
            print("Valid graph is created.")

        return accessor
    else:
        raise ValueError("No collected vertex!")


def connect_coding_graph(observed_length, vertices, threshold, nucleotides=None, verbose=False):
    """
    Connect a coding algorithm by valid vertices and the threshold for minimum out-degree.

    :param observed_length: length of the DNA string in a vertex.
    :type observed_length: int

    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray

    :param threshold: threshold for minimum out-degree.
    :type threshold: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :param verbose: need to print log.
    :type verbose: bool

    :return: coding vertices and coding accessor.
    :rtype: (numpy.ndarray, numpy.ndarray)

    ..example::
        >>> from numpy import array
        >>> from dsw import connect_coding_graph
        >>> vertices = array([False,  True,  True, False,  True, False, False,  True,  True, \
                              False, False,  True, False,  True,  True, False])
        >>> vertices, accessor = connect_coding_graph(observed_length=2, vertices=vertices, threshold=2)
        >>> accessor
        array([[-1, -1, -1, -1],
               [ 4, -1, -1,  7],
               [ 8, -1, -1, 11],
               [-1, -1, -1, -1],
               [-1,  1,  2, -1],
               [-1, -1, -1, -1],
               [-1, -1, -1, -1],
               [-1, 13, 14, -1],
               [-1,  1,  2, -1],
               [-1, -1, -1, -1],
               [-1, -1, -1, -1],
               [-1, 13, 14, -1],
               [-1, -1, -1, -1],
               [ 4, -1, -1,  7],
               [ 8, -1, -1, 11],
               [-1, -1, -1, -1]])
        >>> vertices
        array([False,  True,  True, False,  True, False, False,  True,  True,
               False, False,  True, False,  True,  True, False])
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    times = 1

    while True:
        if verbose:
            print("Check the vertex collection requirement in round " + str(times) + ".")

        new_vertices, monitor = zeros(shape=(int(len(nucleotides) ** observed_length),), dtype=bool), Monitor()
        saved_indices = where(vertices != 0)[0]
        for current, vertex_index in enumerate(saved_indices):
            latter_indices = obtain_latters(current=vertex_index, observed_length=observed_length)
            new_vertices[vertex_index] = sum(vertices[latter_indices]) >= threshold

            if verbose:
                monitor.output(current + 1, len(saved_indices))

        changed = sum(vertices) - sum(new_vertices)

        if verbose:
            print(str(round(sum(new_vertices) / len(vertices) * 100, 2)) + "% (" + str(sum(new_vertices)) + ") " +
                  "valid vertices are saved.")

        if sum(new_vertices) < 1:
            raise ValueError("No algorithmic graph is created!")

        if not changed:
            break

        vertices = new_vertices
        times += 1

    valid_rate = sum(vertices) / len(vertices)
    if valid_rate > 0:
        accessor = -ones(shape=(int(len(nucleotides) ** observed_length), len(nucleotides)), dtype=int)
        for vertex_index in range(int(len(nucleotides) ** observed_length)):
            if vertices[vertex_index]:
                latters = obtain_latters(current=vertex_index, observed_length=observed_length)
                for position, latter_vertex_index in enumerate(latters):
                    if vertices[latter_vertex_index]:
                        accessor[vertex_index][position] = latter_vertex_index

            if verbose:
                monitor.output(vertex_index + 1, len(vertices))

        if verbose:
            print("The coding graph is created.")

        return vertices, accessor
    else:
        raise ValueError("The coding graph cannot be created!")


def create_random_shuffles(observed_length, nucleotides=None, random_seed=None, verbose=False):
    """
    Create the shuffles for accessor through the random mechanism.

    :param observed_length: length of the DNA string in a vertex.
    :type observed_length: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :param random_seed: random seed for creating shuffles.
    :type random_seed: int

    :param verbose: need to print log.
    :type verbose: bool

    :return: shuffles for accessor.
    :rtype: numpy.ndarray

    ..note::
        The shuffle strategy disrupts the bit-to-nucleotide mapping order, so it can be used as an encryption policy.

        Value 0 ~ 3 in each line shuffle only describes the relationship between progressive order and position.
        The original position is [0, 1, 2, 3].
        If the follow-up nucleotides in a vertex are ["A", "G"] when the line shuffle is [3, 2, 1, 0],
        This line shuffle can be regarded as [1, -, 0, -] for ["A", -, "G", -].

        Ths shuffle will not disclose the topology information of accessor.

    ..example::
        >>> from dsw import create_random_shuffles
        >>> create_random_shuffles(observed_length=2, random_seed=2021)
        array([[3, 2, 1, 0],
               [2, 3, 1, 0],
               [3, 1, 0, 2],
               [0, 3, 1, 2],
               [3, 2, 0, 1],
               [1, 0, 3, 2],
               [0, 3, 1, 2],
               [2, 0, 1, 3],
               [2, 3, 0, 1],
               [1, 0, 3, 2],
               [2, 0, 1, 3],
               [0, 1, 3, 2],
               [2, 3, 1, 0],
               [2, 0, 3, 1],
               [0, 1, 3, 2],
               [0, 3, 2, 1]])
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    shuffles = zeros(shape=(4 ** observed_length, len(nucleotides)), dtype=int)
    shuffles[:, 1] = 1
    shuffles[:, 2] = 2
    shuffles[:, 3] = 3

    random.seed(random_seed)

    monitor = Monitor()
    for index in range(4 ** observed_length):
        card = shuffles[index]
        random.shuffle(card)
        shuffles[index] = card

        if verbose:
            monitor.output(index + 1, 4 ** observed_length)

    random.seed(None)

    return shuffles


def remove_random_edges(accessor, remove_number, random_seed=None, verbose=False):
    """
    Remove directed edges in accessor through the random mechanism.

    :param accessor: accessor of the coding algorithm.
    :type accessor: numpy.ndarray

    :param remove_number: number of directed edge removed in the accessor.
    :type remove_number: int

    :param random_seed: random seed for removing directed edges in accessor.
    :type random_seed: int

    :param verbose: need to print log.
    :type verbose: bool

    :return: records of removed directed edges.
    :rtype: list

    ..example::
        >>> from dsw import remove_random_edges, get_complete_accessor
        >>> accessor = get_complete_accessor(observed_length=2)
        >>> remove_random_edges(accessor, remove_number=2, random_seed=2021)
        [(4, 1), (9, 0)]
    """
    outs = zeros(shape=(len(accessor),))
    for index, vertex in enumerate(accessor):
        used_indices = where(vertex >= 0)[0]
        outs[index] = len(used_indices)

    usage_3, usage_4 = len(outs[outs == 3]), len(outs[outs == 4])

    if remove_number > usage_3 + 2 * usage_4:
        raise ValueError("More edges are required to be deleted than can be deleted!")

    random.seed(random_seed)

    records, monitor = [], Monitor()

    while True:
        used_index = random.choice(where(outs >= 3)[0])
        remove_choice = random.choice(where(accessor[used_index] >= 0)[0])
        records.append((used_index, remove_choice))
        outs[used_index] -= 1

        if verbose:
            monitor.output(len(records), remove_number)

        if len(records) == remove_number:
            break

    random.seed(None)

    return records
