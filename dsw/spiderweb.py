from collections import Counter
from itertools import product
from networkx import DiGraph, find_cycle
from numpy import zeros, ones, array, random, log, sum, max, argmax, argsort, unique, intersect1d, where

from dsw.operation import Monitor, calculus_addition, calculus_multiplication, calculus_division
from dsw.operation import bit_to_number, number_to_bit, number_to_dna, dna_to_number
from dsw.graphized import obtain_vertices, obtain_formers, obtain_latters, path_matching, calculate_intersection_score


def encode(binary_message, accessor, start_index,
           is_faster=False, vt_length=0, shuffles=None, need_path=False, verbose=False):
    """
    Encode a bit array by the specific accessor.

    :param binary_message: binary message.
    :type binary_message: numpy.ndarray

    :param accessor: accessor of the coding algorithm.
    :type accessor: numpy.ndarray

    :param start_index: virtual vertex to start encoding.
    :type start_index: int

    :param is_faster: encode in a faster way.
    :type is_faster: bool

    :param vt_length: length of Varshamov-Tenengolts code for error-correction.
    :type vt_length: int or None

    :param shuffles: shuffle relationships for bit-nucleotide mapping.
    :type: numpy.ndarray

    :param need_path: need to record the restricted state transition path.
    :type need_path: bool

    :param verbose: need to print log.
    :type verbose: bool

    :return: DNA sequence encoded by this graph (and VT check sequence if required).
    :rtype: str or (str, str)

    .. note::
        If the parameter "is_faster" is set as True, the out-degree of coding digraph cannot contains 3.


    Example
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
        >>> encode(accessor=accessor, binary_message=binary_message, start_index=1, vt_length=5)
        ('TCTCTCT', 'TAAGC')
        >>> shuffles = array([[3, 2, 1, 0], [2, 3, 1, 0], [3, 1, 0, 2], [0, 3, 1, 2], \
                              [3, 2, 0, 1], [1, 0, 3, 2], [0, 3, 1, 2], [2, 0, 1, 3], \
                              [2, 3, 0, 1], [1, 0, 3, 2], [2, 0, 1, 3], [0, 1, 3, 2], \
                              [2, 3, 1, 0], [2, 0, 3, 1], [0, 1, 3, 2], [0, 3, 2, 1]])
        >>> encode(accessor=accessor, binary_message=binary_message, shuffles=shuffles, start_index=1)
        'AGAGAGA'
    """
    monitor, record_path, vertex_index, dna_sequence, nucleotides = Monitor(), [], start_index, "", "ACGT"

    if not is_faster:
        quotient = bit_to_number(binary_message, verbose=verbose)
        total_state = len(quotient)  # number of symbol.

        while quotient != "0":
            used_indices = where(accessor[vertex_index] >= 0)[0]

            if len(used_indices) > 1:  # current vertex contains information.
                quotient, remainder = calculus_division(number=quotient, base=str(len(used_indices)))
                remainder = int(remainder)

                if shuffles is not None:  # shuffle remainder based on the inputted shuffles.
                    remainder = argsort(shuffles[vertex_index, used_indices])[remainder]

                value = used_indices[remainder]

                if need_path:
                    record_path.append([vertex_index, 1])

            elif len(used_indices) == 1:  # current vertex does not contain information.
                value = used_indices[0]

                if need_path:
                    record_path.append([vertex_index, 0])

            else:  # current vertex is wrong.
                raise ValueError("Current vertex doesn't have an out-degree, "
                                 + "the accessor or the start vertex is wrong!")

            nucleotide, vertex_index = nucleotides[value], accessor[vertex_index][value]

            dna_sequence += nucleotide

            if verbose:
                if quotient != "0":
                    monitor(total_state - len(quotient), total_state)
                else:
                    monitor(total_state, total_state)

    else:
        location = 0
        while location < len(binary_message):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            radix = len(used_indices)

            if radix == 4:  # current vertex contains information.
                remainder = binary_message[location] * 2 + binary_message[location + 1]

                if shuffles is not None:  # shuffle remainder based on the inputted shuffles.
                    remainder = argsort(shuffles[vertex_index, used_indices])[remainder]

                value = used_indices[remainder]
                location += 2
            elif radix == 2:
                remainder = binary_message[location]

                if shuffles is not None:  # shuffle remainder based on the inputted shuffles.
                    remainder = argsort(shuffles[vertex_index, used_indices])[remainder]

                value = used_indices[remainder]
                location += 1
            elif radix == 1:  # current vertex does not contain information.
                value = used_indices[0]
            elif radix == 3:
                raise ValueError("Not implementation!")
            else:  # current vertex is wrong.
                raise ValueError("Current vertex doesn't have an out-degree, "
                                 + "the accessor or the start vertex is wrong!")

            nucleotide, vertex_index = nucleotides[value], accessor[vertex_index][value]
            dna_sequence += nucleotides[value]

            if need_path:
                record_path.append([vertex_index, int(radix > 1)])

            if verbose:
                monitor(location + 1, len(binary_message))

    if need_path:
        record_path = array(record_path)

    if vt_length > 0:
        vt_check = set_vt(dna_sequence=dna_sequence, vt_length=vt_length)
        if need_path:
            return dna_sequence, set_vt(dna_sequence=dna_sequence, vt_length=vt_length), record_path
        else:
            return dna_sequence, vt_check
    else:
        if need_path:
            return dna_sequence, record_path
        else:
            return dna_sequence


def decode(dna_sequence, bit_length, accessor, start_index,
           is_faster=False, vt_check=None, shuffles=None, verbose=False):
    """
    Decode a DNA sequence by the specific accessor.

    :param dna_sequence: DNA sequence encoded by this graph.
    :type dna_sequence: str

    :param bit_length: length of the bit array.
    :type bit_length: int

    :param accessor: accessor of the coding algorithm (consistent with the encoding process).
    :type accessor: numpy.ndarray

    :param start_index: virtual vertex to start decoding (consistent with the encoding process).
    :type start_index: int

    :param is_faster: encode in a faster way.
    :type is_faster: bool

    :param vt_check: check sequence of Varshamov-Tenengolts code.
    :type vt_check: str or None

    :param shuffles: shuffle relationships for bit-nucleotide mapping.
    :type shuffles: numpy.ndarray

    :param verbose: need to print log.
    :type verbose: bool

    :return: binary message decoded by this graph.
    :rtype: numpy.ndarray

    :raise ValueError: if one or more errors are found.

    .. note::
        If the parameter "is_faster" is set as True, the out-degree of coding digraph cannot contains 3.

    Example
        >>> from numpy import array
        >>> from dsw import decode
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> dna_sequence = "TCTCTCT"
        >>> decode(accessor=accessor, dna_sequence=dna_sequence, start_index=1, bit_length=8)
        array([0, 1, 0, 1, 0, 1, 0, 1])
        >>> vt_check = "TAAGC"
        >>> decode(accessor=accessor, dna_sequence=dna_sequence, start_index=1, vt_check=vt_check, bit_length=8)
        array([0, 1, 0, 1, 0, 1, 0, 1])
        >>> shuffles = array([[3, 2, 1, 0], [2, 3, 1, 0], [3, 1, 0, 2], [0, 3, 1, 2], \
                              [3, 2, 0, 1], [1, 0, 3, 2], [0, 3, 1, 2], [2, 0, 1, 3], \
                              [2, 3, 0, 1], [1, 0, 3, 2], [2, 0, 1, 3], [0, 1, 3, 2], \
                              [2, 3, 1, 0], [2, 0, 3, 1], [0, 1, 3, 2], [0, 3, 2, 1]])
        >>> decode(accessor=accessor, dna_sequence="TCTCTCT", start_index=1, shuffles=shuffles, bit_length=8)
        array([0, 0, 0, 0, 0, 0, 0, 0])
    """
    vertex_index, nucleotides, monitor = start_index, "ACGT", Monitor()

    if vt_check is not None:
        if vt_check != set_vt(dna_sequence=dna_sequence, vt_length=len(vt_check)):
            raise ValueError("At least one error is found in this DNA sequence!")

    if not is_faster:
        quotient, saved_values = "0", []

        for location, nucleotide in enumerate(dna_sequence):
            used_indices = where(accessor[vertex_index] >= 0)[0]

            if len(used_indices) > 1:  # current vertex contains information.
                used_nucleotides = [nucleotides[used_index] for used_index in used_indices]

                if nucleotide in used_nucleotides:  # check whether the DNA sequence is right currently.
                    remainder = used_nucleotides.index(nucleotide)
                else:
                    raise ValueError("At least one error is found in this DNA sequence!")

                if shuffles is not None:  # shuffle remainder based on the inputted shuffles.
                    remainder = where(argsort(shuffles[vertex_index, used_indices]) == remainder)[0][0]

                saved_values.append((len(used_indices), remainder))
                vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]

            elif len(used_indices) == 1:  # current vertex does not contain information.
                used_nucleotide = nucleotides[used_indices[0]]
                if nucleotide == used_nucleotide:
                    vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                else:
                    raise ValueError("At least one error is found in this DNA sequence!")

            else:  # current vertex is wrong.
                raise ValueError("Current vertex doesn't have an out-degree, "
                                 + "the accessor, the start vertex, or DNA sequence is wrong!")

            if verbose:
                monitor(location + 1, len(dna_sequence))

        for location, (out_degree, number) in enumerate(saved_values[::-1]):
            quotient = calculus_multiplication(number=quotient, base=str(out_degree))
            quotient = calculus_addition(number=quotient, base=str(number))

        binary_message = array(number_to_bit(decimal_number=quotient, bit_length=bit_length), dtype=int)

    else:
        message_location, binary_message = 0, zeros(shape=(bit_length,), dtype=int)

        for location, nucleotide in enumerate(dna_sequence):
            used_indices = where(accessor[vertex_index] >= 0)[0]
            radix = len(used_indices)

            used_nucleotides = [nucleotides[used_index] for used_index in used_indices]
            if nucleotide in used_nucleotides:  # check whether the DNA sequence is right currently.
                remainder = used_nucleotides.index(nucleotide)
            else:
                raise ValueError("At least one error is found in this DNA sequence!")

            if shuffles is not None:  # shuffle remainder based on the inputted shuffles.
                remainder = where(argsort(shuffles[vertex_index, used_indices]) == remainder)[0][0]

            vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]

            if radix == 4:
                binary_message[message_location] = remainder // 2
                binary_message[message_location + 1] = remainder % 2
                message_location += 2
            elif radix == 2:
                binary_message[message_location] = remainder % 2
                message_location += 1
            elif radix == 1:
                pass  # do nothing.
            elif radix == 3:
                raise ValueError("Not implementation!")
            else:
                raise ValueError("Current vertex doesn't have an out-degree, "
                                 + "the accessor, the start vertex, or DNA sequence is wrong!")

    return binary_message


def set_vt(dna_sequence, vt_length):
    """
    Set Varshamov-Tenengolts-based path check string ('salt-protected') from DNA (payload) sequence.

    :param dna_sequence: DNA sequence encoded through SPIDER-WEB.
    :type dna_sequence: str

    :param vt_length: length of DNA sequence (path check).
    :type vt_length: int or None

    :return: path check DNA sequence with required length.

    Example
        >>> from dsw import set_vt
        >>> dna_sequence = "TCTCTCT"
        >>> vt_length = 5
        >>> set_vt(dna_sequence=dna_sequence, vt_length=vt_length)
        'TAAGC'

    .. note::
        Reference [1] Rom R. Varshamov and Grigory M. Tenengolts (1965) Avtomat. i Telemekh

        Reference [2] Grigory Tenengolts (1984) IEEE Transactions on Information Theory

        Reference [3] William H. Press et al. (2020) Proceedings of the National Academy of Sciences

        Reference [4] A. Xavier Kohll et al. (2020) Chemical Communications
    """
    nucleotides = "ACGT"

    values = array([nucleotides.index(nucleotide) for nucleotide in dna_sequence])
    vt_value = sum(where((values[1:] - values[:-1]) > 0)[0]) % (len(nucleotides) ** (vt_length - 1))
    vt_flag = sum(values) % len(nucleotides)

    return nucleotides[vt_flag] + number_to_dna(decimal_number=int(vt_value), dna_length=vt_length - 1)


def repair_dna(dna_sequence, accessor, start_index, observed_length, vt_check=None, has_indel=False, heap_size=1e3):
    """
    Repair the DNA sequence containing one (or more) errors.

    :param dna_sequence: DNA sequence waiting for recovery.
    :type dna_sequence: str

    :param accessor: accessor of the coding algorithm.
    :type accessor: numpy.ndarray

    :param start_index: virtual vertex to start encoding.
    :type start_index: int

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :param vt_check: check sequence of Varshamov-Tenengolts code.
    :type vt_check: str or None

    :param has_indel: consider insertion and/or deletion error.
    :type has_indel: bool

    :param heap_size: maximum heap size.
    :type heap_size: int

    :return: repaired DNA sequence set and additional information.
    :rtype: (list, (bool, bool, int))

    Example
        >>> from numpy import array
        >>> from dsw import repair_dna
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> dna_sequence = "TCTCTATCTCTC"  # "TCTCTCTCTCTC" is original DNA sequence
        >>> repair_dna(dna_sequence=dna_sequence, accessor=accessor, start_index=1, observed_length=2, has_indel=True)
        (['TCTCTCTCTCTC', 'TCTCTGTCTCTC'], (1, False, 2, 14))
        >>> vt_check = "AACGC"  # check list of Varshamov-Tenengolts code.
        >>> repair_dna(dna_sequence=dna_sequence, accessor=accessor, start_index=1, observed_length=2, \
                       vt_check=vt_check, has_indel=True)
        (['TCTCTCTCTCTC'], (1, True, 2, 14))
    """
    nucleotides = "ACGT"

    location, vertex_index, index_queue = 0, start_index, -ones(shape=(len(dna_sequence),), dtype=int)
    split_sequences, chuck_sequences, index_markers = [""], [], []
    detected_count, chuck_flag, visited_times = 0, False, 0

    while location < len(dna_sequence):
        used_indices, nucleotide = where(accessor[vertex_index] >= 0)[0], dna_sequence[location]
        if dna_sequence[location] in [nucleotides[used_index] for used_index in used_indices]:
            split_sequences[-1] += nucleotide
            vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
            index_queue[location] = vertex_index
            visited_times += 1
            location += 1
        elif len(split_sequences[-1]) > 0:
            detected_count += 1
            split_sequences[-1] = split_sequences[-1][: - observed_length + 1]
            vertex_index = dna_to_number(dna_sequence[location + 1: location + observed_length + 1], is_string=False)
            split_sequences.append(nucleotides[vertex_index % 4])
            index_markers.append(index_queue[location - observed_length: location])
            chuck_sequences.append(dna_sequence[location - observed_length + 1: location + observed_length])
            location += observed_length + 1

    repaired_fragment_set = [set() for _ in range(len(index_markers))]
    for index, (chuck_sequence, index_marker) in enumerate(zip(chuck_sequences, index_markers)):
        for recall, vertex_index in enumerate(index_marker[::-1]):
            record, times = path_matching(dna_sequence=chuck_sequence, accessor=accessor, has_indel=has_indel,
                                          previous_index=vertex_index, occur_location=observed_length - recall - 1)
            visited_times += times
            for _, fragment in record:
                if dna_sequence not in repaired_fragment_set[index]:
                    repaired_fragment_set[index].add(fragment)
        repaired_fragment_set[index] = list(repaired_fragment_set[index])

    repaired_results, count = set(), 1
    for fragments in repaired_fragment_set:
        count *= len(fragments)

    if count == 0 or count > heap_size:
        if vt_check is not None:
            if vt_check == set_vt(dna_sequence=dna_sequence, vt_length=len(vt_check)):
                return [dna_sequence], (0, False, 0, visited_times)
            else:
                return [], (0, True, 0, visited_times)
        else:
            return [dna_sequence], (0, False, 0, visited_times)

    for fragments in product(*repaired_fragment_set):
        repaired_dna_sequence = ""
        for index in range(len(split_sequences) - 1):
            repaired_dna_sequence += split_sequences[index] + fragments[index]
        repaired_dna_sequence += split_sequences[-1]
        if vt_check is not None:
            if vt_check == set_vt(dna_sequence=repaired_dna_sequence, vt_length=len(vt_check)):
                repaired_results.add(repaired_dna_sequence)
            else:
                chuck_flag = True
        else:
            repaired_results.add(repaired_dna_sequence)

    return sorted(list(repaired_results)), (detected_count, chuck_flag, count, visited_times)


def find_vertices(observed_length, bio_filter, verbose=False):
    """
    Find valid vertices based on the given the biochemical constraints.

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :param bio_filter: screening operation for identifying the valid DNA sequence (required the given constraints).
    :type bio_filter: dsw.biofilter.DefaultBioFilter

    :param verbose: need to print log.
    :type verbose: bool

    :return: available vertices.
    :rtype: numpy.ndarray

    :raise ValueError: if no valid vertices are available under the required biochemical constraints.

    Example
        >>> from dsw import LocalBioFilter, find_vertices
        >>> bio_filter = LocalBioFilter(observed_length=2, max_homopolymer_runs=2, gc_range=[0.5, 0.5])
        >>> # "1" refers to the available index of vertex, otherwise "0" refers to unavailable index of vertex.
        >>> find_vertices(observed_length=2, bio_filter=bio_filter).astype(int)
        array([0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0])

    .. note::
        Reference [1] Florent Capelli and Yann Strozecki (2019) Discrete Applied Mathematics
    """
    nucleotides = "ACGT"

    vertices, monitor = zeros(shape=(int(len(nucleotides) ** observed_length),), dtype=bool), Monitor()

    if verbose:
        print("Find valid vertices in this observed length of DNA sequence.")

    for vertex_index in range(len(vertices)):
        dna_sequence = number_to_dna(decimal_number=vertex_index, dna_length=observed_length)
        vertices[vertex_index] = bio_filter.valid(dna_sequence=dna_sequence)

        if verbose:
            monitor(vertex_index + 1, len(vertices), extra={"valid": sum(vertices[: vertex_index + 1])})

    valid_rate = sum(vertices) / len(vertices)

    if valid_rate == 0:
        raise ValueError("No vertex is collected!")

    if verbose:
        print(str(round(valid_rate * 100, 2)) + "% (" + str(sum(vertices)) + ") valid vertices are collected.")

    return vertices


def connect_valid_graph(observed_length, vertices, verbose=False):
    """
    Connect a valid graph by valid vertices.

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray

    :param verbose: need to print log.
    :type verbose: bool

    :return: accessor of the valid graph.
    :rtype: numpy.ndarray

    :raise ValueError: if no valid vertices are available.

    Example
        >>> from numpy import array
        >>> from dsw import connect_valid_graph
        >>> vertices = array([0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0])
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

    .. note::
        Reference [1] Nicolaas Govert de Bruijn (1946) Indagationes Mathematicae
    """
    nucleotides = "ACGT"

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
                monitor(vertex_index + 1, len(vertices))

        if verbose:
            print("Valid graph is created.")

        return accessor
    else:
        raise ValueError("No collected vertex!")


def connect_coding_graph(observed_length, vertices, threshold, verbose=False):
    """
    Connect a coding algorithm by valid vertices and the threshold for minimum out-degree.

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray

    :param threshold: threshold for minimum out-degree.
    :type threshold: int

    :param verbose: need to print log.
    :type verbose: bool

    :return: coding vertices and coding accessor.
    :rtype: (numpy.ndarray, numpy.ndarray)

    :raise ValueError: if no coding graph are created because of the vertex set and/or the trimming requirement.

    Example
        >>> from numpy import array
        >>> from dsw import connect_coding_graph
        >>> vertices = array([0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0])
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
        >>> vertices.astype(int)
        array([0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0])

    .. note::
        Reference [1] Nicolaas Govert de Bruijn (1946) Indagationes Mathematicae
    """
    times, nucleotides = 1, "ACGT"

    while True:
        if verbose:
            print("Check the vertex collection requirement in round " + str(times) + ".")

        new_vertices, monitor = zeros(shape=(int(len(nucleotides) ** observed_length),), dtype=bool), Monitor()
        saved_indices = where(vertices != 0)[0]
        for current, vertex_index in enumerate(saved_indices):
            latter_indices = obtain_latters(current=vertex_index, observed_length=observed_length)
            new_vertices[vertex_index] = sum(vertices[latter_indices]) >= threshold

            if verbose:
                monitor(current + 1, len(saved_indices))

        changed = sum(vertices) - sum(new_vertices)

        if verbose:
            print(str(round(sum(new_vertices) / len(vertices) * 100, 2)) + "% (" + str(sum(new_vertices)) + ") "
                  + "valid vertices are saved.")

        if sum(new_vertices) < 1:
            raise ValueError("No coding graph is created!")

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
                monitor(vertex_index + 1, len(vertices))

        if threshold == 1:
            while True:
                vertices = obtain_vertices(accessor)
                graph = DiGraph()
                for former_index, latter_indices in enumerate(accessor):
                    for latter_index in latter_indices:
                        if latter_index >= 0:
                            graph.add_edge(u_of_edge=former_index, v_of_edge=latter_index)
                useless_vertices, cycle = [], find_cycle(graph)
                for former_index, latter_index in cycle:
                    if len(where(accessor[former_index] >= 0)[0]) == 1:
                        useless_vertices.append(former_index)
                if len(useless_vertices) == len(cycle):
                    for useless_vertex in useless_vertices:
                        accessor[useless_vertex] = -1
                        pairs = [(i, useless_vertex) for i in obtain_formers(useless_vertex, 10)]
                        while len(pairs) > 0:
                            new_pairs = []
                            for former_index, latter_index in pairs:
                                previous = len(where(accessor[former_index] >= 0)[0])
                                accessor[former_index, latter_index % 4] = -1
                                current = len(where(accessor[former_index] >= 0)[0])
                                if previous > current == 0:
                                    new_pairs += [(i, former_index) for i in obtain_formers(former_index, 10)]
                            pairs = new_pairs
                else:
                    break

        if verbose:
            print("The coding graph is created.")

        return vertices, accessor
    else:
        raise ValueError("The coding graph cannot be created!")


def remove_nasty_arc(accessor, latter_map, iteration=0, has_insertion=True, has_deletion=True, verbose=False):
    """
    Remove the nasty arc.

    :param accessor: accessor of the coding algorithm.
    :type accessor: numpy.ndarray

    :param latter_map: latter map of the coding algorithm.
    :type latter_map: dict

    :param iteration: current round if required.
    :type iteration: int

    :param has_insertion: need to repair insertion error.
    :type has_insertion: bool

    :param has_deletion: need to repair deletion error.
    :type has_deletion: bool

    :param verbose: need to print log.
    :type verbose: bool

    :return: adjusted accessor, adjusted latter map, removed arc, and maximum intersection score.
    :rtype: (numpy.ndarray, dict, tuple, int)

    Example
        >>> from numpy import array
        >>> from dsw import accessor_to_latter_map, remove_nasty_arc
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> latter_map = accessor_to_latter_map(accessor)
        >>> result = remove_nasty_arc(accessor=accessor, latter_map=latter_map, iteration=0, verbose=False, \
                                      has_insertion=True, has_deletion=True)
        >>> result[0]
        array([[-1, -1, -1, -1],
               [-1, -1, -1,  7],
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
        >>> result[1]
        {1: [7], 2: [8, 11], 4: [1, 2], 7: [13, 14], 8: [1, 2], 11: [13, 14], 13: [4, 7], 14: [8, 11]}
        >>> result[2]
        (1, 4)
        >>> result[3]
        [16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16]

    .. note::
        The graph information contained in "accessor" and "latter_map" is consistent.

        It is a gift for the follow-up investigation.
        That is, removing arc to improve the capability of the probabilistic error correction.
    """
    nucleotides = "ACGT"

    observed_length = int(log(len(accessor)) / log(len(nucleotides)))

    if verbose:
        if iteration > 0:
            print("Calculate the intersection score for each remained arc in " + str(iteration) + " round(s).")
        else:
            print("Calculate the intersection score for each remained arc.")

    scores = calculate_intersection_score(latter_map=latter_map, has_insertion=has_insertion, has_deletion=has_deletion,
                                          observed_length=observed_length, verbose=verbose)

    vertex_indices = unique(where(scores == max(scores))[0])

    former = intersect1d(obtain_vertices(accessor=accessor), vertex_indices)[0]
    former_value, latter_value = former % len(nucleotides), argmax(scores[former])
    latter = int((former * len(nucleotides) + latter_value) % (len(nucleotides) ** observed_length))

    accessor[former, latter_value] = -1
    del latter_map[former][latter_map[former].index(latter)]
    if len(latter_map[former]) == 0:
        del latter_map[former]

    scores = scores.reshape(-1)
    scores = scores[scores > 0].tolist()
    score_record = array(list(Counter(scores).items())).T
    score_record = score_record[:, argsort(score_record[0])[::-1]]

    if verbose:
        print("Current scores are:")
        print(score_record)
        print("Remove arc " + number_to_dna(int(former), dna_length=observed_length) + " -> "
              + number_to_dna(int(latter), dna_length=observed_length)
              + " with the maximum intersection score " + str(max(scores)) + ".")

    return accessor, latter_map, (former, latter), scores


def create_random_shuffles(observed_length, random_seed=None, verbose=False):
    """
    Create the shuffles for accessor through the random mechanism.

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :param random_seed: random seed for creating shuffles.
    :type random_seed: int

    :param verbose: need to print log.
    :type verbose: bool

    :return: shuffles for accessor.
    :rtype: numpy.ndarray

    Example
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

    .. note::
        The mapping shuffle strategy disrupts the digit-to-nucleotide mapping order,
        so it can be used as an privacy protection mechanism.

        Value 0 ~ 3 in each line shuffle only describes the relationship between progressive order and position.
        The original position is [0, 1, 2, 3].
        If the follow-up nucleotides in a vertex are ["A", "G"] when the line shuffle is [3, 2, 1, 0],
        This line shuffle can be regarded as [1, -, 0, -] for ["A", -, "G", -].

        Ths shuffle will not disclose the topology information of accessor.
    """
    nucleotides = "ACGT"

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
            monitor(index + 1, 4 ** observed_length)

    random.seed(None)

    return shuffles
