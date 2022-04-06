from copy import deepcopy
from numpy import zeros, ones, array, random, sum, argsort, where

from dsw.operation import Monitor, calculus_addition, calculus_multiplication, calculus_division
from dsw.operation import bit_to_number, number_to_bit, number_to_dna
from dsw.graphized import obtain_latters, path_matching


def encode(binary_message, accessor, start_index, nucleotides=None,
           vt_length=0, shuffles=None, verbose=False):
    """
    Encode a bit array by the specific accessor.

    :param binary_message: binary message.
    :type binary_message: numpy.ndarray

    :param accessor: accessor of the coding algorithm.
    :type accessor: numpy.ndarray

    :param start_index: virtual vertex to start encoding.
    :type start_index: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list or None

    :param vt_length: length of Varshamov-Tenengolts code for error-correction.
    :type vt_length: int or None

    :param shuffles: shuffle relationships for bit-nucleotide mapping.
    :type: numpy.ndarray

    :param verbose: need to print log.
    :type verbose: bool

    :return: DNA string encoded by this graph (and VT check string if required).
    :rtype: str or (str, str)

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
        ('TCTCTCT', 'TAATA')
        >>> shuffles = array([[3, 2, 1, 0], [2, 3, 1, 0], [3, 1, 0, 2], [0, 3, 1, 2], \
                              [3, 2, 0, 1], [1, 0, 3, 2], [0, 3, 1, 2], [2, 0, 1, 3], \
                              [2, 3, 0, 1], [1, 0, 3, 2], [2, 0, 1, 3], [0, 1, 3, 2], \
                              [2, 3, 1, 0], [2, 0, 3, 1], [0, 1, 3, 2], [0, 3, 2, 1]])
        >>> encode(accessor=accessor, binary_message=binary_message, shuffles=shuffles, start_index=1)
        'AGAGAGA'
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    quotient, dna_string, vertex_index, monitor = bit_to_number(binary_message), "", start_index, Monitor()
    total_state = len(quotient)  # number of symbol.

    while quotient != "0":
        used_indices = where(accessor[vertex_index] >= 0)[0]

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

        nucleotide, vertex_index = nucleotides[value], accessor[vertex_index][value]

        dna_string += nucleotide

        if verbose:
            if quotient != "0":
                monitor.output(total_state - len(quotient), total_state)
            else:
                monitor.output(total_state, total_state)

    if vt_length > 0:
        vt_check = set_vt(dna_string=dna_string, vt_length=vt_length, nucleotides=nucleotides)
        return dna_string, vt_check
    else:
        return dna_string


def decode(dna_string, bit_length, accessor, start_index, nucleotides=None,
           vt_check=None, shuffles=None, verbose=False):
    """
    Decode a DNA string by the specific accessor.

    :param dna_string: DNA string encoded by this graph.
    :type dna_string: str

    :param bit_length: length of the bit array.
    :type bit_length: int

    :param accessor: accessor of the coding algorithm (consistent with the encoding process).
    :type accessor: numpy.ndarray

    :param start_index: virtual vertex to start decoding (consistent with the encoding process).
    :type start_index: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list or None

    :param vt_check: check string of Varshamov-Tenengolts code.
    :type vt_check: str or None

    :param shuffles: shuffle relationships for bit-nucleotide mapping.
    :type shuffles: numpy.ndarray

    :param verbose: need to print log.
    :type verbose: bool

    :return: binary message decoded by this graph.
    :rtype: numpy.ndarray

    :raise ValueError: if one or more errors are found.

    Example
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
        >>> vt_check = "TAATA"
        >>> decode(accessor=accessor, dna_string=dna_string, start_index=1, vt_check=vt_check, bit_length=8)
        array([0, 1, 0, 1, 0, 1, 0, 1])
        >>> shuffles = array([[3, 2, 1, 0], [2, 3, 1, 0], [3, 1, 0, 2], [0, 3, 1, 2], \
                              [3, 2, 0, 1], [1, 0, 3, 2], [0, 3, 1, 2], [2, 0, 1, 3], \
                              [2, 3, 0, 1], [1, 0, 3, 2], [2, 0, 1, 3], [0, 1, 3, 2], \
                              [2, 3, 1, 0], [2, 0, 3, 1], [0, 1, 3, 2], [0, 3, 2, 1]])
        >>> decode(accessor=accessor, dna_string="TCTCTCT", start_index=1, shuffles=shuffles, bit_length=8)
        array([0, 0, 0, 0, 0, 0, 0, 0])
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    quotient, saved_values, vertex_index, monitor = "0", [], start_index, Monitor()

    for location, nucleotide in enumerate(dna_string):
        used_indices = where(accessor[vertex_index] >= 0)[0]

        if len(used_indices) > 1:  # current vertex contains information.
            used_nucleotides = [nucleotides[used_index] for used_index in used_indices]

            if nucleotide in used_nucleotides:  # check whether the DNA string is right currently.
                remainder = used_nucleotides.index(nucleotide)
            else:
                raise ValueError("At least one error is found in this DNA string!")

            if shuffles is not None:  # shuffle remainder based on the inputted shuffles.
                remainder = where(argsort(shuffles[vertex_index, used_indices]) == remainder)[0][0]

            saved_values.append((len(used_indices), remainder))
            vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]

        elif len(used_indices) == 1:  # current vertex does not contain information.
            used_nucleotide = nucleotides[used_indices[0]]
            if nucleotide == used_nucleotide:
                vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
            else:
                raise ValueError("At least one error is found in this DNA string!")

        else:  # current vertex is wrong.
            raise ValueError("Current vertex doesn't have an out-degree, "
                             + "the accessor, the start vertex, or DNA string is wrong!")

        if verbose:
            monitor.output(location + 1, len(dna_string))

    if vt_check is not None:
        if vt_check != set_vt(dna_string=dna_string, vt_length=len(vt_check), nucleotides=nucleotides):
            raise ValueError("At least one error is found in this DNA string!")

    for location, (out_degree, number) in enumerate(saved_values[::-1]):
        quotient = calculus_multiplication(number=quotient, base=str(out_degree))
        quotient = calculus_addition(number=quotient, base=str(number))

    return array(number_to_bit(decimal_number=quotient, bit_length=bit_length), dtype=int)


def set_vt(dna_string, vt_length, nucleotides=None):
    """
    Set Varshamov-Tenengolts-based path check string ('salt-protected') from DNA (payload) string.

    :param dna_string: DNA string encoded through SPIDER-WEB.
    :type dna_string: str

    :param vt_length: length of DNA string (path check).
    :type vt_length: int or None

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list or None

    :return: path check DNA string with required length.

    Example
        >>> from dsw import set_vt
        >>> dna_string = "TCTCTCT"
        >>> vt_length = 5
        >>> set_vt(dna_string=dna_string, vt_length=vt_length)
        'TAATA'

    .. note::
        Reference [1] Rom R. Varshamov and Grigory M. Tenengolts (1965) Avtomat. i Telemekh

        Reference [2] Grigory Tenengolts (1984) IEEE Transactions on Information Theory

        Reference [3] William H. Press et al. (2020) Proceedings of the National Academy of Sciences

        Reference [4] A. Xavier Kohll et al. (2020) Chemical Communications
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    vt_flag, vt_value = 0, 0
    for location, nucleotide in enumerate(dna_string):
        value = nucleotides.index(nucleotide)
        if location > 0 and value >= nucleotides.index(dna_string[location - 1]):
            vt_value += location
        vt_flag += value
    vt_value %= 4 ** (vt_length - 1)
    vt_flag %= 4

    former = nucleotides[vt_flag % len(nucleotides)]
    latter = number_to_dna(decimal_number=vt_value % len(nucleotides) ** (vt_length - 1), dna_length=vt_length - 1)

    return former + latter


def repair_dna(dna_string, accessor, start_index, observed_length, check_iterations=1, heap_size=1e6, vt_check=None,
               has_insertion=False, has_deletion=False, nucleotides=None, verbose=False):
    """
    Repair the DNA string containing one (or more) errors.

    :param dna_string: DNA string waiting for recovery.
    :type dna_string: str

    :param accessor: accessor of the coding algorithm.
    :type accessor: numpy.ndarray

    :param start_index: virtual vertex to start encoding.
    :type start_index: int

    :param observed_length: length of the DNA string in a vertex.
    :type observed_length: int

    :param check_iterations: repair check iterations, which greater than or equal to the estimated number of errors.
    :type check_iterations: int

    :param heap_size: available set size of repaired DNA strings.
    :type heap_size: int or None

    :param vt_check: check string of Varshamov-Tenengolts code.
    :type vt_check: str or None

    :param has_insertion: consider insertion error.
    :type has_insertion: bool

    :param has_deletion: consider deletion error.
    :type has_deletion: bool

    :param verbose: need to print log.
    :type verbose: bool

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :return: repaired DNA string set and additional information (includes detect flags and initial repaired strings).
    :rtype: (list, (bool, bool, int))

    Example
        >>> from numpy import array
        >>> from dsw import repair_dna
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> dna_string = "TCTCTATCTCTC"  # "TCTCTCTCTCTC" is original DNA string
        >>> repair_dna(dna_string=dna_string, accessor=accessor, start_index=1, observed_length=2, \
                       check_iterations=1, has_insertion=True, has_deletion=True, verbose=False)
        (['TCTCTCTCTCTC', 'TCTCTGTCTCTC'], (True, False, 2))
        >>> vt_check = "AACTG"  # check list of Varshamov-Tenengolts code.
        >>> repair_dna(dna_string=dna_string, accessor=accessor, start_index=1, observed_length=2, \
                       check_iterations=1, vt_check=vt_check, has_insertion=True, has_deletion=True, verbose=False)
        (['TCTCTCTCTCTC'], (True, True, 2))
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    original_dna_string = deepcopy(dna_string)
    dna_strings, repaired_dna_strings, detected_flag_1, detected_flag_2 = [dna_string], [], False, False
    for _ in range(check_iterations + 1):
        new_dna_strings = []
        for detected_dna_string in dna_strings:
            found_wrong, vertex_index, index_queue = False, start_index, [start_index]
            for location, nucleotide in enumerate(detected_dna_string):
                used_indices = where(accessor[vertex_index] >= 0)[0]
                if nucleotide in [nucleotides[used_index] for used_index in used_indices]:
                    vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                    index_queue.append(vertex_index)
                else:
                    if verbose:
                        used_indices = where(accessor[index_queue[-1]] >= 0)[0]
                        o_options = [nucleotides[used_index] for used_index in used_indices]
                        b_options = [accessor[index_queue[-1]][used_index] for used_index in used_indices]
                        print("error occur in location " + str(location) + ", found " + detected_dna_string[location]
                              + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))
                        print("-" * 50)
                        print("Search error in current location.")
                        print("Check location = " + detected_dna_string[location] + " | " + str(location) + ", found "
                              + detected_dna_string[location] + " | "
                              + str(o_options).replace("\'", "") + ", " + str(b_options))

                    found_wrong, detected_flag_1 = True, True
                    for _, _, _, dna_string in path_matching(dna_string=detected_dna_string,
                                                             accessor=accessor, index_queue=index_queue,
                                                             occur_location=location,
                                                             has_insertion=has_insertion,
                                                             has_deletion=has_deletion,
                                                             verbose=verbose):
                        if dna_string not in new_dna_strings:
                            new_dna_strings.append(dna_string)

                    recall_queue = index_queue[::-1][:observed_length - 1]
                    for recall_index, original_vertex_index in enumerate(recall_queue):
                        error_index = location - recall_index - 1
                        if verbose:
                            print("-" * 50)
                            print("Search the error in the previous (" + str((recall_index + 1)) + ") location.")
                            print("Check location = " + detected_dna_string[error_index] + " | " + str(error_index))

                        for _, _, _, dna_string in path_matching(dna_string=detected_dna_string,
                                                                 accessor=accessor, index_queue=index_queue,
                                                                 occur_location=error_index,
                                                                 has_insertion=has_insertion,
                                                                 has_deletion=has_deletion,
                                                                 verbose=verbose):
                            if dna_string not in new_dna_strings:
                                new_dna_strings.append(dna_string)
                    break

            if (not found_wrong) and (detected_dna_string not in repaired_dna_strings):
                repaired_dna_strings.append(detected_dna_string)

        if 0 <= len(new_dna_strings) <= heap_size:
            dna_strings = new_dna_strings
        else:
            break

    repaired_results = []
    if len(repaired_dna_strings) > 0:
        for repaired_dna_string in repaired_dna_strings:
            if vt_check is not None:
                if vt_check == set_vt(dna_string=repaired_dna_string, vt_length=len(vt_check), nucleotides=nucleotides):
                    repaired_results.append(repaired_dna_string)
                else:
                    detected_flag_2 = True
            else:
                repaired_results.append(repaired_dna_string)

    if len(repaired_results) == 0:  # use original Varshamov-Tenengolts repair (only one error).
        if verbose:
            print()
            print("No result check, we use original Varshamov-Tenengolts repair.")
            print()

        # no mathematical method, using saturated attempt (substitute nucleotide in each location).
        for assume_location in range(len(original_dna_string)):
            for assume_nucleotide in list(filter(lambda n: n != original_dna_string[assume_location], nucleotides)):
                if verbose:
                    print("-" * 50)
                    print("Find the substitution error (" + str(assume_location) + ") location.")
                    print(original_dna_string[assume_location] + " needs to be substituted by "
                          + assume_nucleotide + ".")

                # substitution error checking through Varshamov-Tenengolts principle.
                repaired_dna_string = list(original_dna_string)
                repaired_dna_string[assume_location] = assume_nucleotide
                repaired_dna_string = "".join(repaired_dna_string)

                if vt_check == set_vt(dna_string=repaired_dna_string, vt_length=len(vt_check), nucleotides=nucleotides):
                    found_wrong, vertex_index = False, start_index
                    for location, nucleotide in enumerate(repaired_dna_string):
                        used_indices = where(accessor[vertex_index] >= 0)[0]
                        if nucleotide in [nucleotides[used_index] for used_index in used_indices]:
                            vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                        else:
                            found_wrong = True
                            break
                    if not found_wrong:
                        repaired_results.append(repaired_dna_string)

        if has_insertion:
            # insertion error checking through Varshamov-Tenengolts principle.
            for assume_location in range(len(original_dna_string)):
                if verbose:
                    nucleotide = original_dna_string[assume_location]
                    print("-" * 50)
                    print("Find the insertion error (" + str(assume_location) + ") location.")
                    print(nucleotide + " needs to be deleted.")

                repaired_dna_string = list(original_dna_string)
                del repaired_dna_string[assume_location]
                repaired_dna_string = "".join(repaired_dna_string)

                if vt_check == set_vt(dna_string=repaired_dna_string, vt_length=len(vt_check), nucleotides=nucleotides):
                    found_wrong, vertex_index = False, start_index
                    for location, nucleotide in enumerate(repaired_dna_string):
                        used_indices = where(accessor[vertex_index] >= 0)[0]
                        if nucleotide in [nucleotides[used_index] for used_index in used_indices]:
                            vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                        else:
                            found_wrong = True
                            break
                    if not found_wrong:
                        repaired_results.append(repaired_dna_string)

        if has_deletion:
            # deletion error checking through Varshamov-Tenengolts principle.
            for assume_location in range(len(original_dna_string) + 1):
                for assume_nucleotide in nucleotides:
                    if verbose:
                        print("-" * 50)
                        print("Find the deletion error (" + str(assume_location) + ") location.")
                        print(assume_nucleotide + " need to be inserted.")

                    repaired_dna_string = list(original_dna_string)
                    repaired_dna_string.insert(assume_location, assume_nucleotide)
                    repaired_dna_string = "".join(repaired_dna_string)

                    if vt_check == set_vt(dna_string=repaired_dna_string, vt_length=len(vt_check),
                                          nucleotides=nucleotides):
                        found_wrong, vertex_index = False, start_index
                        for location, nucleotide in enumerate(repaired_dna_string):
                            used_indices = where(accessor[vertex_index] >= 0)[0]
                            if nucleotide in [nucleotides[used_index] for used_index in used_indices]:
                                vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                            else:
                                found_wrong = True
                                break
                        if not found_wrong:
                            repaired_results.append(repaired_dna_string)

    return repaired_results, (detected_flag_1, detected_flag_2, len(repaired_dna_strings))


def find_vertices(observed_length, bio_filter, nucleotides=None, verbose=False):
    """
    Find valid vertices based on the given the biochemical constraints.

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
