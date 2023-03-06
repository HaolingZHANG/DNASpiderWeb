from itertools import combinations
from numpy import zeros, ones, zeros_like, array, union1d, min, median, max, random, log, log2, sum, abs, all, where

from dsw.operation import Monitor


def get_complete_accessor(observed_length, verbose=False):
    """
    Get a complete accessor with the required observed length.

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :param verbose: need to print log.
    :type verbose: bool

    :return: complete accessor.
    :rtype: numpy.ndarray

    Example
        >>> from dsw import get_complete_accessor
        >>> get_complete_accessor(observed_length=2)
        array([[ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15],
               [ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15],
               [ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15],
               [ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15]])

    .. note::
        The size of accessor is 4 ^ l * 4 and that of corresponding adjacency matrix is 4 ^ l * 4 ^ l.
    """
    accessor, monitor = -ones(shape=(int(4 ** observed_length), 4), dtype=int), Monitor()

    for vertex_index in range(int(4 ** observed_length)):
        latters = obtain_latters(current=vertex_index, observed_length=observed_length)
        for position, latter_vertex_index in enumerate(latters):
            accessor[vertex_index][position] = latter_vertex_index

        if verbose:
            monitor(vertex_index + 1, int(4 ** observed_length))

    return accessor


def accessor_to_adjacency_matrix(accessor, maximum_length=8, verbose=False):
    """
    Convert the accessor (compressed matrix) to its equivalent adjacency matrix.

    :param accessor: accessor (compressed matrix).
    :type: numpy.ndarray

    :param maximum_length: maximum vertex length (like 8 in general) of the adjacency matrix.
    :type maximum_length: int
    :param verbose: need to print log.
    :type verbose: bool

    :raise MemoryError: when you generate a large adjacency matrix that your memory cannot allocate.
    :raise ValueError: when you input a graph of DNA Spider-Web with wrong format.

    :return: adjacency matrix of the uncompressed graph.
    :rtype: numpy.ndarray

    Example
        >>> from dsw import accessor_to_adjacency_matrix, adjacency_matrix_to_accessor, get_complete_accessor
        >>> accessor = array([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15], \
                              [0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15], \
                              [0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15], \
                              [0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15]])
        >>> accessor_to_adjacency_matrix(accessor=accessor)
        array([[1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
               [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
               [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
               [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]])

    .. note::
        The size of accessor is 4 ^ l * 4 and that of corresponding adjacency matrix is 4 ^ l * 4 ^ l.
    """
    nucleotides = "ACGT"

    if len(accessor) >= 4 ** maximum_length:
        raise MemoryError("Unable to allocate adjacency matrix when length of DNA sequence (vertex) is more than 7.")
    if accessor.shape[1] != len(nucleotides) or min(accessor) < -1 or max(accessor) > len(accessor) - 1:
        raise ValueError("Wrong format in the accessor")

    matrix, monitor = zeros(shape=(len(accessor), len(accessor)), dtype=int), Monitor()
    for vertex_index, vertex in enumerate(accessor):
        matrix[vertex_index][vertex[vertex >= 0]] = 1

        if verbose:
            monitor(vertex_index + 1, len(accessor))

    return matrix


def adjacency_matrix_to_accessor(matrix, verbose=False):
    """
    Convert the adjacency matrix to the equivalent accessor (compressed matrix).

    :param matrix: adjacency matrix.
    :type matrix: numpy.ndarray

    :param verbose: need to print log.
    :type verbose: bool

    :raise ValueError: when you input an adjacency matrix with wrong format.

    :return: equivalent compressed matrix (accessor), compress rate = len(n_system) / len(matrix).
    :rtype: numpy.ndarray

    Example
        >>> from dsw import adjacency_matrix_to_accessor
        >>> adjacency_matrix = array([[1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1], \
                                      [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1], \
                                      [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1], \
                                      [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0], \
                                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]])
        >>> adjacency_matrix_to_accessor(matrix=adjacency_matrix)
        array([[ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15],
               [ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15],
               [ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15],
               [ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15]])

    .. note::
        The size of accessor is 4 ^ l * 4 and that of corresponding adjacency matrix is 4 ^ l * 4 ^ l.
    """
    nucleotides = "ACGT"

    accessor, monitor = -ones(shape=(len(matrix), len(nucleotides)), dtype=int), Monitor()
    observed_length = int(log(len(accessor)) / log(len(nucleotides)))

    for vertex_index, vertex in enumerate(matrix):
        next_indices = where(vertex == 1)[0].tolist()
        reference_latters = obtain_latters(current=vertex_index, observed_length=observed_length)
        if list(set(next_indices) | set(reference_latters)) != reference_latters:
            raise ValueError("Wrong format in the adjacency matrix, "
                             + "which cannot be converted to equivalent compressed accessor!")
        saved_information = [index if index in next_indices else -1 for index in reference_latters]
        accessor[vertex_index] = saved_information

        if verbose:
            monitor(vertex_index + 1, len(matrix))

    return accessor


def accessor_to_latter_map(accessor, verbose=False):
    """
    Convert the accessor to its equivalent latter map.

    :param accessor: accessor of graph.
    :type accessor: numpy.ndarray

    :param verbose: need to print log.
    :type verbose: bool

    :return: latter vertex map of graph.
    :rtype: dict

    Example
        >>> from numpy import array
        >>> from dsw import get_complete_accessor, accessor_to_latter_map
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> accessor_to_latter_map(accessor=accessor)
        {1: [4, 7], 2: [8, 11], 4: [1, 2], 7: [13, 14], 8: [1, 2], 11: [13, 14], 13: [4, 7], 14: [8, 11]}

    .. note::
        The size of accessor is 4 ^ l * 4.

        The size of corresponding latter map is further reduced,
        which only retains available information of follow-up vertices.
        However, latter map is not suitable for matrix calculation.
    """
    latter_map, monitor, total = {}, Monitor(), len(accessor)

    locations = where(sum(((accessor + 1).astype(bool)), axis=1).astype(int) > 0)[0]
    for index, location in enumerate(locations):
        vertex = accessor[location]
        latter_map[location] = vertex[vertex >= 0].tolist()

        if verbose:
            monitor(current_state=index + 1, total_state=len(locations))

    return latter_map


def latter_map_to_accessor(latter_map, observed_length, threshold=None, verbose=False):
    """
    Convert the latter map to the equivalent accessor.

    :param latter_map: latter vertex map of graph.
    :type latter_map: dict

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :param threshold: minimum out-degree threshold.
    :type threshold: int or None

    :param verbose: need to print log.
    :type verbose: bool

    :return: equivalent accessor.
    :rtype: numpy.ndarray

    Example
        >>> from dsw import latter_map_to_accessor
        >>> latter_map = {1: [4, 7], 2: [8, 11], 4: [1, 2], 7: [13, 14], \
                          8: [1, 2], 11: [13, 14], 13: [4, 7], 14: [8, 11]}
        >>> latter_map_to_accessor(latter_map=latter_map, observed_length=2)
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
        The size of accessor is 4 ^ l * 4.

        The size of corresponding latter map is further reduced,
        which only retains available information of follow-up vertices.
        However, latter map is not suitable for matrix calculation.
    """
    nucleotides = "ACGT"

    monitor = Monitor()

    if threshold is not None:
        latter_map = remove_useless(latter_map, threshold=threshold, verbose=verbose)

    accessor = -ones(shape=(len(nucleotides) ** observed_length, len(nucleotides)), dtype=int)
    if len(latter_map) > 0:
        if verbose:
            print("Convert the latter map to the accessor.")

        total = len(latter_map.items())
        for current, (former_vertex, latter_vertices) in enumerate(latter_map.items()):
            for latter_vertex in latter_vertices:
                accessor[former_vertex, latter_vertex % len(nucleotides)] = latter_vertex

            if verbose:
                monitor(current + 1, total)

    return accessor


def remove_useless(latter_map, threshold, verbose=False):
    """
    Remove useless vertices (the out-degree of witch less than threshold) in the latter map.

    :param latter_map: latter vertex map of graph.
    :type latter_map: dict

    :param threshold: minimum out-degree threshold.
    :type threshold: int

    :param verbose: need to print log.
    :type verbose: bool

    :return: useful latter map.
    :rtype: dict

    Example
        >>> from dsw import remove_useless
        >>> latter_map_1 = {0: [1, 2], 1: [], 2: [3, 4], 3: [0]}
        >>> remove_useless(latter_map=latter_map_1, threshold=1)
        {0: [2], 2: [3], 3: [0]}
        >>> latter_map_2 = {1: [4, 7], 2: [8, 11], 4: [1, 2], 7: [13, 14], \
                            8: [1, 2], 11: [13, 14], 13: [4, 7], 14: [8, 11]}
        >>> remove_useless(latter_map=latter_map_2, threshold=1)
        {1: [4, 7], 2: [8, 11], 4: [1, 2], 7: [13, 14], 8: [1, 2], 11: [13, 14], 13: [4, 7], 14: [8, 11]}
    """
    monitor = Monitor()

    if verbose:
        print("Remove useless vertex, the out-degree of witch less than " + str(threshold) + ".")

    round_number = 1
    while True:
        if verbose:
            print("Check available vertices.")

        remove_vertices, saved_vertices, total = [], [], len(latter_map)
        for current, (former_vertex, latter_vertices) in enumerate(latter_map.items()):
            if len(latter_vertices) < threshold:
                remove_vertices.append(former_vertex)
            else:
                saved_vertices.append(former_vertex)

            if verbose:
                monitor(current + 1, total, extra={"round": round_number})

        if verbose:
            print("Remove vertices " + str(remove_vertices) + " and load a novel latter map.")

        new_latter_map, remove_flag = {}, False
        for current, (former_vertex, latter_vertices) in enumerate(latter_map.items()):
            if former_vertex not in remove_vertices:
                available_latter_vertices = []
                for latter_vertex in latter_vertices:
                    if (latter_vertex not in remove_vertices) and (latter_vertex in saved_vertices):
                        available_latter_vertices.append(latter_vertex)
                    else:
                        remove_flag = True
                new_latter_map[former_vertex] = available_latter_vertices

            if verbose:
                monitor(current + 1, total, extra={"round": round_number})

        latter_map = new_latter_map

        round_number += 1

        if not remove_flag:
            break

    return latter_map


def obtain_formers(current, observed_length):
    """
    Obtain former vertex given_amino_acids based on the current vertex index.

    :param current: current vertex index.
    :type current: int

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :return: former vertex given_amino_acids.
    :rtype: list

    Example
        >>> from dsw import obtain_formers
        >>> current = 0
        >>> obtain_formers(current=current, observed_length=2)
        [0, 4, 8, 12]
        >>> obtain_formers(current=current, observed_length=3)
        [0, 16, 32, 48]
        >>> obtain_formers(current=current, observed_length=4)
        [0, 64, 128, 192]
        >>> obtain_formers(current=current, observed_length=5)
        [0, 256, 512, 768]
        >>> obtain_formers(current=current, observed_length=6)
        [0, 1024, 2048, 3072]
        >>> obtain_formers(current=current, observed_length=7)
        [0, 4096, 8192, 12288]
        >>> obtain_formers(current=current, observed_length=8)
        [0, 16384, 32768, 49152]
        >>> obtain_formers(current=current, observed_length=9)
        [0, 65536, 131072, 196608]
    """
    nucleotides = "ACGT"

    formers = []
    for former_value in range(len(nucleotides)):
        former = current // len(nucleotides) + former_value * int(len(nucleotides) ** (observed_length - 1))
        formers.append(former)

    return formers


def obtain_latters(current, observed_length):
    """
    Obtain latter vertex given_amino_acids based on the current vertex index.

    :param current: current vertex index.
    :type current: int

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :return: latter vertex given_amino_acids.
    :rtype: list

    Example
        >>> from dsw import obtain_latters
        >>> current = 0
        >>> obtain_latters(current=current, observed_length=2)
        [0, 1, 2, 3]
        >>> obtain_latters(current=current, observed_length=3)
        [0, 1, 2, 3]
        >>> obtain_latters(current=current, observed_length=4)
        [0, 1, 2, 3]
        >>> obtain_latters(current=current, observed_length=5)
        [0, 1, 2, 3]
        >>> obtain_latters(current=current, observed_length=6)
        [0, 1, 2, 3]
        >>> obtain_latters(current=current, observed_length=7)
        [0, 1, 2, 3]
        >>> obtain_latters(current=current, observed_length=8)
        [0, 1, 2, 3]
        >>> obtain_latters(current=current, observed_length=9)
        [0, 1, 2, 3]
    """
    nucleotides = "ACGT"

    latters = []
    for latter_value in range(len(nucleotides)):
        latter = int((current * len(nucleotides) + latter_value) % (len(nucleotides) ** observed_length))
        latters.append(latter)

    return latters


def obtain_vertices(accessor):
    """
    Obtain available vertices from the established graph.

    :param accessor: accessor of graph.
    :type accessor: numpy.ndarray

    :return: vertex list.
    :rtype: numpy.ndarray

    Example
        >>> from numpy import array
        >>> from dsw import obtain_vertices
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> obtain_vertices(accessor=accessor)
        array([ 1,  2,  4,  7,  8, 11, 13, 14])
    """
    return where(sum(((accessor + 1).astype(bool)), axis=1).astype(bool) == 1)[0].astype(int)


def obtain_leaf_vertices(vertex_index, depth, accessor=None, latter_map=None):
    """
    Obtain leaf vertices in required depth of the tree with the rooted vertex.

    :param vertex_index: vertex index in the graph.
    :type vertex_index: int

    :param depth: stride of the breadth-first search.
    :type depth: int

    :param accessor: accessor of graph.
    :type accessor: numpy.ndarray

    :param latter_map: latter vertex map of graph.
    :type latter_map: dict

    :return: given_amino_acids of required leaf vertex.
    :rtype: numpy.ndarray

    Example
        >>> from numpy import array
        >>> from dsw import obtain_leaf_vertices
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> obtain_leaf_vertices(vertex_index=1, depth=1, accessor=accessor)
        array([4, 7])
        >>> obtain_leaf_vertices(vertex_index=1, depth=2, accessor=accessor)
        array([ 1,  2, 13, 14])
        >>> latter_map = {1: [4, 7], 2: [8, 11], 4: [1, 2], 7: [13, 14], \
                          8: [1, 2], 11: [13, 14], 13: [4, 7], 14: [8, 11]}
        >>> obtain_leaf_vertices(vertex_index=1, depth=1, latter_map=latter_map)
        array([4, 7])
        >>> obtain_leaf_vertices(vertex_index=1, depth=2, latter_map=latter_map)
        array([ 1,  2, 13, 14])

    .. note::
        Either the parameter "accessor" or the parameter "latter_map" must occur, but not both.
    """
    if (accessor is not None) and (latter_map is not None):
        raise ValueError("Too many variables (accessor and latter map) are assigned, "
                         + "we do not know which variable needs to be used as a priority!")

    branch = [vertex_index]

    if accessor is not None:
        for step in range(depth):  # do breadth-first search
            level = []
            for former_index in branch:
                vertex = accessor[former_index]
                available_latters = vertex[vertex >= 0].tolist()
                level += available_latters
            branch = level

    elif latter_map is not None:
        for step in range(depth):  # do breadth-first search
            level = []
            for former_index in branch:
                if former_index in latter_map:
                    for latter_index in latter_map[former_index]:
                        level.append(latter_index)
            branch = level

    else:
        raise ValueError("We need to select a data type (accessor and latter map) input for the graph!")

    return array(branch)


# noinspection PyUnresolvedReferences
def approximate_capacity(accessor, tolerance_level=-10, repeats=1, maximum_iteration=500, process=False, verbose=False):
    """
    Approximate the capacity of the specific graph through Perron–Frobenius theorem.

    :param accessor: accessor of graph.
    :type accessor: numpy.ndarray

    :param tolerance_level: error tolerance of power iteration.
    :type tolerance_level: int

    :param repeats: random repeats for approximating the capacity.
    :type repeats: int

    :param maximum_iteration: maximum iteration in the power method.
    :type maximum_iteration: int

    :param process: need eigenvalue in the process.
    :type process: bool

    :param verbose: need to print log.
    :type verbose: bool

    :return: capacity of this graph (accessor) and process values if required.
    :rtype: float or (float, list)

    Example
        >>> from numpy import array, random
        >>> from dsw import approximate_capacity
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> approximate_capacity(accessor=accessor, tolerance_level=-10, repeats=2, process=False)
        1.0
        >>> random.seed(0)
        >>> capacity, processes = approximate_capacity(accessor=accessor, tolerance_level=-10, repeats=2, process=True)
        >>> capacity
        1.0
        >>> ["%.5f" % _ for _ in processes[0]]
        ['0.57779', '0.91175', '1.00000', '1.00000']
        >>> ["%.5f" % _ for _ in processes[1]]
        ['0.81488', '0.68189', '1.00000', '1.00000']

    .. note::
        Reference [1] Oskar Perron (1907) Mathematische Annalen

        Reference [2] Georg Frobenius (1912) Sitzungsberichte der Königlich Preussischen Akademie der Wissenschaften

        Reference [3] Brian H. Marcus et al. (2001) Lecture notes

        Reference [4] Nabil Kahale (1995) Journal of the ACM

        Reference [5] William Ford (2014) Academic Press
    """
    if all(accessor == -1):
        if process:
            return (0.0, [0.0]) if repeats == 1 else (0.0, [[0.0] for _ in range(repeats)])
        else:
            return 0.0

    ignore_positions = where(sum(accessor, axis=1) == -len(accessor[0]))[0]

    results, record = [], []
    for repeat in range(repeats):
        if verbose and repeats > 1:
            print("Approximate capacity in " + str(repeat + 1) + " (" + str(repeats) + ") times.")

        record.append([])
        if repeats > 1:
            last_eigenvector = abs(random.random(size=(len(accessor),)))  # Random initialization for a faster fitness.
        else:
            last_eigenvector = ones(shape=(len(accessor),), dtype=float)

        last_eigenvector[ignore_positions] = 0.0  # refers to the vertex without follow-up vertices.

        monitor, queue, last_eigenvalue, current = Monitor(), [], None, 0
        while True:
            eigenvector = zeros_like(last_eigenvector)
            for positions in accessor.T:
                available = where(positions >= 0)
                eigenvector[available] += last_eigenvector[positions[available]]
            eigenvalue = max(eigenvector)
            if eigenvalue > 0:
                eigenvector = eigenvector / eigenvalue
            else:
                eigenvector = eigenvector * 0.0
            record[-1].append(log2(eigenvalue) if eigenvalue > 10 ** tolerance_level else 0.0)

            if last_eigenvalue is not None:
                if last_eigenvalue > 0.0:
                    relative_error = abs(eigenvalue - last_eigenvalue) / last_eigenvalue
                else:
                    relative_error = 0.0
                queue.append(eigenvalue)

                if verbose and current + 1 < maximum_iteration:
                    monitor(current + 1, maximum_iteration,
                            extra={"largest eigenvalue": "%.5f" % eigenvalue, "error": "%.5f" % relative_error})

                is_finished = False
                if relative_error < 10 ** tolerance_level:
                    results.append(log2(eigenvalue) if eigenvalue > 10 ** tolerance_level else 0.0)
                    is_finished = True

                if len(queue) > maximum_iteration:
                    eigenvalue = median(queue)
                    results.append(log2(eigenvalue) if eigenvalue > 10 ** tolerance_level else 0.0)
                    is_finished = True

                if is_finished:
                    if verbose:
                        monitor(maximum_iteration, maximum_iteration, extra={"capacity": "%.5f" % results[-1]})
                    break

            last_eigenvalue, last_eigenvector, current = eigenvalue, eigenvector, current + 1

    if process:
        return (median(results), record[0]) if repeats == 1 else (median(results), record)
    else:
        return median(results)


def path_matching(dna_sequence, accessor, previous_index, occur_location, has_indel=False, nucleotides=None):
    """
    Perform saturation repair at the selected position and obtain the DNA sequences matching the path of accessor.

    :param dna_sequence: DNA sequence waiting for saturation substitution in the specific location.
    :type dna_sequence: str

    :param accessor: accessor.
    :type accessor: numpy.ndarray

    :param previous_index: previous vertex index before the occurred error location.
    :type previous_index: int

    :param occur_location: the location that may occurring substitution (or specific location).
    :type occur_location: int

    :param has_indel: consider insertion and/or deletion error.
    :type has_indel: bool

    :param nucleotides: usage of nucleotides.
    :type nucleotides: str or None

    :return: repaired DNA sequences (may contain multiple repair show) and visited count.
    :rtype: list, int

    Example
        >>> from numpy import array
        >>> from dsw import path_matching
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> dna_sequence = "TCTCTATCTCT"  # "TCTCTCTCTCT" is original DNA sequence
        >>> path_matching(dna_sequence=dna_sequence, accessor=accessor, previous_index=7, occur_location=5, \
                          has_indel=True)
        ([(('S', 5, 'C'), 'TCTCTCTCTCT'), (('S', 5, 'G'), 'TCTCTGTCTCT')], 12)
    """
    if nucleotides is None:
        nucleotides = "ACGT"

    repair_info, visited_count = [], 0
    original, used_indices = dna_sequence[occur_location], where(accessor[previous_index] >= 0)[0]

    for r_nucleotide in list(filter(lambda n: n != original, [nucleotides[index] for index in used_indices])):
        vertex_index, reliable = accessor[previous_index][nucleotides.index(r_nucleotide)], True
        for index, nucleotide in enumerate(dna_sequence[occur_location + 1:]):
            used_nucleotides = [nucleotides[used_index] for used_index in where(accessor[vertex_index] >= 0)[0]]
            if nucleotide in used_nucleotides:
                vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                visited_count += 1
            else:
                reliable = False
                break

        if reliable:  # "S" refers to repair by substation.
            obtained_dna_sequence = list(dna_sequence)
            obtained_dna_sequence[occur_location] = r_nucleotide
            repair_info.append((("S", occur_location, r_nucleotide), "".join(obtained_dna_sequence)))

    if has_indel:
        for a_nucleotide in [nucleotides[used_index] for used_index in used_indices]:
            vertex_index, reliable = accessor[previous_index][nucleotides.index(a_nucleotide)], True
            for nucleotide in dna_sequence[occur_location:]:
                used_nucleotides = [nucleotides[used_index] for used_index in where(accessor[vertex_index] >= 0)[0]]
                if nucleotide in used_nucleotides:
                    vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                    visited_count += 1
                else:
                    reliable = False
                    break

            if reliable:  # "I" refers to repair by insertion.
                obtained_dna_sequence = list(dna_sequence)
                obtained_dna_sequence.insert(occur_location, a_nucleotide)
                repair_info.append((("I", occur_location, a_nucleotide), "".join(obtained_dna_sequence)))

        d_nucleotide, vertex_index, reliable = original, previous_index, True
        for index, nucleotide in enumerate(dna_sequence[occur_location + 1:]):
            used_nucleotides = [nucleotides[used_index] for used_index in where(accessor[vertex_index] >= 0)[0]]
            if nucleotide in used_nucleotides:
                vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                visited_count += 1
            else:
                reliable = False
                break

        if reliable:   # "D" refers to repair by deletion.
            obtained_dna_sequence = list(dna_sequence)
            del obtained_dna_sequence[occur_location]
            repair_info.append((("D", occur_location, d_nucleotide), "".join(obtained_dna_sequence)))

    return repair_info, visited_count


def calculate_intersection_score(latter_map, observed_length=10, has_insertion=True, has_deletion=True, verbose=False):
    """
    Calculate the intersection score based on the breach-first search.

    :param latter_map: latter vertex map of graph.
    :type latter_map: dict

    :param observed_length: length of the DNA sequence in a vertex.
    :type observed_length: int

    :param has_insertion: consider to repair insertion errors.
    :type has_insertion: bool

    :param has_deletion: consider to repair deletion errors.
    :type has_deletion: bool

    :param verbose: need to print log.
    :type verbose: bool

    :return: intersection scores for each arc in the coding graph (consistent with the shape of accessor).
    :rtype: numpy.ndarray

    Example
        >>> from dsw import calculate_intersection_score
        >>> # latter_map with GC-balanced
        >>> latter_map = {1: [4, 7], 2: [8, 11], 4: [1, 2], 7: [13, 14], \
                          8: [1, 2], 11: [13, 14], 13: [4, 7], 14: [8, 11]}
        >>> calculate_intersection_score(latter_map, observed_length=10, \
                                         has_insertion=True, has_deletion=True, verbose=False)
        array([[ 0,  0,  0,  0],
               [28,  0,  0, 28],
               [28,  0,  0, 28],
               ...,
               [ 0,  0,  0,  0],
               [ 0,  0,  0,  0],
               [ 0,  0,  0,  0]])

    .. note::
        It is a gift for the follow-up investigation.
        That is, removing arc to improve the capability of the probabilistic error correction.
    """
    nucleotides = "ACGT"

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

        if has_insertion:
            for index, former_index in enumerate(latter_map[current_index]):  # insertion
                if former_index in latter_map:
                    for latter_index in latter_map[former_index]:
                        insert_branch = obtain_leaf_vertices(latter_index, depth, latter_map=latter_map)
                        score = len(union1d(mutate_branches[index], insert_branch))
                        scores[current_index, latter_map[current_index][index] % len(nucleotides)] += score

        if has_deletion:
            delete_branch = [obtain_leaf_vertices(current_index, depth, latter_map=latter_map)]  # deletion
            for index in range(len(mutate_branches)):
                score = len(union1d(mutate_branches[index], delete_branch))
                scores[current_index, latter_map[current_index][index] % len(nucleotides)] += score
            del delete_branch

        if verbose:
            monitor(current + 1, len(currents))

    return scores
