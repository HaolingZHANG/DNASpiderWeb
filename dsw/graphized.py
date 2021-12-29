__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from numpy import zeros, ones, zeros_like, array, min, median, max, random, log, log2, log10, abs, all, where

from dsw.operation import Monitor


def get_complete_accessor(observed_length, verbose=False):
    """
    Get a complete accessor with the required observed length.

    :param observed_length: length of the DNA string in a vertex.
    :type observed_length: int

    :param verbose: need to print log.
    :type verbose: bool

    :return: complete accessor.
    :rtype: numpy.ndarray

    ..example::
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
        """
    accessor, monitor = -ones(shape=(int(4 ** observed_length), 4), dtype=int), Monitor()

    for vertex_index in range(int(4 ** observed_length)):
        latters = obtain_latters(current=vertex_index, observed_length=observed_length)
        for position, latter_vertex_index in enumerate(latters):
            accessor[vertex_index][position] = latter_vertex_index

        if verbose:
            monitor.output(vertex_index + 1, int(4 ** observed_length))

    return accessor


def accessor_to_adjacency_matrix(accessor, maximum_length=8, nucleotides=None, verbose=False):
    """
    Convert the accessor (compressed matrix) to its equivalent adjacency matrix.

    :param accessor: accessor (compressed matrix).
    :type: numpy.ndarray

    :param maximum_length: maximum vertex length (like 8 in general) of the adjacency matrix (size is (4^l)^2).
    :type maximum_length: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list or None

    :param verbose: need to print log.
    :type verbose: bool

    :raise MemoryError: when you generate a large adjacency matrix that your memory cannot allocate.
    :raise ValueError: when you input a graph of DNA Spider-Web with wrong format.

    :return: adjacency matrix of the uncompressed graph.
    :rtype: numpy.ndarray

    ..example::
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
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    if len(accessor) >= 4 ** maximum_length:
        raise MemoryError("Unable to allocate adjacency matrix when length of DNA string (vertex) is more than 7.")
    if accessor.shape[1] != len(nucleotides) or min(accessor) < -1 or max(accessor) > len(accessor) - 1:
        raise ValueError("Wrong format in the accessor")

    matrix, monitor = zeros(shape=(len(accessor), len(accessor)), dtype=int), Monitor()
    for vertex_index, vertex in enumerate(accessor):
        matrix[vertex_index][vertex[vertex >= 0]] = 1

        if verbose:
            monitor.output(vertex_index + 1, len(accessor))

    return matrix


def adjacency_matrix_to_accessor(matrix, nucleotides=None, verbose=False):
    """
    Convert the adjacency matrix to the equivalent accessor (compressed matrix).

    :param matrix: adjacency matrix.
    :type matrix: numpy.ndarray

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list or None

    :param verbose: need to print log.
    :type verbose: bool

    :raise ValueError: when you input an adjacency matrix with wrong format.

    :return: equivalent compressed matrix (accessor), compress rate = len(n_system) / len(matrix).
    :rtype: numpy.ndarray

    ..example::
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
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

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
            monitor.output(vertex_index + 1, len(matrix))

    return accessor


def accessor_to_latter_map(accessor, verbose=False):
    """
    Convert the accessor to the equivalent latter map.

    :param accessor: accessor of graph.
    :type accessor: numpy.ndarray

    :param verbose: need to print log.
    :type verbose: bool

    :return: latter vertices map of graph.
    :rtype: dict
    
    ..example::
        >>> from numpy import array
        >>> from dsw import get_complete_accessor, accessor_to_latter_map
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> accessor_to_latter_map(accessor=accessor)
        {1: [4, 7], 2: [8, 11], 4: [1, 2], 7: [13, 14], 8: [1, 2], 11: [13, 14], 13: [4, 7], 14: [8, 11]}
    """
    latter_map, monitor, total = {}, Monitor(), len(accessor)
    for current, vertex in enumerate(accessor):
        if not all(vertex == -1):
            latters = vertex[vertex >= 0].tolist()
            if len(latters) > 0:
                latter_map[current] = latters

        if verbose:
            monitor.output(current_state=current + 1, total_state=total)

    return latter_map


def latter_map_to_accessor(latter_map, observed_length, threshold=None, nucleotides=None, verbose=False):
    """
    Convert the latter map to the equivalent accessor.

    :param latter_map: latter vertices map of graph.
    :type latter_map: dict

    :param observed_length: length of the DNA string in a vertex.
    :type observed_length: int

    :param threshold: minimum out-degree threshold.
    :type threshold: int or None

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list or None

    :param verbose: need to print log.
    :type verbose: bool

    :return: equivalent accessor.
    :rtype: numpy.ndarray

    ..example::
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
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

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
                monitor.output(current + 1, total)

    return accessor


def remove_useless(latter_map, threshold, verbose=False):
    """
    Remove useless vertices (the out-degree of witch less than threshold) in the latter map.

    :param latter_map: latter vertices map of graph.
    :type latter_map: dict

    :param threshold: minimum out-degree threshold.
    :type threshold: int

    :param verbose: need to print log.
    :type verbose: bool

    :return: useful latter map.
    :rtype: dict

    ..example::
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
                monitor.output(current + 1, total, extra={"round": round_number})

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
                monitor.output(current + 1, total, extra={"round": round_number})

        latter_map = new_latter_map

        round_number += 1

        if not remove_flag:
            break

    return latter_map


def obtain_formers(current, observed_length, nucleotides=None):
    """
    Obtain former vertex indices based on the current vertex index.

    :param current: current vertex index.
    :type current: int

    :param observed_length: length of the DNA string in a vertex.
    :type observed_length: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :return: former vertex indices.
    :rtype: list

    ..example::
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
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    formers = []
    for former_value in range(len(nucleotides)):
        former = current // len(nucleotides) + former_value * int(len(nucleotides) ** (observed_length - 1))
        formers.append(former)

    return formers


def obtain_latters(current, observed_length, nucleotides=None):
    """
    Obtain latter vertex indices based on the current vertex index.

    :param current: current vertex index.
    :type current: int

    :param observed_length: length of the DNA string in a vertex.
    :type observed_length: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :return: latter vertex indices.
    :rtype: list

    ..example::
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
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    latters = []
    for latter_value in range(len(nucleotides)):
        latter = int((current * len(nucleotides) + latter_value) % (len(nucleotides) ** observed_length))
        latters.append(latter)

    return latters


def obtain_leaf_vertices(vertex_index, depth, accessor=None, latter_map=None):
    """
    Obtain leaf vertices in required depth of the tree with the rooted vertex.

    :param vertex_index: vertex index in the graph.
    :type vertex_index: int

    :param depth: stride of the breadth-first search.
    :type depth: int

    :param accessor: accessor of graph.
    :type accessor: numpy.ndarray

    :param latter_map: latter vertices map of graph.
    :type latter_map: dict

    :return: indices of required leaf vertex.
    :rtype: numpy.ndarray

    ..note::
        Either the parameter "accessor" or the parameter "latter_map" must occur, but not both.

    ..example::
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


def approximate_capacity(accessor, tolerance_level=-10, repeats=1, maximum_iteration=500,
                         need_process=False, verbose=False):
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

    :param need_process: need eigenvalue in the process.
    :type need_process: bool

    :param verbose: need to print log.
    :type verbose: bool

    :return: capacity of this graph (accessor) and process values if required.
    :rtype: float or (float, list)

    ..note::
        Reference [1] Oskar Perron (1907) Mathematische Annalen.
        Reference [2] Georg Frobenius (1912) Sitzungsberichte der Königlich Preussischen Akademie der Wissenschaften.
        Reference [3] Brian H. Marcus et al. (2001) Lecture notes.

    ..example::
        >>> from numpy import array, random
        >>> from dsw import approximate_capacity
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], \
                              [ 4, -1, -1,  7], \
                              [ 8, -1, -1, 11], \
                              [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], \
                              [-1, -1, -1, -1], \
                              [-1, -1, -1, -1], \
                              [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], \
                              [-1, -1, -1, -1], \
                              [-1, -1, -1, -1], \
                              [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], \
                              [ 4, -1, -1,  7], \
                              [ 8, -1, -1, 11], \
                              [-1, -1, -1, -1]])
        >>> approximate_capacity(accessor=accessor, tolerance_level=-10, repeats=10)
        1.0
        >>> random.seed(0)
        >>> capacity, processes = approximate_capacity(accessor=accessor, tolerance_level=-10, repeats=2, \
                                                       need_process=True)
        >>> capacity
        1.0
        >>> processes
        [[0.5777866858675372, 0.9117487898585336, 1.0, 1.0], [0.814876469343339, 0.6818879690650359, 1.0, 1.0]]
    """
    if all(accessor == -1):
        if need_process:
            return (0.0, [0.0]) if repeats == 1 else (0.0, [[0.0] for _ in range(repeats)])
        else:
            return 0.0

    results, process = [], []
    for repeat in range(repeats):
        if verbose and repeats > 1:
            print("Approximate capacity in " + str(repeat + 1) + " (" + str(repeats) + ") times.")

        process.append([])
        if repeats > 1:
            last_eigenvector, last_eigenvalue = abs(random.random(size=(len(accessor),))), None
        else:
            last_eigenvector, last_eigenvalue = ones(shape=(len(accessor),), dtype=float), None

        for vertex_index, vertex in enumerate(accessor):  # remove ignored vertex.
            if all(vertex == -1):
                last_eigenvector[vertex_index] = 0.0

        monitor, queue = Monitor(), []
        while True:
            eigenvector = zeros_like(last_eigenvector)
            for positions in accessor.T:
                available = where(positions >= 0)
                eigenvector[available] += last_eigenvector[positions[available]]
            eigenvalue = max(eigenvector)
            eigenvector = eigenvector / eigenvalue
            process[-1].append(log2(eigenvalue))

            if last_eigenvalue is not None:
                relative_error = abs(eigenvalue - last_eigenvalue) / last_eigenvalue
                queue.append(eigenvalue)

                if verbose:
                    if 0.1 <= relative_error:
                        monitor.output(1, -tolerance_level,
                                       extra={"largest eigenvalue": "%.5f" % eigenvalue,
                                              "relative error": "%.5f" % relative_error})
                    elif 10 ** tolerance_level < relative_error < 0.1:
                        monitor.output(int(-log10(relative_error)), -tolerance_level,
                                       extra={"largest eigenvalue": "%.5f" % eigenvalue,
                                              "relative error": "%.5f" % relative_error})

                if relative_error < 10 ** tolerance_level:
                    if eigenvalue < 1.0:
                        eigenvalue = 1.0
                    results.append(log2(eigenvalue))
                    if verbose:
                        monitor.output(-tolerance_level, -tolerance_level,
                                       extra={"capacity": "%.5f" % results[-1]})
                    break

                if len(queue) > maximum_iteration:
                    eigenvalue = median(queue)
                    if eigenvalue < 1.0:
                        eigenvalue = 1.0
                    results.append(log2(eigenvalue))
                    if verbose:
                        monitor.output(-tolerance_level, -tolerance_level,
                                       extra={"capacity": "%.5f" % results[-1]})
                    break

            last_eigenvalue, last_eigenvector = eigenvalue, eigenvector

    if need_process:
        return (median(results), process[0]) if repeats == 1 else (median(results), process)
    else:
        return median(results)


def repair_by_match(dna_string, accessor, index_queue, occur_location, nucleotides=None,
                    has_insertion=False, has_deletion=False, verbose=False):
    """
    Perform saturation repair at the selected position and obtain the DNA strings meeting conditions of accessor.

    :param dna_string: DNA string waiting for saturation substitution in the specific location.
    :type dna_string: str

    :param accessor: accessor.
    :type accessor: numpy.ndarray

    :param index_queue: queue of vertex indices (transcoding path in the graph).
    :type index_queue: list

    :param occur_location: the location that may occurring substitution (or specific location).
    :type occur_location: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :param has_insertion: consider insertion error.
    :type has_insertion: bool

    :param has_deletion: consider deletion error.
    :type has_deletion: bool

    :param verbose: need to print log.
    :type verbose: bool

    :return: repaired DNA strings (may contain multiple repair results).
    :rtype: list

    ..example::
        >>> from numpy import array
        >>> from dsw import repair_by_match
        >>> # accessor with GC-balanced
        >>> accessor = array([[-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1,  1,  2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1], \
                              [-1, -1, -1, -1], [ 4, -1, -1,  7], [ 8, -1, -1, 11], [-1, -1, -1, -1]])
        >>> dna_string = "TCTCTATCTCTC"  # "TCTCTCTCTCTC" is original DNA string
        >>> repair_by_match(dna_string=dna_string, accessor=accessor, index_queue=[1, 7, 13, 7, 13, 7], \
                            occur_location=4, has_insertion=True, has_deletion=True)
        [('D', 4, 'T', 'TCTCATCTCTC')]
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    repair_info, observed_length = [], int(log(len(accessor)) / log(len(nucleotides)))
    original = dna_string[occur_location]

    if occur_location + observed_length < len(dna_string):
        check_segment = dna_string[occur_location: occur_location + observed_length]
    else:
        check_segment = dna_string[occur_location:]

    if verbose:
        print("Original = [" + dna_string + "]")
        if occur_location + observed_length < len(dna_string):
            remain = len(dna_string) - (occur_location + observed_length)
            print("Checking = [" + "-" * occur_location + check_segment + "-" * remain + "]")
        else:
            print("Checking = [" + "-" * occur_location + check_segment + "]")

    # find substitution error.
    used_indices = where(accessor[index_queue[occur_location]] >= 0)[0]
    o_options = [nucleotides[used_index] for used_index in used_indices]
    b_options = accessor[index_queue[occur_location]][~(accessor[index_queue[occur_location]] == -1)]
    for replace_nucleotide in list(filter(lambda n: n != original, o_options)):
        segment = "".join([replace_nucleotide] + list(check_segment)[1:])

        if verbose:
            print("Replace " + original + " to " + replace_nucleotide)
            if occur_location + observed_length < len(dna_string):
                remain = len(dna_string) - (occur_location + observed_length)
                print("         = [" + "-" * occur_location + segment + "-" * remain + "]"
                      + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))
            else:
                print("         = [" + "-" * occur_location + segment + "]"
                      + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))

        vertex_index, reliable = index_queue[occur_location], True
        for index, nucleotide in enumerate(segment):
            inner_used_indices = where(accessor[vertex_index] >= 0)[0]
            used_nucleotides = [nucleotides[used_index] for used_index in inner_used_indices]
            if nucleotide in used_nucleotides:
                if verbose:
                    show_used_indices = where(accessor[vertex_index] >= 0)[0]
                    o_options = [nucleotides[used_index] for used_index in show_used_indices]
                    b_options = accessor[vertex_index][~(accessor[vertex_index] == -1)]
                    if occur_location + observed_length < len(dna_string):
                        remain = len(dna_string) - (occur_location + observed_length)
                        print("           [" + "-" * (occur_location + index) + segment[index]
                              + "-" * (remain + len(segment) - index - 1) + "]"
                              + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))
                    else:
                        print("           [" + "-" * (occur_location + index) + segment[index] + "]"
                              + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))

                vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
            else:
                reliable = False
                break

        if reliable:
            obtained_dna_string = list(dna_string)
            obtained_dna_string[occur_location] = replace_nucleotide
            obtained_dna_string = "".join(obtained_dna_string)
            # "S" refers to repair by substation.
            repair_info.append(("S", occur_location, replace_nucleotide, obtained_dna_string))

    if has_insertion:  # find insertion error if required.
        used_indices = where(accessor[index_queue[occur_location]] >= 0)[0]
        o_options = [nucleotides[used_index] for used_index in used_indices]
        b_options = accessor[index_queue[occur_location]][~(accessor[index_queue[occur_location]] == -1)]
        for add_nucleotide in o_options:
            segment = "".join([add_nucleotide] + list(check_segment))
            if verbose:
                print("Insert " + add_nucleotide + " in location " + str(occur_location))
                if occur_location + observed_length < len(dna_string):
                    remain = len(dna_string) - (occur_location + observed_length)
                    print("         = [" + "-" * occur_location + segment + "-" * remain + "]"
                          + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))
                else:
                    print("         = [" + "-" * occur_location + segment + "]"
                          + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))

            vertex_index, reliable = index_queue[occur_location], True
            for index, nucleotide in enumerate(segment):
                inner_used_indices = where(accessor[vertex_index] >= 0)[0]
                used_nucleotides = [nucleotides[used_index] for used_index in inner_used_indices]
                if nucleotide in used_nucleotides:
                    if verbose:
                        show_used_indices = where(accessor[vertex_index] >= 0)[0]
                        o_options = [nucleotides[used_index] for used_index in show_used_indices]
                        b_options = accessor[vertex_index][~(accessor[vertex_index] == -1)]
                        if occur_location + observed_length < len(dna_string):
                            remain = len(dna_string) - (occur_location + observed_length)
                            print("           [" + "-" * (occur_location + index) + segment[index]
                                  + "-" * (remain + len(segment) - index - 1) + "]"
                                  + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))
                        else:
                            print("           [" + "-" * (occur_location + index) + segment[index] + "]"
                                  + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))

                    vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
                else:
                    reliable = False
                    break

            if reliable:
                obtained_dna_string = list(dna_string)
                obtained_dna_string.insert(occur_location, add_nucleotide)
                obtained_dna_string = "".join(obtained_dna_string)
                # "I" refers to repair by insertion.
                repair_info.append(("I", occur_location, add_nucleotide, obtained_dna_string))

    if has_deletion:  # find deletion error if required.
        used_indices = where(accessor[index_queue[occur_location]] >= 0)[0]
        o_options = [nucleotides[used_index] for used_index in used_indices]
        b_options = accessor[index_queue[occur_location]][~(accessor[index_queue[occur_location]] == -1)]

        delete_nucleotide = check_segment[0]
        segment = "".join(list(check_segment)[1:])
        if verbose:
            print("Delete " + check_segment[0] + " in location " + str(occur_location))
            if occur_location + observed_length < len(dna_string):
                remain = len(dna_string) - (occur_location + observed_length)
                print("         = [" + "-" * occur_location + segment + "-" * remain + "]"
                      + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))
            else:
                print("         = [" + "-" * occur_location + segment + "]"
                      + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))

        vertex_index, reliable = index_queue[occur_location], True
        for index, nucleotide in enumerate(segment):
            inner_used_indices = where(accessor[vertex_index] >= 0)[0]
            used_nucleotides = [nucleotides[used_index] for used_index in inner_used_indices]
            if nucleotide in used_nucleotides:
                if verbose:
                    show_used_indices = where(accessor[vertex_index] >= 0)[0]
                    o_options = [nucleotides[used_index] for used_index in show_used_indices]
                    b_options = accessor[vertex_index][~(accessor[vertex_index] == -1)]
                    if occur_location + observed_length < len(dna_string):
                        remain = len(dna_string) - (occur_location + observed_length)
                        print("           [" + "-" * (occur_location + index) + segment[index]
                              + "-" * (remain + len(segment) - index - 1) + "]"
                              + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))
                    else:
                        print("           [" + "-" * (occur_location + index) + segment[index] + "]"
                              + " | " + str(o_options).replace("\'", "") + ", " + str(b_options))

                vertex_index = accessor[vertex_index][nucleotides.index(nucleotide)]
            else:
                reliable = False
                break

        if reliable:
            obtained_dna_string = list(dna_string)
            del obtained_dna_string[occur_location]
            obtained_dna_string = "".join(obtained_dna_string)
            # "D" refers to repair by deletion.
            repair_info.append(("D", occur_location, delete_nucleotide, obtained_dna_string))

    return repair_info
